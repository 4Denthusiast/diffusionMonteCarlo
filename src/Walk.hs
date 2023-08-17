{-# LANGUAGE ExistentialQuantification #-}
{-# LANGUAGE TupleSections #-}

module Walk (
    PopulationState(..),
    population,
    totalAmplitude,
    initialPopState,
    step
) where

import Particle
import Configuration
import Ansatz
import Variance

import qualified Data.Array as A
import Data.Complex
import Data.List
import Data.Maybe
import Debug.Trace
import Control.Monad.Random
import Control.Monad.Trans.State
import Math.List.FFT

data PopulationState = forall a. Ansatz a => PopState {walkerSet :: [Walker], energy :: Double, deltaTime :: Double, variance :: Variance, setPoint :: Double, requiredError :: Double, dataSeries :: [Double], ansatz :: a}

population :: PopulationState -> Int
population = length . walkerSet

totalAmplitude :: PopulationState -> Double
totalAmplitude = sum . map amplitude . walkerSet

initialPopState :: Configuration -> Double -> Int -> Double -> String -> PopulationState
initialPopState c dt n e a = case a of
        "" -> incomplete ()
        "c3d" -> incomplete cuspsJastrow3d
        "h3d" -> incomplete hydrogenStyle3d
        _ -> error ("Unrecognised ansatz name: \""++a++"\"")
    where incomplete :: forall x. Ansatz x => x -> PopulationState
          incomplete = PopState (replicate n (Walker c 1 0)) 0 dt [] (fromIntegral n) e []

type M x = RandT StdGen (State PopulationState) x

getEnergy :: M Double
getEnergy = energy <$> lift get

getRequiredError :: M Double
getRequiredError = requiredError <$> lift get

uniformVar :: M Double
uniformVar = liftRandT (pure . random)

normalVar :: M Double
normalVar = do
    d <- uniformVar
    t <- uniformVar
    return $ sqrt (-2*log d) * cos (pi*t)

stepWalkers :: [Walker] -> M [Walker]
stepWalkers = (concat <$>) . mapM (\w -> if localTime w >= 0 then return [w] else stepWalkers =<< splitWalker =<< moveWalker w)

splitWalker :: Walker -> M [Walker]
splitWalker (Walker c a t) =
    if abs a < 0.5 then (\x -> if abs a < x then pure [] else pure [Walker c 1 t]) =<< uniformVar
    else if abs a > 2 then pure [Walker c (a/2) t, Walker c (a/2) t]
    else pure [Walker c a t]

moveWalker :: Walker -> M Walker
moveWalker (Walker c a t) = do
        let Conf ps = c
        dt <- (min (-t) . flip suitableStepSize c) <$> getRequiredError
        PopState{ansatz=an} <- lift get
        c' <- Conf <$> mapM (moveParticle dt) (zip ps (drift an c))
        let v = (potentialEnergy c - aEnergy an c + potentialEnergy c' - aEnergy an c')/2
        e <- getEnergy
        let a' = a * exp (-dt * (v-e))
        return $ Walker c' a' (t+dt)
    where moveParticle :: Double -> ((Position, Particle),[Double]) -> M (Position, Particle)
          moveParticle dt ((r, p), d) = let dt' = dt/particleMass p in (,p) <$> zipWithM (\x dx -> ((x+dt'*dx+) . (sqrt dt'*)) <$> normalVar) r d

loopSteps :: M ()
loopSteps = do
    ps@PopState{walkerSet = ws} <- lift get
    ws' <- stepWalkers ws
    lift $ put ps{walkerSet = ws'}

shiftWalkerTimes :: M ()
shiftWalkerTimes = do
    dt <- deltaTime <$> lift get
    let shift (Walker c a t) = (Walker c a (t-dt))
    lift $ modify (\ps -> ps{walkerSet = map shift $ walkerSet ps})

doubleWalkerSet :: M ()
doubleWalkerSet = lift get >>= (\ps -> lift (put ps{walkerSet = walkerSet ps ++ walkerSet ps}))

trimWalkerSet :: M ()
trimWalkerSet = lift $ do
        ps <- get
        put ps{walkerSet = halve $ walkerSet ps}
    where halve (x:x':xs) = x:halve xs
          halve xs = xs

step :: M Double
step = do
    ps0 <- lift get
    if (population ps0 == 0) then error "Sample population extinct." else return ()
    let pop = totalAmplitude ps0
    loopSteps
    ps <- lift $ get
    let pop' = totalAmplitude ps
    let growthEstimate = (\x -> traceShow (-log x/deltaTime ps) x) $ traceShowId $ exp (-energy ps * deltaTime ps) * pop'/pop
    let variance' = addDataPoint growthEstimate (variance ps)
    let relaxationTime = traceShowId $ deltaTime ps * getTimeAtBound 0.5 variance'
    let (averagedGrowth, randomError) = getMeanAndStddev variance'
    let energy' = -log averagedGrowth / deltaTime ps + log (setPoint ps / pop) / relaxationTime
    lift $ put ps{energy = energy', variance = variance' {-, dataSeries = energyEstimate : dataSeries ps-}}
    shiftWalkerTimes
    if pop' < setPoint ps / 4 then doubleWalkerSet else return ()
    if pop' > setPoint ps * 4 then trimWalkerSet else return ()
    return (randomError / deltaTime ps)


-- For debugging purposes. Potentially also useful for identifying excited states.
showCorrelation :: M ()
showCorrelation = do
        l' <- dataSeries <$> lift get
        let l = map (\x -> x-average l') l'
        let s = intercalate "," $ take 300 $ map (show . realPart) $ ifft $ map (\x -> x*conjugate x/(fromIntegral (length l)^2)) $ fft $ map (:+ 0) l
        trace (show (length l) ++ "\n" ++ s ++ "\n\n") $ return ()
    where averages n [] = []
          averages n xs = average (take n xs) : averages n (drop n xs)
          average xs = sum xs / fromIntegral (length xs)
