{-# LANGUAGE ExistentialQuantification #-}
{-# LANGUAGE TupleSections #-}

module Walk (
    PopulationState(..),
    population,
    totalAmplitude,
    initialPopState,
    step,
    normalVar
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

data PopulationState = forall a. Ansatz a => PopState {walkerSet :: [Walker], energy :: Double, deltaTime :: Double, totalTime :: Double, variance :: Variance, setPoint :: Double, requiredError :: Double, dataSeries :: [Double], measurementsReqd :: [Measurement], ansatz :: a}

population :: PopulationState -> Int
population = length . walkerSet

totalAmplitude :: PopulationState -> Double
totalAmplitude = sum . map amplitude . walkerSet

initialPopState :: [Configuration] -> Double -> Int -> Double -> String -> String -> PopulationState
initialPopState cs dt n e a m = case a of
        "" -> incomplete ()
        "c3d" -> incomplete cuspsJastrow3d
        "h3d" -> incomplete hydrogenStyle3d
        "p" -> incomplete (PictureAnsatz 1)
        _ -> error ("Unrecognised ansatz name: \""++a++"\"")
    where incomplete :: forall x. Ansatz x => x -> PopulationState
          incomplete = PopState ws 0 dt 0 emptyVariance (fromIntegral n) e [] mrs
          ws = map (\c -> Walker c 1 0 (map (const 0) mrs)) cs
          mrs = map (\l -> case l of
                  'r' -> MDistance
                  'd' -> MDipole
                  '1' -> MOne
              ) m

type R m x = RandT StdGen m x
type M x = R (State PopulationState) x

getEnergy :: M Double
getEnergy = energy <$> lift get

getRequiredError :: M Double
getRequiredError = requiredError <$> lift get

uniformVar :: (Monad m) => R m Double
uniformVar = liftRandT (pure . random)

normalVar :: (Monad m) => R m Double
normalVar = do
    d <- uniformVar
    t <- uniformVar
    return $ sqrt (-2*log d) * cos (pi*t)

stepWalkers :: [Measurement] -> [Walker] -> M [Walker]
stepWalkers mrs = (concat <$>) . mapM (\w -> if localTime w >= 0 then return [w] else stepWalkers mrs =<< splitWalker =<< moveWalker mrs w)

splitWalker :: Walker -> M [Walker]
splitWalker (Walker c a t ms) =
    if abs a < 0.5 then (\x -> if abs a < x then pure [] else pure [Walker c 1 t ms]) =<< uniformVar
    else if abs a > 2 then pure [Walker c (a/2) t ms, Walker c (a/2) t ms]
    else pure [Walker c a t ms]

moveWalker :: [Measurement] -> Walker -> M Walker
moveWalker mrs (Walker c a t ms) = do
        let Conf ps = c
        dt <- (min (-t) . flip suitableStepSize c) <$> getRequiredError
        PopState{ansatz=an} <- lift get
        c' <- Conf <$> mapM (moveParticle dt) (zip ps (drift an c))
        let v = (potentialEnergy c - aEnergy an c + potentialEnergy c' - aEnergy an c')/2
        e <- getEnergy
        let a' = a * exp (-dt * (v-e))
        let ms' = strictenList $ zipWith (+) ms $ map ((*dt) . measure c') mrs
        return $! Walker c' a' (t+dt) $! ms'
    where moveParticle :: Double -> ((Position, Particle),[Double]) -> M (Position, Particle)
          moveParticle dt ((r, p), d) = let dt' = dt/particleMass p in (,p) <$> zipWithM (\x dx -> ((x+dt'*dx+) . (sqrt dt'*)) <$> normalVar) r d

strictenList :: [a] -> [a]
strictenList [] = []
strictenList (x:xs) = seq x $ seq xs' $ x:xs'
    where xs' = strictenList xs

loopSteps :: M ()
loopSteps = do
    ps@PopState{walkerSet = ws, measurementsReqd = mrs} <- lift get
    ws' <- stepWalkers mrs ws
    lift $ put ps{walkerSet = ws'}

shiftWalkerTimes :: M ()
shiftWalkerTimes = do
    dt <- deltaTime <$> lift get
    let shift (Walker c a t ms) = (Walker c a (t-dt) ms)
    lift $ modify (\ps -> ps{walkerSet = map shift $ walkerSet ps})

doubleWalkerSet :: M ()
doubleWalkerSet = lift get >>= (\ps -> lift (put ps{walkerSet = walkerSet ps ++ walkerSet ps}))

trimWalkerSet :: M ()
trimWalkerSet = lift $ do
        ps <- get
        put ps{walkerSet = halve $ walkerSet ps}
    where halve (x:x':xs) = x:halve xs
          halve xs = xs

measurementTotals :: PopulationState -> [Double]
measurementTotals (PopState{walkerSet = ws}) = foldr1 (zipWith (+)) $ map (\(Walker{amplitude = a, measurementValues = ms}) -> map (a*) ms) ws

step :: M Double
step = do
    ps0 <- lift get
    if (population ps0 == 0) then error "Sample population extinct." else return ()
    let pop = totalAmplitude ps0
    loopSteps
    ps <- lift $ get
    let pop' = totalAmplitude ps
    let growthEstimate = (\x -> traceShow (-log x/deltaTime ps) x) $ traceShowId $ exp (-energy ps * deltaTime ps) * pop'/pop
    let variance' = addDataPoint growthEstimate (pop'/setPoint ps) (variance ps) -- If I want to reduce population-control bias further by keeping a separate estimate of what the actual population should be, that should be used in the weight parameter here, but it doesn't seem like population control is actually a major issue.
    let relaxationTime = traceShowId $ deltaTime ps * getTimeAtBound 0.5 variance'
    let energy' = -log (mean variance') / deltaTime ps + log (setPoint ps / pop') / relaxationTime
    lift $ put ps{energy = energy', variance = variance', totalTime = totalTime ps + deltaTime ps {-, dataSeries = energyEstimate : dataSeries ps-}}
    shiftWalkerTimes
    if pop' < setPoint ps / 4 then doubleWalkerSet else return ()
    if pop' > setPoint ps * 4 then trimWalkerSet else return ()
    randomSample <- (walkerSet ps !!) <$> liftRandT (pure . randomR (0,length (walkerSet ps) - 1))
    trace (showConfig $ configuration $ randomSample) $ return ()
    trace (unlines $ zipWith (\mr mv -> show mr ++ ": " ++ show (mv / totalTime ps / pop')) (measurementsReqd ps) (measurementTotals ps)) $ return ()
    return (stdDev variance' / deltaTime ps)


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
