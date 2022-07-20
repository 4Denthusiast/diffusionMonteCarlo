{-# LANGUAGE TupleSections #-}

module Walk (
    PopulationState(..),
    population,
    totalAmplitude,
    energyUncertainty,
    initialPopState,
    step
) where

import Particle
import Configuration
import qualified PriorityQueue as PQ

import Data.Maybe
import Debug.Trace
import Control.Monad.Random
import Control.Monad.Trans.State

type WalkerSet = PQ.PriorityQueue Double Walker

data PopulationState = PopState {walkerSet :: WalkerSet, energy :: Double, popVariance :: Double, variance :: Double, deltaTime :: Double, setPoint :: Double, energySum :: Double, iterations :: Int}

energyUncertainty :: PopulationState -> Double
energyUncertainty ps = if variance ps == 0 then 1 else sqrt (variance ps) / fromIntegral (iterations ps)

population :: PopulationState -> Int
population = PQ.size . walkerSet

totalAmplitude :: PopulationState -> Double
totalAmplitude = sum . map (amplitude . snd) . PQ.toAscList . walkerSet

initialPopState :: Configuration -> Double -> Int -> PopulationState
initialPopState c dt n = PopState (PQ.fromList (replicate n (0,Walker c 1 0))) 0 0 0 dt (fromIntegral n) 0 0

type M x = RandT StdGen (State PopulationState) x

getEnergy :: M Double
getEnergy = energy <$> lift get

getEnergyUncertainty :: M Double
getEnergyUncertainty = energyUncertainty <$> lift get

incVariance :: M ()
incVariance = lift $ modify (\ps -> ps{popVariance = popVariance ps + 1})

uniformVar :: M Double
uniformVar = liftRandT (pure . random)

normalVar :: M Double
normalVar = do
    d <- uniformVar
    t <- uniformVar
    return $ sqrt (-2*log d) * cos (pi*t)

stepWalker :: M ()
stepWalker = pushWalkers =<< splitWalker =<< moveWalker =<< popWalker

-- TODO
splitWalker :: Walker -> M [Walker]
splitWalker (Walker c a t) =
    if abs a < 1 then (\x -> if abs a < x then incVariance >> pure [] else pure [Walker c 1 t]) =<< uniformVar
    else if abs a > 2 then incVariance >> pure [Walker c (a/2) t, Walker c (a/2) t]
    else pure [Walker c a t]

moveWalker :: Walker -> M Walker
moveWalker (Walker c a t) = do
        let Conf ps = c
        dt <- (min (-t) . flip suitableStepSize c . const 0.1) <$> getEnergyUncertainty
        c' <- Conf <$> mapM (moveParticle dt) ps
        let v = (potentialEnergy c + potentialEnergy c')/2
        e <- getEnergy
        let a' = a * exp (-dt * (v-e))
        return (Walker c' a' (t+dt))
    where moveParticle :: Double -> (Position, Particle) -> M (Position, Particle)
          moveParticle dt (r, p) = (,p) <$> mapM (\x -> ((x+) . (sqrt (dt/particleMass p)*)) <$> normalVar) r

pushWalkers :: [Walker] -> M ()
pushWalkers ws = lift $ modify (\ps -> ps{walkerSet = foldr (\w -> PQ.insert (localTime w) w) (walkerSet ps) ws})

popWalker :: M Walker
popWalker = lift $ do
    ps <- get
    let Just (w, ws') = PQ.pop $ walkerSet ps
    put ps{walkerSet = ws'}
    return w

nextTime :: M Double
nextTime = fromMaybe (1/0) <$> PQ.next <$> walkerSet <$> lift get

loopSteps :: M ()
loopSteps = do
    nt <- nextTime
    if nt < 0 then stepWalker >> loopSteps
    else return ()

shiftWalkerTimes :: M ()
shiftWalkerTimes = do
    dt <- deltaTime <$> lift get
    let shift (_, Walker c a t) = (t-dt, Walker c a (t-dt))
    lift $ modify (\ps -> ps{walkerSet = PQ.fromList $ map shift $ PQ.toAscList $ walkerSet ps})

doubleWalkerSet :: M ()
doubleWalkerSet = lift get >>= (\ps -> lift (put ps{walkerSet = PQ.union (walkerSet ps) (walkerSet ps)}))

trimWalkerSet :: M ()
trimWalkerSet = lift $ do
        ps <- get
        put ps{walkerSet = PQ.fromList $ halve $ PQ.toAscList $ walkerSet ps}
    where halve (x:x':xs) = x:halve xs
          halve xs = xs

step :: M ()
step = do
    ps0 <- lift get
    if (population ps0 == 0) then error "Sample population extinct." else return ()
    let pop = totalAmplitude ps0
    loopSteps
    ps <- lift $ get
    let pop' = totalAmplitude ps
    let energyEstimate = traceShowId $ energy ps + log (pop/pop') / deltaTime ps
    let relaxationTime = deltaTime ps * sqrt (pop * pop / (popVariance ps + 1))
    let energy' = (energySum ps + energyEstimate) / fromIntegral (iterations ps + 1) + log (setPoint ps / pop) * 1 / relaxationTime
    lift $ put ps{energy = energy', energySum = energySum ps + energyEstimate, popVariance = 0, variance = variance ps + popVariance ps / (pop * deltaTime ps)^2, iterations = iterations ps + 1}
    shiftWalkerTimes
    if pop' < setPoint ps / 4 then doubleWalkerSet else return ()
    if pop' > setPoint ps * 4 then trimWalkerSet else return ()
