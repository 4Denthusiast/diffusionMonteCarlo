{-# LANGUAGE TupleSections #-}

module Lib (
    stepAndTrace,
    initialPopState,
    initSystem,
    resetIteration,
) where

import Configuration
import Walk
import Variance

import Control.Monad.Random
import Control.Monad.Trans.State
import Text.Printf


stepAndTrace :: (StdGen, PopulationState) -> IO (Double, (StdGen, PopulationState))
stepAndTrace (r, ps) = putStrLn (showPopulationState ps) >> pure ((\((x,y),z) -> (x,(y,z))) $ runState (runRandT step r) ps)

-- TODO: merge most of these arguments into one data structure that's passed straight from argument parsing to the PopulationState
initSystem :: Double -> Int -> Double -> Rand StdGen Configuration -> String -> String -> Double -> IO (StdGen, PopulationState)
initSystem dt n e c a m f = do
    (cs,g) <- runRand (replicateM n c) <$> getStdGen
    pure (g,initialPopState cs dt n e a m f)

showPopulationState :: PopulationState -> String
showPopulationState ps = printf (if abs e >= 0.1 then "pop=%d, pop(w)=%.2f, E=%+.5f +- %.1e" else "pop=%d, pop(w)=%.2f, E=%+.5e +- %.1e") (population ps) (totalAmplitude ps) e (error / deltaTime ps)
    where error = stdDev $ variance ps
          e = -log (mean $ variance ps) / deltaTime ps

resetIteration :: (StdGen, PopulationState) -> (StdGen, PopulationState)
resetIteration (r, ps) = (r,ps{variance=emptyVariance, measureVariances=map (const emptyVariance) $ measurementsReqd ps, dataSeries=[], walkerSet=map resetWalker (walkerSet ps)})
    where resetWalker w = w{measurementValues = map (const 0) $ measurementsReqd ps}
