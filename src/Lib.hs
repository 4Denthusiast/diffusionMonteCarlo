{-# LANGUAGE TupleSections #-}

module Lib (
    stepAndTrace,
    initialPopState,
    hydrogen,
    hydride,
    initSystem,
    resetIteration,
) where

import Configuration
import Walk

import Control.Monad.Random
import Control.Monad.Trans.State
import Text.Printf


stepAndTrace :: (StdGen, PopulationState) -> IO (StdGen, PopulationState)
stepAndTrace (r, ps) = putStrLn (showPopulationState ps) >> pure (runState (execRandT step r) ps)

initSystem :: Double -> Int -> Configuration -> IO (StdGen, PopulationState)
initSystem dt n c = (,initialPopState c dt n) <$> getStdGen

showPopulationState :: PopulationState -> String
showPopulationState ps = printf "pop=%d, pop(w)=%.2f, E=%+.5f +- %.1e" (population ps) (totalAmplitude ps) (energySum ps / fromIntegral (iterations ps)) (energyUncertainty ps)

resetIteration :: (StdGen, PopulationState) -> (StdGen, PopulationState)
resetIteration (r, ps) = (r,ps{variance=0,energySum=0,iterations=0})
