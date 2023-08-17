{-# LANGUAGE TupleSections #-}

module Lib (
    stepAndTrace,
    initialPopState,
    createAtom,
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

initSystem :: Double -> Int -> Double -> Configuration -> String -> IO (StdGen, PopulationState)
initSystem dt n e c a = (,initialPopState c dt n e a) <$> getStdGen

showPopulationState :: PopulationState -> String
showPopulationState ps = printf "pop=%d, pop(w)=%.2f, E=%+.5f +- %.1e" (population ps) (totalAmplitude ps) (-log growth / deltaTime ps) (error / deltaTime ps)
    where (growth, error) = getMeanAndStddev $ variance ps

resetIteration :: (StdGen, PopulationState) -> (StdGen, PopulationState)
resetIteration (r, ps) = (r,ps{variance=[],dataSeries=[]})
