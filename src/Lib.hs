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


stepAndTrace :: (([Walker], StdGen), PopulationState) -> IO (([Walker], StdGen), PopulationState)
stepAndTrace ((ws, r), ps) = putStrLn (showPopulationState ws ps) >> pure (runState (runRandT (step ws) r) ps)

initSystem :: Double -> Int -> Configuration -> IO (([Walker], StdGen), PopulationState)
initSystem dt n c = (,initialPopState dt n) <$> (replicate n (Walker c 1),) <$> getStdGen

showPopulationState :: [Walker] -> PopulationState -> String
showPopulationState ws ps = printf "pop=%d, E=%+.5f +- %.1e" (length ws) (energySum ps / fromIntegral (iterations ps)) (sqrt (variance ps) / fromIntegral (iterations ps))

resetIteration :: (([Walker], StdGen), PopulationState) -> (([Walker], StdGen), PopulationState)
resetIteration ((ws, r), ps) = ((ws,r), initialPopState (deltaTime ps) (floor $ setPoint ps))
