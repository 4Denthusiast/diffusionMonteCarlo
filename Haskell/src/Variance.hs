module Variance (
    Variance(..),
    emptyVariance,
    addDataPoint,
    getTimeAtBound,
    layerVariance
) where
-- Since my efforts to find a nice expression for the variance based on the properties of the diffusion algorithm have failed, I have to just make a generic solution based on the observed properties of the data instead.

import Data.List
import Debug.Trace

data Layer = Layer {
    latest :: Maybe (Double,Double),
    lTotalMean :: Double,
    lTotalSquare :: Double,
    lTotalWeight :: Double
} deriving (Show)

data Variance = Variance {
    mean :: Double,
    stdDev :: Double,
    layers :: [Layer]
}

emptyVariance :: Variance
emptyVariance = Variance 0 (1/0) []

addDataPoint :: Double -> Double -> Variance -> Variance
addDataPoint w x v@(Variance _ _ ls) = if isNaN x then v else getMeanAndStddev $ addDataPoint0 x w ls

addDataPoint0 :: Double -> Double -> [Layer] -> [Layer]
addDataPoint0 x w [] = [Layer (Just (x,w)) (x*w) (x*x*w) w]
addDataPoint0 x w (Layer Nothing tm ts tw : ls) = Layer (Just (x,w)) (tm + x*w) (ts + x*x*w) (tw+w) : ls
addDataPoint0 x w (Layer (Just (x',w')) tm ts tw : ls) = Layer Nothing (tm + x*w) (ts + x*x*w) (tw+w) : addDataPoint0 ((x*w+x'*w')/(w+w')) ((w+w')/2) ls

-- An approximation to the greatest time t such that the standard deviation of the total of the data points over a time t is less than the bound given.
getTimeAtBound :: Double -> Variance -> Double
getTimeAtBound _ (Variance _ _ []) = 1
getTimeAtBound b (Variance _ _ (l : ls)) = if layerVariance l < b*b then 2*getTimeAtBound (b/2) (Variance 0 0 ls) else 1

getMeanAndStddev :: [Layer] -> Variance
getMeanAndStddev [] = Variance 0 (1/0) [] -- This is very much based on where this will be used rather than trying to make sense.
getMeanAndStddev ls@(Layer _ tm0 _ tw0 : _) = Variance m (if length ls < 5 then 1/0 else sd) ls
    where m = tm0/tw0
          lastFew = init $ drop (div (length ls) 2) ls
          regressionPoints = snd $ mapAccumL (\l (Layer _ tm ts tw) -> (l/2, (log l, log ((ts-tm*tm/tw)/(tw-1)), max 0 $ tw-1))) (lTotalWeight $ head lastFew) lastFew --This returns a list of triples (scale, variance at that scale, confidence). I'm not sure the calculation for the confidence is correct though.
          totalWeight = sum $ map (\(_,_,c) -> c) regressionPoints
          averageScale = sum (map (\(s,_,c) -> s*c) regressionPoints) / totalWeight
          averageVariance = sum (map (\(_,v,c) -> v*c) regressionPoints) / totalWeight
          gradient = sum (map (\(s,v,c) -> (v-averageVariance)*(s-averageScale)*c) regressionPoints) / sum (map (\(s,_,c) -> (s-averageScale)**2 *c) regressionPoints)
          intercept = averageVariance - gradient*averageScale
          sd = exp (intercept/2)

-- In order to get the right sample variance, tw is used as an approximation for the number of samples, which requires that all of the weights are of order 1.
layerVariance :: Layer -> Double
layerVariance (Layer _ tm ts tw) = (ts-tm*tm/tw)/(max 0 $ tw-1)
