module Variance (
    Variance,
    addDataPoint,
    getTimeAtBound,
    getMeanAndStddev,
    layerVariance
) where
-- Since my efforts to find a nice expression for the variance based on the properties of the diffusion algorithm have failed, I have to just make a generic solution based on the observed properties of the data instead.

import Data.List
import Debug.Trace

data Layer = Layer {
    latest :: Maybe Double,
    lTotalMean :: Double,
    lTotalSquare :: Double,
    lTotalIncluded :: Double -- This will be an integer, but I'm always using it in a floating-point context.
} deriving (Show)

type Variance = [Layer]

addDataPoint :: Double -> Variance -> Variance
addDataPoint x [] = [Layer (Just x) x (x*x) 1]
addDataPoint x (Layer Nothing tm ts ti : ls) = Layer (Just x) (tm + x) (ts + x*x) (ti+1) : ls
addDataPoint x (Layer (Just y) tm ts ti : ls) = Layer Nothing (tm + x) (ts + x*x) (ti+1) : addDataPoint ((x+y)/2) ls

-- An approximation to the greatest time t such that the standard deviation of the total of the data points over a time t is less than the bound given.
getTimeAtBound :: Double -> Variance -> Double
getTimeAtBound _ [] = 1
getTimeAtBound b (l : ls) = if layerVariance l < b*b then 2*getTimeAtBound (b/2) ls else 1

getMeanAndStddev :: Variance -> (Double, Double)
getMeanAndStddev [] = (0,1/0) -- This is very much based on where this will be used rather than trying to make sense.
getMeanAndStddev ls@(Layer _ tm0 _ ti0 : _) = (mean, if length ls < 5 then 1/0 else stddev)
    where mean = tm0/ti0
          lastFew = init $ drop (div (length ls) 2) ls
          regressionPoints = snd $ mapAccumL (\l (Layer _ tm ts ti) -> (l/2, (log l, log ((ts-tm*tm/ti)/(ti-1)), ti-1))) (lTotalIncluded $ head lastFew) lastFew --This returns a list of triples (scale, variance at that scale, confidence). I'm not sure the calculation for the confidence is correct though.
          totalWeight = sum $ map (\(_,_,c) -> c) regressionPoints
          averageScale = sum (map (\(s,_,c) -> s*c) regressionPoints) / totalWeight
          averageVariance = sum (map (\(_,v,c) -> v*c) regressionPoints) / totalWeight
          gradient = sum (map (\(s,v,c) -> (v-averageVariance)*(s-averageScale)*c) regressionPoints) / sum (map (\(s,_,c) -> (s-averageScale)**2 *c) regressionPoints)
          intercept = averageVariance - gradient*averageScale
          stddev = exp (intercept/2)

layerVariance :: Layer -> Double
layerVariance (Layer _ tm ts ti) = (ts-tm*tm/ti)/(ti-1)
