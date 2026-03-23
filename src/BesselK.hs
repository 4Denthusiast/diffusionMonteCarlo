module BesselK (
  besselK1x
) where

import Data.Array

-- Following https://arxiv.org/abs/1209.1547
besselK :: Double -> Double -> Double
besselK ν 0 = 1/0
besselK ν z = tryH 1 (withH 2)
    where tryH h prev = let new = withH h in if abs (new - prev) <= 1e-11 * prev then new else tryH (h/2) new
          withH h = (*h) $ sum $ (exp (-z) / 2 :) $ takeWhile (> exp (-z) * 1e-20) $ map (\n -> cosh (ν*n*h) * exp (-z * cosh (n*h))) [1..]

data Quartic = Quartic Double Double Double Double Double

evalQuartic :: Quartic -> Double -> Double
evalQuartic (Quartic a b c d e) x = a + x*(b + x*(c + x*(d + x*e)))

-- taylor series of x^2 K_1(x^2)
besselK1x2TaylorSeries :: Double -> Quartic
besselK1x2TaylorSeries x = Quartic
    (x ^ 2 * k1)
    (-2 * x^3 * k0)
    ((4 * x^4 * k1 - 6 * x^2 * k0)/2)
    ((20 * x^3 * k1 - (8 * x^5 + 12 * x) * k0)/6)
    (((16 * x^6 + 44*x^2) * k1 - (80 * x^4 + 12) * k0)/24)
  where k0 = besselK 0 (x ^ 2)
        k1 = besselK 1 (x ^ 2)

gridFineness = 1e-3
gridMax = 6

-- I thought using an unboxed array of doubles might be quicker here, but when I benchmarked it in GHCI, the array of `Quartic`s turned out to actually be a little faster.
memoized :: Array Int Quartic
memoized = listArray (0, limit) [besselK1x2TaylorSeries ((fromIntegral x + 0.5)*gridFineness) | x <- [0..limit]]
  where limit = ceiling (gridMax/gridFineness) - 1

-- x K1(x), using the memoized table of taylor expansions.
besselK1x :: Double -> Double
besselK1x x | x >= gridMax * gridMax = 0
besselK1x x | otherwise    = evalQuartic (memoized ! i) dx
  where i = round (sqrt x/gridFineness - 0.5)
        dx = sqrt x - (fromIntegral i + 0.5) * gridFineness
