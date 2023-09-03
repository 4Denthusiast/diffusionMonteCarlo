module Configuration (
    Position,
    Configuration(..),
    Walker(..),
    confDimension,
    potentialEnergy,
    suitableStepSize,
    dist,
    dist2,
    createAtom
) where

import Particle

import Data.List
import Debug.Trace


type Position = [Double]

data Configuration = Conf [(Position, Particle)] deriving Show

data Walker = Walker {configuration :: Configuration, amplitude :: Double, localTime :: Double} deriving Show

confDimension :: Configuration -> Int
confDimension (Conf ps) = length $ fst $ head ps

potentialEnergy :: Configuration -> Double
potentialEnergy c@(Conf ps) = sum $ concat $ zipWith (\p t -> map (pairEnergy p) t) ps (tail $ tails ps)
    where pairEnergy (r0,p0) (r1,p1) = coulumb (dist r0 r1) * particleCharge p0 * particleCharge p1 * if d < 4 || (p0 == Electron) == (p1 == Electron) then 1 else (1 - 1.8 * exp (-dist r0 r1))
          d = confDimension c
          coulumb = case d of
              2 -> negate . log
              3 -> recip
              4 -> (coulumb4D*) . (^^(-2)) -- In this case, it can't be normalised to 1.
          coulumb4D = 0.9

-- Decrease the step-size when the potential energy is locally highly variable (i.e. when particles are nearby), and adapt to the desired energy error.
suitableStepSize :: Double -> Configuration -> Double
suitableStepSize de c@(Conf ps) = min ordinaryBound radiusBound
    where ordinaryBound = (de/(sum $ concat $ zipWith (\p t -> map (dv2 p) t) ps (tail $ tails ps)))**(1/3) -- t^3 < E/|V|^2
          dv2 (r0,p0) (r1,p1) = dist2 r0 r1 ^^ (1-confDimension c)
          radiusBound = (0.07*) $ minimum $ concat $ zipWith (\p t -> map (rb p) t) ps (tail $ tails ps) -- the distance that may be travelled is less than the distance between any two particles
          rb (r0,p0) (r1,p1) = dist2 r0 r1

dist2 :: Position -> Position -> Double
dist2 r0 r1 = sum $ map (^2) $ zipWith (-) r0 r1

dist :: Position -> Position -> Double
dist r0 r1 = sqrt $ dist2 r0 r1

-- Fill with 0s up to the specified dimension.
pad :: Int -> Position -> Position
pad d r = take d $ r ++ repeat 0

createAtom :: Int -> Int -> Int -> Configuration
createAtom d z q = Conf $ (pad d [0], Nucleus z) : map (\n -> (pad d [fromIntegral n], Electron)) [1..(z-q)]
