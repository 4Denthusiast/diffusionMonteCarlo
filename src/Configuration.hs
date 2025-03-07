module Configuration (
    Position,
    Configuration(..),
    Measurement(..),
    Walker(..),
    confDimension,
    potentialEnergy,
    suitableStepSize,
    measure,
    dist,
    dist2,
    showConfig
) where

import Particle

import Data.List
import Debug.Trace
import Text.Printf


type Position = [Double]

data Configuration = Conf [(Position, Particle)] deriving Show

data Measurement = MDistance | MDipole | MPotential | MOne

instance Show Measurement where
    show MDistance = "Particle separation"
    show MDipole = "Dipole moment"
    show MPotential = "Potential energy"
    show MOne = "Constant one"

data Walker = Walker {configuration :: Configuration, amplitude :: Double, localTime :: Double, measurementValues :: [Double]} deriving Show

confDimension :: Configuration -> Int
confDimension (Conf ps) = length $ fst $ head ps

potentialEnergy :: Configuration -> Double
potentialEnergy c@(Conf ps) = sum $ concat $ zipWith (\p t -> map (pairEnergy p) t) ps (tail $ tails ps)
    where pairEnergy (r0,p0) (r1,p1) = coulumb (dist r0 r1) * case (p0,p1) of
              (Electron,Electron) -> 1
              (Nucleus q q' _,Electron) -> q'*exp(-dist r0 r1)-q
              (Electron,Nucleus q q' _) -> q'*exp(-dist r0 r1)-q
              (Nucleus q1 q'1 _,Nucleus q2 q'2 _) -> q1*q2 - sChargeRatio*q'1*q'2*exp(-dist r0 r1)
          d = confDimension c
          coulumb = case d of
              2 -> negate . log
              3 -> recip
              4 -> (coulumb4D*) . (^^(-2)) -- In this case, it can't be normalised to 1.
          coulumb4D = 0.9
          sChargeRatio = 0 --The ratio of the s charges of a nucleon and an electron. Realistically this should be some high value, but it doesn't matter all that much.

measure :: Configuration -> Measurement -> Double
measure (Conf ((r0,p0):(r1,p1):_)) MDistance = dist r0 r1
measure (Conf ps) MDipole = sum $ map pDipole ps
    where pDipole (r,p) = particleCharge p * (head r - centre)
          centre = if length fixedParticles > 0 then (sum $ map (head . fst) fixedParticles)/fromIntegral (length fixedParticles) else (sum $ map (\(x:_,p) -> x * particleMass p) ps) / (sum $ map (particleMass . snd) ps)
          fixedParticles = filter ((== 1/0) . particleMass . snd) ps
measure x MPotential = potentialEnergy x
measure _ MOne = 1

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

-- Sort of a lossy prettyprint thing. Shows the distance between each pair of particles.
showConfig :: Configuration -> String
showConfig (Conf ps) = show (map snd ps) ++ show (potentialEnergy $ Conf ps) ++ concatMap distsFrom ps
    where distsFrom (x,_) = "\n[" ++ intercalate "," (map (distFrom x) ps) ++ "]"
          distFrom x (y,_) = printf "%.2e" $ dist x y
