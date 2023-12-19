module Particle (
    Particle(..),
    particleMass,
    particleCharge
) where

data Particle = Electron | Nucleus Double Double Double deriving Eq

particleMass :: Particle -> Double
particleMass Electron = 1
particleMass (Nucleus _ _ m) = m

particleCharge :: Particle -> Double
particleCharge Electron = -1
particleCharge (Nucleus n _ _) = n

instance Show Particle where
    show Electron = "e"
    show (Nucleus z a m) = "Z"++show z++"A"++show a++(if m<1/0 then "m"++show m else "")
