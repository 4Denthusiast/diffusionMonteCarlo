module Particle (
    Particle(..),
    particleMass,
    particleCharge
) where

data Particle = Electron | Nucleus Int deriving Eq

particleMass :: Particle -> Double
particleMass Electron = 1
particleMass (Nucleus _) = 1/0

particleCharge :: Particle -> Double
particleCharge Electron = -1
particleCharge (Nucleus n) = fromIntegral n
