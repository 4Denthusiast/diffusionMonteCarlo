{-# LANGUAGE TupleSections #-}

module Walk (
    PopulationState(..),
    initialPopState,
    step
) where

import Particle
import Configuration

import Control.Monad.Random
import Control.Monad.Trans.State

data PopulationState = PopState {energy :: Double, popVariance :: Double, variance :: Double, deltaTime :: Double, setPoint :: Double, energySum :: Double, iterations :: Int, energyUncertainty :: Double}

initialPopState :: Double -> Int -> PopulationState
initialPopState dt n = PopState 0 0 0 dt (fromIntegral n) 0 0 0

type M x = RandT StdGen (State PopulationState) x

getEnergy :: M Double
getEnergy = energy <$> lift get

getDeltaTime :: M Double
getDeltaTime = deltaTime <$> lift get

incVariance :: M ()
incVariance = lift $ modify (\ps -> ps{popVariance = popVariance ps + 1})

uniformVar :: M Double
uniformVar = liftRandT (pure . random)

normalVar :: M Double
normalVar = do
    d <- uniformVar
    t <- uniformVar
    return $ sqrt (-2*log d) * cos (pi*t)

stepWalker :: Walker -> M [Walker]
stepWalker w = splitWalker =<< moveWalker w

-- TODO
splitWalker :: Walker -> M [Walker]
splitWalker (Walker c a) =
    if abs a < 1 then (\x -> if abs a < x then incVariance >> pure [] else pure [Walker c 1]) =<< uniformVar
    else if abs a > 2 then incVariance >> pure [Walker c (a/2), Walker c (a/2)]
    else pure [Walker c a]

moveWalker :: Walker -> M Walker
moveWalker (Walker c a) = do
        let Conf ps = c
        dt <- getDeltaTime
        c' <- Conf <$> mapM (moveParticle dt) ps
        let v = (potentialEnergy c + potentialEnergy c')/2
        e <- getEnergy
        let a' = a * exp (-dt * (v-e))
        return (Walker c' a')
    where moveParticle :: Double -> (Position, Particle) -> M (Position, Particle)
          moveParticle dt (r, p) = (,p) <$> mapM (\x -> ((x+) . (sqrt (dt/particleMass p)*)) <$> normalVar) r

step :: [Walker] -> M [Walker]
step ws = do
    let pop = (fromIntegral $ length ws) :: Double
    ws' <- concat <$> mapM stepWalker ws
    let pop' = (fromIntegral $ length ws') :: Double
    ps <- lift $ get
    let energyEstimate = energy ps + log (pop/pop') / deltaTime ps
    let relaxationTime = deltaTime ps * sqrt (pop * pop / (popVariance ps + 1))
    let energy' = (energySum ps + energyEstimate) / fromIntegral (iterations ps + 1) + log (setPoint ps / pop) * 1 / relaxationTime
    lift $ put ps{energy = energy', energySum = energySum ps + energyEstimate, popVariance = 0, variance = variance ps + popVariance ps / (pop * deltaTime ps)^2, iterations = iterations ps + 1, energyUncertainty = 1/relaxationTime}
    let ws'' = take (floor $ setPoint ps * 4) ws' -- Stop the population exploding while the energy equilibriates.
    let ws''' = if (pop < setPoint ps/4) then take (floor $ setPoint ps/4) (cycle ws'') else ws'' -- Likewise, stop it dying out.
    return ws'''
