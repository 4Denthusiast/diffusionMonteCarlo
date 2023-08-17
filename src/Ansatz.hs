module Ansatz (
    Ansatz(..),
    JastrowAnsatz(..),
    cuspsJastrow3d,
    hydrogenStyle3d
) where

import Particle
import Configuration

import Data.List
import Debug.Trace

{-
For an ansatz Ψ(x0,x1,...), the function drift just gives (∇_xi Ψ)/Ψ for each particle.
The change in position this implies per time step is drift*dt/m.
aEnergy gives Σ_i (∇_xi^2 Ψ)/2mΨ, which should be used as the energy offset.
-}
class Ansatz a where
    aValue :: a -> Configuration -> Double
    drift :: a -> Configuration -> [[Double]]
    aEnergy :: a -> Configuration -> Double

instance Ansatz () where
    aValue () _ = 1
    drift () (Conf ps) = map (\_ -> replicate (confDimension (Conf ps)) 0) ps
    aEnergy () _ = 0

-- Jastrow f has value exp(Σ_p1 Σ_p2 fst f(p1,p2,d(p1,p2))). The second and third components of f's result are the first and second derivatives of (fst . f p1 f2).
data JastrowAnsatz = Jastrow (Particle -> Particle -> Double -> (Double, Double, Double))

{-
∇^2 exp(f)
= ∇. (exp(f) ∇f)
= (∇f exp(f) . ∇f + exp(f) ∇^2f)
= exp(f) (∇f^2 + ∇^2f)

∇^2 f(r)
= ∇. (f'(r) x/r)
= f''(r) (x/r)^2 + f'(r) ∇x/r + f'(r) x.∇(1/r)
= f''(r) + d f'(r)/r - f'(r)/r
= f''(r) + (d-1) f'(r)/r

With a potential of -Z/r in 3D, this gives f'(r) -> -f(0)Z/2 as r->0, and for two dissimilar-spin electrons it's f'(r) -> f(0)/4.
In 4D, if there are no neutrons, it's just f'(r) -> -f(0)kZ/3 instead because of the (d-1) factor and electric constant k. If the neutron number (and thus excess s-charge) N>0, f(r) ~ r^x for small r, where x(x-1) + 3x = kN therefore x=sqrt(1+kN)-1 (there are more solutions but they're divergent). For opposite-spin electrons, it's x=sqrt(1+k/2)-1.
In 2D, there is no cusp. f'(r) -> 0 as r -> 0, though not smoothly.
-}

-- There are a lot of re-calculations involved in this, but for now I'd rather keep it neat than make it efficient. Maybe the compiler will spot some of them.
instance Ansatz JastrowAnsatz where
    aValue (Jastrow f) (Conf ps) = exp $ sum $ concat $ zipWith (\(x0,p0) ps' -> map (\(x1,p1) -> fst3 $ f p0 p1 (dist x0 x1)) ps') ps (tail (tails ps))
    drift (Jastrow f) (Conf ps) = flip map ps (\(x0,p0) -> foldr1 (zipWith (+)) $ flip map ps (\(x1,p1) -> if x0 == x1 then repeat 0.0 else let
            x = zipWith (-) x0 x1
            r = sqrt $ sum $ map (^2) x
            (_,f',_) = f p0 p1 r
        in map (*(f'/r)) x))
    aEnergy (Jastrow f) (Conf ps) = sum $ zipWith (\(x0,p0) dr -> (sum (map (^2) dr) / (2*particleMass p0)) + sum (flip map ps (\(x1,p1) -> if x0 == x1 then 0 else let
                x = zipWith (-) x0 x1
                r = sqrt $ sum $ map (^2) x
                (_,f',f'') = f p0 p1 r
            in (f'' + (fromIntegral (confDimension (Conf ps)) - 1) * f'/r)/(2*particleMass p0)
        ))) ps (drift (Jastrow f) (Conf ps))

fst3 (x,_,_) = x

cuspsJastrow3d :: JastrowAnsatz
cuspsJastrow3d = Jastrow $ \p0 p1 r -> let s = particleCharge p0 * particleCharge p1/(1/particleMass p0 + 1/particleMass p1) in (-s/(r+1), s/(r+1)^2, -2*s/(r+1)^3)

-- Exact for hydrogen-like atoms (in 3D only)
hydrogenStyle3d :: JastrowAnsatz
hydrogenStyle3d = Jastrow $ \p0 p1 r -> let x = particleCharge p0 * particleCharge p1/(1/particleMass p0 + 1/particleMass p1) in (x*r, x, 0)
