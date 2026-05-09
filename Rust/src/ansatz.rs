use dyn_clone::DynClone;

use crate::particle::*;
use crate::configuration::*;
use crate::walk::*;

/*
For an ansatz Ψ(x0,x1,...), the function drift just gives (∇_xi Ψ)/Ψ for each particle.
The change in position this implies per time step is drift*dt/m.
aEnergy gives Σ_i (∇_xi^2 Ψ)/2mΨ, which should be used as the energy offset.
*/
pub trait Ansatz : DynClone + Send + Sync + 'static {
  fn value(&self, configuration: &[Position], context: &Context) -> f64;
  fn drift(&self, configuration: &[Position], context: &Context, index: usize) -> Position; // index must be < configuration.len()
  fn energy(&self, configuration: &[Position], context: &Context) -> f64;
  
  // This has to be an Ansatz method in order for dynamic dispatch to work.
  fn step_walkers_dynamic(&self, pop : PopulationState, ctx : &Context) -> PopulationState;
}

dyn_clone::clone_trait_object!(Ansatz);

// This separate trait exists to allow Ansatz to have a default implementation for step_walkers_dynamic, which requires Self : Sized, without making Ansatz a sub-trait of Sized (which would make it not dyn compatible).
trait AnsatzBase : DynClone + Send + Sync + 'static {
  fn value(&self, configuration: &[Position], context: &Context) -> f64;
  fn drift(&self, configuration: &[Position], context: &Context, index: usize) -> Position;
  fn energy(&self, configuration: &[Position], context: &Context) -> f64;
}

impl<T> Ansatz for T where T : AnsatzBase + Sized{
  fn value(&self, configuration: &[Position], context: &Context) -> f64 {self.value(configuration, context)}
  fn drift(&self, configuration: &[Position], context: &Context, index: usize) -> Position {self.drift(configuration, context, index)}
  fn energy(&self, configuration: &[Position], context: &Context) -> f64 {self.energy(configuration, context)}
  fn step_walkers_dynamic(&self, pop : PopulationState, ctx : &Context) -> PopulationState {
    step_walkers(pop, ctx, self)
  }
}

#[derive(Clone)]
pub struct TrivialAnsatz;

impl AnsatzBase for TrivialAnsatz {
  fn value(&self, _: &[Position], _: &Context) -> f64 {1.0}
  fn drift(&self, _: &[Position], _: &Context, _: usize) -> Position {Position::default()}
  fn energy(&self, _: &[Position], _: &Context) -> f64 {0.0}
}

// Just something very easy to calculate the derivatives of, to make sure the ansatz system in general is working.
#[derive(Clone)]
pub struct SineTestAnsatz;

impl AnsatzBase for SineTestAnsatz {
  fn value(&self, configuration: &[Position], _: &Context) -> f64 {configuration[1].x.sin() + 2.0}
  fn drift(&self, configuration: &[Position], _: &Context, index: usize) -> Position {
    if index != 1 {
      return Position::default()
    }
    let x = configuration[1].x;
    Position {
      x : x.cos() / (x.sin() + 2.0),
      y : 0.0, z : 0.0, w : 0.0
    }
  }
  fn energy(&self, configuration: &[Position], ctx: &Context) -> f64 {
    let x = configuration[1].x;
    - x.sin() / (2.0 * particle_mass(&ctx.molecule[1]) * (x.sin() + 2.0))
  }
}

// exp (Σ_p1≠p2 dist(p1,p2) q(p1) q(p2))
// This gives correct cusps for infinitely massive nuclei, and the correct limit for cases with 1 distant electron, but incorrect cusps for electron-electron interactions. Still, hopefully it's a good enough approximation to be useful.
#[derive(Clone)]
pub struct Cusps3DSimple;

impl AnsatzBase for Cusps3DSimple {
  fn value(&self, configuration: &[Position], ctx: &Context) -> f64 {
    let mut log_value = 0.0;
    for i in 0..configuration.len() {
      for j in 0..i {
        let r = configuration[i].dist(&configuration[j]);
        log_value += r * particle_charge(&ctx.molecule[i]) * particle_charge(&ctx.molecule[j]);
      }
    }
    log_value.exp()
  }
  fn drift(&self, configuration: &[Position], ctx: &Context, index: usize) -> Position {
    let mut result = Position::default();
    for i in 0..configuration.len() {
      if i != index {
        let pos = configuration[index] - configuration[i];
        let r = pos.dist(&Position::default());
        result = result + pos * (particle_charge(&ctx.molecule[i]) / r);
      }
    }
    result * particle_charge(&ctx.molecule[index])
  }
  fn energy(&self, configuration: &[Position], ctx: &Context) -> f64 {
    let mut result = 0.0; // The total laplacian of ln(value)
    for i in 0..configuration.len() {
      if !particle_mass(&ctx.molecule[i]).is_finite() {
        continue;
      }
      let mut particle_gradient = Position::default(); // the gradient of ln(value)
      let mut particle_laplacian = 0.0; // the laplacian of ln(value)
      for j in 0..configuration.len() { // TODO: refactor to only use the j < i part.
        if i != j {
          let pos = configuration[i] - configuration[j];
          let r = pos.dist(&Position::default());
          particle_gradient = particle_gradient + pos * (particle_charge(&ctx.molecule[j]) / r);
          particle_laplacian += (ctx.dimension - 1) as f64 * particle_charge(&ctx.molecule[j]) / r; // This ansatz is only designed to be good in 3D, but I'm keeping it general just in case.
        }
      }
      particle_gradient = particle_gradient * particle_charge(&ctx.molecule[i]);
      particle_laplacian *= particle_charge(&ctx.molecule[i]);
      result += (particle_gradient.dist2(&Position::default()) + particle_laplacian) / (2.0 * particle_mass(&ctx.molecule[i]));
    }
    result
  }
}
/*
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

-- Optimised for getting good pictures rather than good energy values.
data PictureAnsatz = PictureAnsatz Double

instance Ansatz PictureAnsatz where
    aValue (PictureAnsatz r0) (Conf ps) = product $ map (\(p,_) -> recip $ sqrt $ (r0+) $ sum $ map (^2) p) ps
    drift (PictureAnsatz r0) (Conf ps) = map (\(p,_) -> map (\x -> -x/(r0 + sum (map (^2) p))) p) ps
    aEnergy (PictureAnsatz r0) (Conf ps) = sum $ map (\(p,t) -> let
            d = fromIntegral $ length p
            r = sqrt $ sum $ map (^2) p
        in (3*r*r/(r0+r*r)-d)/(r0+r*r)/(2*particleMass t)) ps
--d/dr 1/sqrt(r0+r^2) = -1/2 (r0+r^2)^-3/2 2r = -r/(r0+r^2)^3/2
--d/dr p(r0,r) / p(r0,r) = -r/(r0+r^2)
--d^2/dr^2 1/sqrt(r0+r^2) = d/dr (-r/(r0+r^2)^3/2) = -1/(r0+rr)^3/2 + 3rr/(r0+rr)^5/2
-- ∇^2 p(r0,r)/p(r0,r) = -1/(r0+rr) + 3rr/(r0+rr)^2 - (d-1)/(r0+rr)
*/
