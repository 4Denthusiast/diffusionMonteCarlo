use std::fmt::{self, Display};
use std::iter::zip;

use crate::particle::*;
use crate::bessel_k::*;
use crate::walk::Context;

#[derive(Copy, Clone)]
pub struct Position {pub x: f64, pub y: f64, pub z: f64, pub w: f64}

impl Position {
  fn dist2(&self, other: &Position) -> f64 {
    (self.x - other.x) * (self.x - other.x) + 
    (self.y - other.y) * (self.y - other.y) + 
    (self.z - other.z) * (self.z - other.z) + 
    (self.w - other.w) * (self.w - other.w)
  }
  
  fn dist(&self, other: &Position) -> f64 {
    self.dist2(other).sqrt()
  }
}

impl Default for Position {
  fn default() -> Self {
    Position{x:0.0, y:0.0, z:0.0, w:0.0}
  }
}

pub enum Measurement {
  Distance,
  Dipole,
  Potential,
  One,
}
pub const MEASUREMENT_COUNT : usize = 4; // mem::variant_count is not const.

impl Display for Measurement {
  fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    match &self {
      Measurement::Distance => write!(f, "Particle separation"),
      Measurement::Dipole => write!(f, "Dipole moment"),
      Measurement::Potential => write!(f, "Potential energy"),
      Measurement::One => write!(f, "Constant one"),
    }
  }
}

pub fn potential_energy(configuration: &[Position], context: &Context) -> f64 {
  const COULOMB_4D : f64 = 0.9;
  const S_CHARGE_RATIO : f64 = 0.0; // The ratio of the s charges of a nucleon and an electron. Realistically this should be some high value, but it doesn't matter all that much.
  assert!(configuration.len() == context.molecule.len());
  let mut energy = 0.0;
  for i in 0..configuration.len() {
    for j in 0..i {
      let dist = configuration[i].dist(&configuration[j]);
      let coulomb = match context.dimension {
        2 => -dist.ln(),
        3 => 1.0 / dist,
        4 => COULOMB_4D / (dist * dist),
        _ => panic!("unsupported dimension"),
      };
      energy += coulomb * match (&context.molecule[i], &context.molecule[j]) {
        (Particle::Electron, Particle::Electron) => 1.0,
        (Particle::Nucleus{z,a,..}, Particle::Electron) => (if context.dimension > 3 {a * bessel_K1_x(dist)} else {0.0}) - z,
        (Particle::Electron, Particle::Nucleus{z,a,..}) => (if context.dimension > 3 {a * bessel_K1_x(dist)} else {0.0}) - z,
        (Particle::Nucleus{z:z0,a:a0,..}, Particle::Nucleus{z:z1,a:a1,..}) => (if context.dimension > 3 {a0 * a1 * bessel_K1_x(dist) * S_CHARGE_RATIO} else {0.0}) + z0 * z1,
      };
    }
  }
  energy
}

impl Measurement {
  pub fn measure(&self, configuration: &[Position], context: &Context) -> f64 {
    match &self {
      Measurement::Distance => configuration[0].dist(&configuration[1]),
      Measurement::Dipole => zip(configuration, context.molecule).map(|(x,p)| x.x * particle_charge(p)).sum(), // Unlike the Haskell version, this requires the total charge to be 0.
      Measurement::Potential => potential_energy(configuration, context),
      Measurement::One => 1.0,
    }
  }
  
  pub fn from_char(c : char) -> Measurement {
    match c {
      'r' => Measurement::Distance,
      'd' => Measurement::Dipole,
      'v' => Measurement::Potential,
      '1' => Measurement::One,
      _ => panic!("Unrecognised measurement name: {}", c),
    }
  }
}

pub fn suitable_step_size(configuration: &[Position], context: &Context) -> f64 {
  let mut ordinary_bound_inv = 0.0; // This bound limits the first term in the energy error's power series
  let mut min_radius2 = f64::INFINITY; // This is to ensure that particles don't travel as far as the distance between them in one step.
  for i in 0..configuration.len() {
    for j in 0..i {
      let r2 = configuration[i].dist2(&configuration[j]);
      min_radius2 = min_radius2.min(r2);
      ordinary_bound_inv += r2.powi((1 - (context.dimension as i32)).into());
    }
  }
  (context.required_error / ordinary_bound_inv).powf(1.0/3.0).min(0.07 * min_radius2)
}

pub fn show_config(configuration : &[Position]) -> String {
  let mut result = String::from("");
  for p0 in configuration {
    result.push_str("[");
    for p1 in configuration {
      result.push_str(&format!("{:.2e}, ", p0.dist(p1)));
    }
    result.truncate(result.len() - 2); // Remove ", "
    result.push_str("]\n");
  }
  result.pop(); // Remove the trailing newline.
  result
}
