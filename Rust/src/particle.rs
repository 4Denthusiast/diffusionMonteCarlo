use std::fmt::{self, Display};

#[derive(Clone, PartialEq)]
pub enum Particle {
  Electron,
  Nucleus {z:f64, a:f64, m:f64},
}

pub fn particle_mass(particle : &Particle) -> f64 {
  match particle {
    Particle::Electron => 1.0,
    Particle::Nucleus{m,..} => *m,
  }
}

pub fn particle_charge(particle : &Particle) -> f64 {
  match particle {
    Particle::Electron => -1.0,
    Particle::Nucleus{z,..} => *z,
  }
}

impl Display for Particle {
  fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    match &self {
      Particle::Electron => write!(f, "e"),
      Particle::Nucleus{z,a,m} => if m.is_infinite() {
          write!(f, "Z{}A{}", z, a)
        } else {
          write!(f, "Z{}A{}m{}", z, a, m)
        },
    }
  }
}
