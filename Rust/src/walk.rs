use crate::particle::*;
use crate::configuration::*;
use crate::variance::*;

use rand::{self, RngExt};
use rand_distr::{self, Distribution};
use std::iter;
use std::iter::zip;

#[derive(Copy, Clone)]
pub struct Walker {
  pub amplitude: f64,
  pub local_time: f64,
  pub measurement_values: [f64; MEASUREMENT_COUNT],
}

pub struct Context<'a> {
  pub dimension: u8,
  pub molecule: &'a [Particle],
  pub delta_time: f64,
  pub set_point: f64,
  pub required_error: f64,
  pub measurements_required: Vec<Measurement>,
  pub measurement_fade: f64,
}

pub struct PopulationState {
  pub walker_set: Vec<Walker>,
  pub positions: Vec<Position>,
  pub energy: f64, // includes an adjustment term to shift the population back towards the set point
  pub variance: Variance,
  pub measurement_variances: Vec<Variance>,
}

impl PopulationState {
  pub fn population(&self) -> u32 {
    self.walker_set.len() as u32
  }

  pub fn total_amplitude(&self) -> f64 {
    self.walker_set.iter().map(|w| w.amplitude).sum()
  }

  pub fn energy(&self, ctx : &Context) -> f64 {
    -self.variance.mean.ln() / ctx.delta_time
  }
  
  pub fn energy_std_dev(&self, ctx : &Context) -> f64 {
    self.variance.std_dev / ctx.delta_time
  }
}

pub fn initial_pop_state(positions: &Vec<Vec<Position>>, context : &Context) -> PopulationState {
  let n_measurements = context.measurements_required.len();
  PopulationState {
    walker_set: vec![ Walker {
        amplitude: 1.0,
        local_time: 0.0,
        measurement_values: [0.0; MEASUREMENT_COUNT],
      };
      positions.len()
    ],
    positions: positions.iter().flatten().map(|p| p.clone()).collect(),
    energy: 0.0,
    variance: Variance::empty(),
    measurement_variances: vec![Variance::empty(); n_measurements],
  }
}

pub fn random_position<R : RngExt>(start : Position, stddev : f64, dimension : u8, rng : &mut R) -> Position {
  let dist = rand_distr::Normal::new(0.0, stddev).unwrap();
  Position {
    x : start.x + dist.sample(rng),
    y : start.y + dist.sample(rng),
    z : start.z + if dimension >= 3 {dist.sample(rng)} else {0.0},
    w : start.w + if dimension >= 4 {dist.sample(rng)} else {0.0},
  }
}

fn step_walkers<R : RngExt>(mut pop : PopulationState, ctx : &Context, rng : &mut R) -> PopulationState {
  // TODO: don't allocate new vectors every time.
  let old_walker_set = pop.walker_set;
  let old_positions = pop.positions;
  pop.walker_set = Vec::new();
  pop.positions = Vec::new();
  for (walker, configuration) in zip(
    old_walker_set.iter(),
    old_positions.as_slice().chunks_exact(ctx.molecule.len())
  ) {
    step_walker(*walker, configuration, &mut pop, ctx, rng, 0);
  }
  pop
}

fn step_walker<R : RngExt>(mut walker : Walker, configuration : &[Position], pop : &mut PopulationState, ctx : &Context, rng : &mut R, depth : u32) {
  if walker.local_time >= 0.0 {
    walker.local_time -= ctx.delta_time;
    for mv in &mut walker.measurement_values {
      *mv *= (-ctx.delta_time * ctx.measurement_fade).exp();
    }
    pop.walker_set.push(walker);
    pop.positions.extend(configuration);
  } else {
    // move walker
    let mut new_configuration = configuration.to_vec(); // TODO: get rid of this allocation.
    let dt = suitable_step_size(configuration, ctx).min(-walker.local_time);
    for (r, p) in zip(&mut new_configuration, ctx.molecule) {
      // TODO: check if omitting infinite-mass particles entirely speeds this step up much.
      *r = random_position(*r, (dt / particle_mass(p)).sqrt(), ctx.dimension, rng);
    }
    let v = (potential_energy(new_configuration.as_slice(), ctx) + potential_energy(configuration, ctx))/2.0;
    walker.amplitude *= (-dt * (v - pop.energy)).exp();
    walker.local_time += dt;
    for (mv, m) in zip(&mut walker.measurement_values, iter::once(&Measurement::One).chain(&ctx.measurements_required)) {
      *mv += dt * m.measure(&new_configuration[..], ctx); // TODO: consider using the average of new_configuration and old_configuration or something.
    }
    
    // split walker
    const MIN_AMPLITUDE : f64 = 0.5;
    if walker.amplitude < MIN_AMPLITUDE {
      if rng.random_bool(walker.amplitude) {
        walker.amplitude = 1.0;
      } else {
        return; // delete the walker.
      }
    }
    let split_into = if walker.amplitude > 2.0 { 2 } else { 1 };
    walker.amplitude /= split_into as f64;
    
    // recurse
    for _ in 0..split_into {
      step_walker(walker, new_configuration.as_slice(), pop, ctx, rng, depth + 1);
    }
  }
}

fn double_walker_set(pop : &mut PopulationState) {
  pop.walker_set.extend_from_within(..);
  pop.positions.extend_from_within(..);
}

fn trim_walker_set(pop : &mut PopulationState, ctx : &Context) {
  pop.walker_set.truncate(pop.walker_set.len()/2);
  pop.positions.truncate(pop.walker_set.len() * ctx.molecule.len());
}

fn measurement_totals(pop : &PopulationState, ctx : &Context) -> Vec<f64> {
  let mut result = vec![0.0; ctx.measurements_required.len()];
  for walker in &pop.walker_set {
    let norm = walker.measurement_values[0]; // Fixed as Measurement::One
    for i in 0..result.len() {
      result[i] += walker.measurement_values[i+1] / norm;
    }
  }
  result
}

pub fn step<R : RngExt>(mut pop : PopulationState, ctx : &Context, rng : &mut R) -> PopulationState{
  let start_pop = pop.total_amplitude();
  if start_pop == 0.0 {
    panic!("Sample population extinct.");
  }
  pop = step_walkers(pop, ctx, rng);
  let end_pop = pop.total_amplitude();
  let growth_estimate = (end_pop / start_pop) * (-pop.energy * ctx.delta_time).exp();
  pop.variance.add_data_point(end_pop / ctx.set_point, growth_estimate);
  let measurements = measurement_totals(&pop, ctx);
  for (v, m) in zip(&mut pop.measurement_variances, measurements) {
    v.add_data_point(end_pop / ctx.set_point, m / end_pop);
  }
  let relaxation_time = ctx.delta_time * pop.variance.time_at_bound(0.5);
  pop.energy = pop.energy(ctx) + (ctx.set_point / end_pop).ln() / relaxation_time;
  if end_pop < ctx.set_point / 4.0 {
    double_walker_set(&mut pop);
  } else if end_pop > ctx.set_point * 4.0 {
    trim_walker_set(&mut pop, ctx);
  }
  pop
}

pub fn display_popstate<R : RngExt>(pop : &PopulationState, ctx : &Context, rng : &mut R) {
  for (v, m) in zip(&pop.measurement_variances, &ctx.measurements_required) {
    println!("{}: {:.5e} ± {:.2e}", m, v.mean, v.std_dev);
  }
  let config_index = rng.random_range(0..pop.walker_set.len());
  let random_configuration = &pop.positions[(config_index * ctx.molecule.len())..((config_index+1) * ctx.molecule.len())];
  println!("pop: {}, pop(w): {:.1}, energy: {:.5e} ± {:.2e}\n{}", pop.population(), pop.total_amplitude(), pop.energy(ctx), pop.energy_std_dev(ctx), show_config(random_configuration));
}
