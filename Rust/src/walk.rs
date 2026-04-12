use crate::particle::*;
use crate::configuration::*;
use crate::variance::*;
use crate::thread_pool::*;

use rand::{self, RngExt, rngs::SmallRng};
use rand_distr::{self, Distribution};
use std::iter;
use std::iter::zip;
use std::sync::mpsc;

#[derive(Copy, Clone)]
pub struct Walker {
  pub amplitude: f64,
  pub local_time: f64,
  pub measurement_values: [f64; MEASUREMENT_COUNT],
  prev_potential: f64,
}

#[derive(Clone)]
pub struct Context {
  pub dimension: u8,
  pub molecule: Vec<Particle>,
  pub delta_time: f64,
  pub set_point: f64,
  pub required_error: f64,
  pub measurements_required: Vec<Measurement>,
  pub measurement_fade: f64,
}

#[derive(Clone)]
struct PopChunk {
  walker_set: Vec<Walker>,
  positions: Vec<Position>,
}

const CHUNK_SIZE : usize = 100;

pub struct PopulationState {
  chunks: Vec<PopChunk>,
  pub energy: f64, // includes an adjustment term to shift the population back towards the set point
  pub variance: Variance,
  pub measurement_variances: Vec<Variance>,
  pub thread_pool: ThreadPool<SmallRng>,
}

impl PopulationState {
  pub fn population(&self) -> u32 {
    self.chunks.iter().map(|c| c.walker_set.len() as u32).sum()
  }

  pub fn total_amplitude(&self) -> f64 {
    self.chunks.iter().map(|c| c.walker_set.iter().map(|w| w.amplitude).sum::<f64>()).sum()
  }

  pub fn energy(&self, ctx : &Context) -> f64 {
    -self.variance.mean.ln() / ctx.delta_time
  }
  
  pub fn energy_std_dev(&self, ctx : &Context) -> f64 {
    self.variance.std_dev / ctx.delta_time
  }
  
  pub fn random_configuration<R : RngExt>(&self, ctx : &Context, rng : &mut R) -> &[Position] {
    let chunk = &self.chunks[rng.random_range(0..self.chunks.len())];
    let config_index = rng.random_range(0..chunk.walker_set.len());
    &chunk.positions[(config_index * ctx.molecule.len())..((config_index+1) * ctx.molecule.len())]
  }
}

pub fn initial_pop_state(positions: &Vec<Vec<Position>>, context : &Context, thread_pool : ThreadPool<SmallRng>) -> PopulationState {
  let n_measurements = context.measurements_required.len();
  PopulationState {
    chunks: positions.chunks(CHUNK_SIZE).map(|ps| PopChunk {
      walker_set: vec![ Walker {
          amplitude: 1.0,
          local_time: 0.0,
          measurement_values: [0.0; MEASUREMENT_COUNT],
          prev_potential: 0.0, // This is kind of incorrect, but it doesn't matter to the final result.
        };
        ps.len()
      ],
      positions: ps.iter().flatten().map(|p| p.clone()).collect(),
    }).collect(),
    energy: 0.0,
    variance: Variance::empty(),
    measurement_variances: vec![Variance::empty(); n_measurements],
    thread_pool,
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

fn step_walkers(mut pop : PopulationState, ctx : &Context) -> PopulationState {
  let (in_send, in_rcv) = mpsc::channel();
  for chunk in pop.chunks {
    let thread_send = in_send.clone();
    let thread_ctx = ctx.clone();
    pop.thread_pool.submit(move |mut rng| {
      let mut configuration_stack = Vec::new();
      let mut new_chunk = PopChunk {
        walker_set : vec![],
        positions : vec![],
      };
      for (walker, configuration) in zip(
        chunk.walker_set.iter(),
        chunk.positions.as_slice().chunks_exact(thread_ctx.molecule.len())
      ) {
        configuration_stack.extend_from_slice(configuration);
        step_walker(*walker, &mut configuration_stack, &mut new_chunk, pop.energy, &thread_ctx, &mut rng, 0);
      }
      thread_send.send(new_chunk).unwrap();
    });
  }
  drop(in_send);
  pop.chunks = rebalance_chunks(in_rcv.iter(), ctx);
  pop
}

fn step_walker<R : RngExt>(mut walker : Walker, stack : &mut Vec<Position>, out_chunk : &mut PopChunk, energy : f64, ctx : &Context, rng : &mut R, depth : u32) {
  let n_particles = ctx.molecule.len();
  let config_range = (stack.len()-n_particles)..stack.len();
  let configuration = &mut stack[config_range.clone()];
  if walker.local_time >= 0.0 {
    walker.local_time -= ctx.delta_time;
    for mv in &mut walker.measurement_values {
      *mv *= (-ctx.delta_time * ctx.measurement_fade).exp();
    }
    out_chunk.walker_set.push(walker);
    out_chunk.positions.extend_from_slice(&configuration);
    stack.truncate(stack.len() - n_particles);
  } else {
    // move walker
    let dt = suitable_step_size(&configuration, ctx).min(-walker.local_time);
    for (r, p) in zip(&mut *configuration, &ctx.molecule) {
      if particle_mass(p).is_finite() {
        *r = random_position(*r, (dt / particle_mass(p)).sqrt(), ctx.dimension, rng);
      }
    }
    let v = potential_energy(&configuration, ctx);
    walker.amplitude *= (-dt * ((v + walker.prev_potential)/2.0 - energy)).exp();
    walker.prev_potential = v;
    walker.local_time += dt;
    for (mv, m) in zip(&mut walker.measurement_values, iter::once(&Measurement::One).chain(&ctx.measurements_required)) {
      *mv += dt * m.measure(&configuration, ctx); // TODO: consider using the average of new_configuration and old_configuration or something.
    }
    
    // split walker
    const MIN_AMPLITUDE : f64 = 0.5;
    if walker.amplitude < MIN_AMPLITUDE {
      if rng.random_bool(walker.amplitude) {
        walker.amplitude = 1.0;
      } else {
        stack.truncate(stack.len() - n_particles);
        return; // delete the walker.
      }
    }
    let split_into = if walker.amplitude > 2.0 {
      stack.extend_from_within(config_range.clone());
      2
    } else { 1 };
    walker.amplitude /= split_into as f64;
    
    // recurse
    for _ in 0..split_into {
      step_walker(walker, stack, out_chunk, energy, ctx, rng, depth + 1);
    }
  }
}

fn rebalance_chunks(iter : impl Iterator<Item = PopChunk>, ctx : &Context) -> Vec<PopChunk> {
  let mut short_chunk = None;
  let mut result = vec![];
  for mut chunk in iter {
    if chunk.walker_set.len() < CHUNK_SIZE / 2 {
      match short_chunk {
        None => short_chunk = Some(chunk),
        Some(other) => {
          short_chunk = None;
          chunk.walker_set.extend(other.walker_set);
          chunk.positions.extend(other.positions);
          result.push(chunk);
        }
      }
    } else if chunk.walker_set.len() > CHUNK_SIZE * 2 {
      let other = PopChunk {
        walker_set : chunk.walker_set.drain(CHUNK_SIZE..).collect(),
        positions  : chunk.positions.drain((CHUNK_SIZE * ctx.molecule.len())..).collect(),
      };
      result.push(chunk);
      result.push(other);
    } else {
      result.push(chunk);
    }
  }
  match short_chunk {
    None => (),
    Some(chunk) => result.push(chunk),
  };
  result
}

fn double_walker_set(pop : &mut PopulationState) {
  pop.chunks.extend_from_within(..);
}

fn trim_walker_set(pop : &mut PopulationState) {
  pop.chunks.truncate((pop.chunks.len()+1)/2); // +1 so that in the unlikely case that there's only one chunk, it isn't deleted.
}

fn measurement_totals(pop : &PopulationState, ctx : &Context) -> Vec<f64> {
  let mut result = vec![0.0; ctx.measurements_required.len()];
  for chunk in &pop.chunks {
    for walker in &chunk.walker_set {
      let norm = walker.measurement_values[0]; // Fixed as Measurement::One
      for i in 0..result.len() {
        result[i] += walker.amplitude * walker.measurement_values[i+1] / norm;
      }
    }
  }
  result
}

pub fn step(mut pop : PopulationState, ctx : &Context) -> PopulationState{
  let start_pop = pop.total_amplitude();
  if start_pop == 0.0 {
    panic!("Sample population extinct.");
  }
  pop = step_walkers(pop, ctx);
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
    trim_walker_set(&mut pop);
  }
  pop
}

pub fn reset_iteration(population_state : &mut PopulationState) {
  population_state.variance = Variance::empty();
  for measurement in &mut population_state.measurement_variances {
    *measurement = Variance::empty();
  }
  for chunk in &mut population_state.chunks {
    for walker in &mut chunk.walker_set {
      walker.measurement_values = [0.0; MEASUREMENT_COUNT];
    }
  }
}
