#![allow(mixed_script_confusables)]
#![allow(confusable_idents)]

mod particle;
mod bessel_k;
mod configuration;
mod variance;
mod walk;

use crate::particle::*;
use crate::configuration::*;
use crate::variance::*;
use crate::walk::*;

use std::env;
use std::io::{stdout, Stdout, Write};
use std::iter::zip;
use crossterm::{self, QueueableCommand, terminal::ClearType};
use rand::{self, RngExt, rngs::SmallRng};
use regex::Regex;

enum AtomSep {
  Atom {z : u32, zf : f64, a : u32, m : f64, q : i32},
  Fixed (f64),
  Random (f64),
  Offset {mean : f64, std_dev : f64},
}

#[derive(Clone, Copy, PartialEq, PartialOrd)]
enum Verbosity {
  Quiet,
  Normal,
  Verbose,
}

struct ExecutionParameters {
  dimension : u8,
  atom_seps : Vec<AtomSep>,
  molecule_name : String,
  required_error : f64,
  time_step : f64,
  walker_count : usize,
  prep_steps : u32,
  measurements : Vec<Measurement>,
  measurement_fade : f64,
  verbosity : Verbosity,
}

impl Default for ExecutionParameters {
  fn default() -> ExecutionParameters {
    ExecutionParameters {
      dimension : 3,
      atom_seps : vec![],
      molecule_name : String::from(""),
      required_error : 0.01,
      time_step : 0.05,
      walker_count : 200,
      prep_steps : 1000,
      measurements : vec![],
      measurement_fade : 0.0,
      verbosity : Verbosity::Normal,
   }
 }
}

fn main() {
  let mut rng : SmallRng = rand::make_rng();
  let params = parse_arguments();
  let molecule : Vec<Particle> = make_molecule(&params.atom_seps, params.dimension, &mut rng).into_iter().map(|(p,_)| p).collect();
  let ctx = Context {
    dimension : params.dimension,
    molecule : &molecule[..],
    delta_time : params.time_step,
    set_point : params.walker_count as f64,
    required_error : params.required_error,
    measurements_required : params.measurements,
    measurement_fade : params.measurement_fade,
  };
  let configurations = (0..params.walker_count).map(|_| make_molecule(&params.atom_seps, params.dimension, &mut rng).into_iter().map(|(_,r)| r).collect()).collect();
  let mut population_state = initial_pop_state(&configurations, &ctx);
  let mut out = stdout();
  if params.verbosity < Verbosity::Verbose {
    out.queue(crossterm::cursor::SavePosition).unwrap();
  }
  let mut n_temp_lines : u16 = 0;
  for i in 0..params.prep_steps {
    population_state = step(population_state, &ctx, &mut rng);
    display_popstate(&population_state, &ctx, &mut rng, params.verbosity, &mut out, &mut n_temp_lines);
    if params.verbosity >= Verbosity::Normal {
      out.write_all(format!("Preparing...{}/{}\n", i, params.prep_steps).as_bytes()).unwrap();
      n_temp_lines += 1;
    }
  }
  reset_iteration(&mut population_state);
  while population_state.energy_std_dev(&ctx) > params.required_error {
    population_state = step(population_state, &ctx, &mut rng);
    display_popstate(&population_state, &ctx, &mut rng, params.verbosity, &mut out, &mut n_temp_lines);
  }
  if params.verbosity == Verbosity::Quiet {
    out.write_all(params.molecule_name.as_bytes()).unwrap();
    out.write_all(": ".as_bytes()).unwrap();
    for v in &population_state.measurement_variances {
      out.write_all(format!("{:.5e} ± {:.2e}, ", v.mean, v.std_dev).as_bytes()).unwrap();
    }
    out.write_all(format!("{:.5e} ± {:.2e}\n", population_state.energy(&ctx), population_state.energy_std_dev(&ctx)).as_bytes()).unwrap();
  }
}

fn parse_arguments() -> ExecutionParameters {
  let mut params = ExecutionParameters::default();
  let mut args = env::args();
  args.next(); // args[0] is the program name.
  while let Some(arg) = args.next() {
    match &arg[..] {
      "-d" => params.dimension        = args.next().expect("missing argument after -d").parse().expect("couldn't understand -d argument"),
      "-r" => params.required_error   = args.next().expect("missing argument after -r").parse().expect("couldn't understand -r argument"),
      "-t" => params.time_step        = args.next().expect("missing argument after -t").parse().expect("couldn't understand -t argument"),
      "-w" => params.walker_count     = args.next().expect("missing argument after -w").parse().expect("couldn't understand -w argument"),
      "-p" => params.prep_steps       = args.next().expect("missing argument after -p").parse().expect("couldn't understand -p argument"),
      "-f" => params.measurement_fade = args.next().expect("missing argument after -f").parse().expect("couldn't understand -f argument"),
      "-m" => params.measurements     = args.next().expect("missing argument after -m").chars().map(Measurement::from_char).collect(),
      "--verbose" => params.verbosity = Verbosity::Verbose,
      "--quiet" => params.verbosity = Verbosity::Quiet,
      arg => {
        params.atom_seps.push(parse_atom_sep(arg).expect(&format!("couldn't understand argument \"{}\"", &arg[..])));
        if params.molecule_name.len() > 0 {
          params.molecule_name.push(' ');
        }
        params.molecule_name.push_str(arg)
      },
    };
  }
  if params.atom_seps.len() == 0 {
    params.atom_seps.push(AtomSep::Atom{z:1,zf:0.0,a:1,m:f64::INFINITY,q:0});
    params.molecule_name = String::from("H");
  }
  params
}

fn parse_atom_sep(s : &str) -> Option<AtomSep> {
  if let Ok(as_number) = s.parse::<f64>() {
    return Some(AtomSep::Fixed(as_number));
  }
  if let Some((s1,s2)) = s.split_once('*') {
    if let Ok(x1) = s1.parse::<f64>() {
      if s2.len() == 0 {
        return Some(AtomSep::Random(x1));
      } else if let Ok(x2) = s2.parse::<f64>() {
        return Some(AtomSep::Offset{mean : x1, std_dev : x2});
      } else {
        return None;
      }
    }
  }
  let regex_matches = Regex::new(r"^(n|H|He|Li|Be)([#♯b♭]?)(\d+)?(?:m([0-9\\.]+))?([\+\-]*)$").unwrap().captures(s)?;
  let z = match &regex_matches[1] {
    "n" => Some(0),
    "H" => Some(1),
    "He" => Some(2),
    "Li" => Some(3),
    "Be" => Some(4),
    _ => None,
  }?;
  let zf = match &regex_matches[2] {
    "#" => 1.0/3.0,
    "♯" => 1.0/3.0,
    "b" => -1.0/3.0,
    "♭" => -1.0/3.0,
    "" => 0.0,
    _ => unreachable!(),
  };
  let a = match regex_matches.get(3) {
    None => if z <= 1 {1} else {2*z},
    Some(a_s) => a_s.as_str().parse::<u32>().ok()?,
  };
  let m = match regex_matches.get(4) {
    None => f64::INFINITY,
    Some(ms) => ms.as_str().parse::<f64>().ok()?,
  };
  let q = regex_matches[5].matches('+').collect::<Vec<_>>().len() as i32 * 2 - regex_matches[5].len() as i32;
  Some(AtomSep::Atom{z,zf,a,m,q})
}

fn make_molecule<R : RngExt>(atom_seps : &Vec<AtomSep>, dimension : u8, rng : &mut R) -> Vec<(Particle,Position)> {
  let mut molecule = vec![];
  let mut position = Position::default();
  let mut had_sep = true;
  for atom_sep in atom_seps {
    match atom_sep {
      AtomSep::Fixed(x) => {
        had_sep = true;
        position.x += *x;
      }
      AtomSep::Random(r) => {
        had_sep = true;
        position = random_position(position, *r, dimension, rng);
      }
      AtomSep::Offset{mean,std_dev} => {
        had_sep = true;
        position = random_position(position, *std_dev, dimension, rng);
        position.x += mean;
      }
      AtomSep::Atom{z,zf,a,m,q} => {
        if !had_sep {
          position.x += 1.0;
        }
        molecule.push((Particle::Nucleus{z:*z as f64 + *zf, a:*a as f64, m:*m}, position));
        for _ in 0..(*z as i32 - q) {
          molecule.push((Particle::Electron, random_position(position, 4.0, dimension, rng)));
        }
        had_sep = false;
      }
    }
  }
  molecule.sort_by_key(|(p,_)| if *p == Particle::Electron {1} else {0});
  molecule
}

fn reset_iteration(population_state : &mut PopulationState) {
  population_state.variance = Variance::empty();
  for measurement in &mut population_state.measurement_variances {
    *measurement = Variance::empty();
  }
  for walker in &mut population_state.walker_set {
    walker.measurement_values = [0.0; MEASUREMENT_COUNT];
  }
}

fn display_popstate<R : RngExt>(pop : &PopulationState, ctx : &Context, rng : &mut R, verbosity : Verbosity, out : &mut Stdout, n_temp_lines : &mut u16) {
  if verbosity <= Verbosity::Quiet {
    return;
  }
  if verbosity < Verbosity::Verbose && *n_temp_lines > 0{
    out.queue(crossterm::cursor::MoveToPreviousLine(*n_temp_lines)).unwrap();
    out.queue(crossterm::terminal::Clear(ClearType::FromCursorDown)).unwrap();
  }
  *n_temp_lines = 0;
  for (v, m) in zip(&pop.measurement_variances, &ctx.measurements_required) {
    out.write_all(format!("{}: {:.5e} ± {:.2e}\n", m, v.mean, v.std_dev).as_bytes()).unwrap();
    *n_temp_lines += 1;
  }
  let config_index = rng.random_range(0..pop.walker_set.len());
  let random_configuration = &pop.positions[(config_index * ctx.molecule.len())..((config_index+1) * ctx.molecule.len())];
  out.write_all(format!("pop: {}, pop(w): {:.1}, energy: {:.5e} ± {:.2e}\n", pop.population(), pop.total_amplitude(), pop.energy(ctx), pop.energy_std_dev(ctx)).as_bytes()).unwrap();
    *n_temp_lines += 1;
  if verbosity >= Verbosity::Verbose {
    out.write_all(show_config(random_configuration).as_bytes()).unwrap();
    out.write_all("\n".as_bytes()).unwrap();
    *n_temp_lines += random_configuration.len() as u16;
  }
  out.flush().unwrap();
}
