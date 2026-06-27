#![allow(mixed_script_confusables)]
#![allow(confusable_idents)]

mod particle;
mod bessel_k;
mod configuration;
mod ansatz;
mod variance;
mod thread_pool;
mod walk;

use crate::particle::*;
use crate::configuration::*;
use crate::ansatz::*;
use crate::thread_pool::*;
use crate::walk::*;

use std::env;
use std::io::{stdout, Stdout, Write};
use std::iter::zip;
use std::thread::available_parallelism;
use std::time::{Duration, Instant};
use crossterm::{self, QueueableCommand, terminal::ClearType};
use rand::{self, RngExt, rngs::SmallRng};
use regex::Regex;

#[derive(PartialEq)]
enum AtomSep {
  Atom {z : u32, zf : f64, a : u32, m : f64, q : i32},
  Fixed (f64),
  Random (f64),
  Offset {mean : f64, std_dev : f64},
}

#[derive(Clone, Copy, PartialEq, PartialOrd)]
enum Verbosity {
  Benchmark,
  Quiet,
  Normal,
  Verbose,
}

struct ExecutionParameters {
  dimension : u8,
  atom_seps : Vec<Vec<AtomSep>>,
  molecule_names : Vec<String>,
  required_error : f64,
  time_step : f64,
  walker_count : usize,
  prep_steps : u32,
  measurements : Vec<Measurement>,
  measurement_fade : f64,
  verbosity : Verbosity,
  thread_count : usize,
  ansatz : Box<dyn Ansatz>,
  timeout : Option<Duration>,
}

impl Default for ExecutionParameters {
  fn default() -> ExecutionParameters {
    ExecutionParameters {
      dimension : 3,
      atom_seps : vec![vec![]],
      molecule_names : vec![String::from("")],
      required_error : 0.001,
      time_step : 0.05,
      walker_count : 10000,
      prep_steps : 1000,
      measurements : vec![],
      measurement_fade : 0.0,
      verbosity : Verbosity::Normal,
      thread_count : available_parallelism().map_or(1, |c| c.get()),
      ansatz : Box::new(TrivialAnsatz),
      timeout : Option::None,
    }
  }
}

fn benchmark_parameters() -> ExecutionParameters {
  ExecutionParameters {
    atom_seps : vec![vec![
      AtomSep::Atom{z:1, zf:0.0, a:1, m:f64::INFINITY, q:0},
      AtomSep::Fixed(1.4),
      AtomSep::Atom{z:1, zf:0.0, a:1, m:1800.0, q:0},
    ]],
    time_step : 1.0,
    prep_steps : 1000,
    measurements : vec![Measurement::Distance, Measurement::Potential],
    verbosity : Verbosity::Benchmark,
    ..ExecutionParameters::default()
  }
}

fn main() {
  let mut rng : SmallRng = rand::make_rng();
  let params = parse_arguments();
  for (atom_seps, molecule_name) in zip(params.atom_seps, params.molecule_names) {
    let start_time = Instant::now();
    let molecule : Vec<Particle> = make_molecule(&atom_seps, params.dimension, &mut rng).into_iter().map(|(p,_)| p).collect();
    let ctx = Context {
      dimension : params.dimension,
      molecule : molecule,
      delta_time : params.time_step,
      set_point : params.walker_count as f64,
      required_error : params.required_error,
      measurements_required : params.measurements.clone(),
      measurement_fade : params.measurement_fade,
      ansatz : params.ansatz.clone(),
    };
    let configurations = (0..params.walker_count).map(|_| make_molecule(&atom_seps, params.dimension, &mut rng).into_iter().map(|(_,r)| r).collect()).collect();
    let mut population_state = initial_pop_state(&configurations, &ctx, ThreadPool::new(&mut rng, params.thread_count));
    let mut out = stdout();
    if params.verbosity == Verbosity::Normal {
      out.queue(crossterm::cursor::SavePosition).unwrap();
    }
    let mut n_temp_lines : u16 = 0;
    for i in 0..params.prep_steps {
      population_state = step(population_state, &ctx);
      display_popstate(&population_state, &ctx, &mut rng, params.verbosity, &mut out, &mut n_temp_lines);
      if params.verbosity >= Verbosity::Normal {
        out.write_all(format!("Preparing...{}/{}\n", i, params.prep_steps).as_bytes()).unwrap();
        n_temp_lines += 1;
      }
    }
    if params.verbosity == Verbosity::Benchmark {
      // This result is somewhat random, which isn't great for a benchmark. However, removing the randomness would give a false sense of security because it's possible for changes being benchmarked to alter how the fixed sequence of samples is used, thereby switching to, essentially, a new random sample. At least this way, you can run it a few times and get an idea what the variance is.
      println!("time: {:?}", start_time.elapsed());
      return;
    }
    reset_iteration(&mut population_state);
    while population_state.energy_std_dev(&ctx) > params.required_error && params.timeout.is_none_or(|t| start_time.elapsed() < t) {
      population_state = step(population_state, &ctx);
      display_popstate(&population_state, &ctx, &mut rng, params.verbosity, &mut out, &mut n_temp_lines);
    }
    if params.verbosity == Verbosity::Quiet {
      out.write_all(molecule_name.as_bytes()).unwrap();
      out.write_all(": ".as_bytes()).unwrap();
      for v in &population_state.measurement_variances {
        out.write_all(format!("{:.5e} ± {:.2e}, ", v.mean, v.std_dev).as_bytes()).unwrap();
      }
      out.write_all(format!("{:.5e} ± {:.2e}\n", population_state.energy(&ctx), population_state.energy_std_dev(&ctx)).as_bytes()).unwrap();
    }
  }
}

fn parse_arguments() -> ExecutionParameters {
  let mut params = ExecutionParameters::default();
  let mut args = env::args();
  args.next(); // args[0] is the program name.
  let mut is_first = true;
  while let Some(arg) = args.next() {
    match &arg[..] {
      "-d" => params.dimension        = args.next().expect("missing argument after -d").parse().expect("couldn't understand -d argument"),
      "-r" => params.required_error   = args.next().expect("missing argument after -r").parse().expect("couldn't understand -r argument"),
      "-t" => params.time_step        = args.next().expect("missing argument after -t").parse().expect("couldn't understand -t argument"),
      "-w" => params.walker_count     = args.next().expect("missing argument after -w").parse().expect("couldn't understand -w argument"),
      "-p" => params.prep_steps       = args.next().expect("missing argument after -p").parse().expect("couldn't understand -p argument"),
      "-f" => params.measurement_fade = args.next().expect("missing argument after -f").parse().expect("couldn't understand -f argument"),
      "-m" => params.measurements     = args.next().expect("missing argument after -m").chars().map(Measurement::from_char).collect(),
      "-þ" => params.thread_count     = args.next().expect("missing argument after -þ").parse().expect("couldn't understand -þ argument"),
      "-a" => params.ansatz = parse_ansatz(&args.next().expect("missing argument after -a")),
      "--verbose" => params.verbosity = Verbosity::Verbose,
      "--quiet" => params.verbosity = Verbosity::Quiet,
      "--benchmark" => {
        if !is_first {
          println!("Warning: --benchmark option will cause all earlier parameters to be ignored.");
        }
        params = benchmark_parameters();
      },
      "--timeout" => params.timeout = Option::Some(parse_timeout(&args.next().expect("missing argument after --timeout"))),
      "," => {
        if params.atom_seps.last().unwrap().len() == 0 {
          panic!("Unexpected comma, a molecule specification is required first.");
        }
        params.atom_seps.push(vec![]);
        params.molecule_names.push(String::from(""));
      },
      arg => {
        params.atom_seps.last_mut().unwrap().push(parse_atom_sep(arg).expect(&format!("couldn't understand argument \"{}\"", &arg[..])));
        let molecule_name = params.molecule_names.last_mut().unwrap();
        if molecule_name.len() > 0 {
          molecule_name.push(' ');
        }
        molecule_name.push_str(arg)
      },
    };
    is_first = false;
  }
  if params.atom_seps.len() == 1 && params.atom_seps[0].len() == 0 {
    params.atom_seps[0].push(AtomSep::Atom{z:1,zf:0.0,a:1,m:f64::INFINITY,q:0});
    params.molecule_names[0] = String::from("H");
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

fn parse_ansatz(string : &str) -> Box<dyn Ansatz> {
  match string {
    "1" => Box::new(TrivialAnsatz),
    "sin" => Box::new(SineTestAnsatz),
    "cusps3d" => Box::new(Cusps3DSimple),
    _ => panic!("Unrecognized ansatz name {string}"),
  }
}

fn parse_timeout(string : &str) -> Duration {
  let last = string.chars().next_back().expect("Empty timeout argument"); // It is probably impossible for a command-like argument to be empty anyway.
  let without_last = &string[0..string.len()-1];
  match last {
    's' => Duration::from_secs(without_last.parse().unwrap()),
    'm' => Duration::from_mins(without_last.parse().unwrap()),
    'h' => Duration::from_hours(without_last.parse().unwrap()),
    'd' => Duration::from_hours(without_last.parse().unwrap()) * 24,
    _ => Duration::from_secs(string.parse().unwrap()),
  }
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
        position = position + random_position(*r, dimension, rng);
      }
      AtomSep::Offset{mean,std_dev} => {
        had_sep = true;
        position = position + random_position(*std_dev, dimension, rng);
        position.x += mean;
      }
      AtomSep::Atom{z,zf,a,m,q} => {
        if !had_sep {
          position.x += 1.0;
        }
        molecule.push((Particle::Nucleus{z:*z as f64 + *zf, a:*a as f64, m:*m}, position));
        for _ in 0..(*z as i32 - q) {
          molecule.push((Particle::Electron, position + random_position(4.0, dimension, rng)));
        }
        had_sep = false;
      }
    }
  }
  molecule.sort_by_key(|(p,_)| if *p == Particle::Electron {1} else {0});
  molecule
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
  let random_configuration = pop.random_configuration(ctx, rng);
  out.write_all(format!("pop: {}, pop(w): {:.1}, energy: {:.5e} ± {:.2e}\n", pop.population(), pop.total_amplitude(), pop.energy(ctx), pop.energy_std_dev(ctx)).as_bytes()).unwrap();
    *n_temp_lines += 1;
  if verbosity >= Verbosity::Verbose {
    out.write_all(show_config(random_configuration).as_bytes()).unwrap();
    out.write_all("\n".as_bytes()).unwrap();
    *n_temp_lines += random_configuration.len() as u16;
  }
  out.flush().unwrap();
}
