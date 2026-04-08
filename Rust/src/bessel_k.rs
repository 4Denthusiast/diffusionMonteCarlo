use std::array;
use std::sync::LazyLock;

// Following https://arxiv.org/abs/1209.1547
fn bessel_K(ν: f64, z: f64) -> f64 {
  #![allow(non_snake_case)]
  if z == 0.0 {
    return f64::INFINITY;
  }
  let mut prev_result : f64 = -2.0;
  let mut result = -1.0;
  let mut h = 1.0;
  while (result - prev_result).abs() > 1e-11 * prev_result.abs() {
    h /= 2.0;
    prev_result = result;
    result = 0.0;
    let start = (1..).filter(|n| (ν*(*n as f64)*h).cosh() * (-z * ((*n as f64)*h).cosh()).exp() < (-z).exp() * 1e-20).next().unwrap();
    for n in (1..start).rev() {
      result += (ν*(n as f64)*h).cosh() * (-z * ((n as f64)*h).cosh()).exp();
    }
    result += (-z).exp() / 2.0;
    result *= h;
  }
  result
}

struct Quartic{a:f64, b:f64, c:f64, d:f64, e:f64}

impl Quartic {
  fn eval(&self, x:f64) -> f64 {
    self.a + x*(self.b + x*(self.c + x*(self.d + x*self.e)))
  }
}

// taylor series of x^2 K_1(x^2)
fn bessel_K1_x2_taylor_series(x: f64) -> Quartic {
  #![allow(non_snake_case)]
  let k0 = bessel_K(0.0, x*x);
  let k1 = bessel_K(1.0, x*x);
  Quartic {
    a: x*x * k1,
    b: -2.0 * x.powi(3) * k0,
    c: (4.0 * x.powi(4) * k1 - 6.0 * x*x * k0)/2.0,
    d: (20.0 * x.powi(3) * k1 - (8.0 * x.powi(5) + 12.0 * x) * k0)/6.0,
    e: ((16.0 * x.powi(6) + 44.0*x*x) * k1 - (80.0 * x.powi(4) + 12.0) * k0)/24.0,
  }
}

const GRID_FINENESS : f64 = 1e-3;
const GRID_MAX : f64 = 6.0;
const GRID_POINTS : usize = (GRID_MAX/GRID_FINENESS).ceil() as usize;

static SERIES_TABLE : LazyLock<[Quartic; GRID_POINTS]> = LazyLock::new(|| array::from_fn(|n| bessel_K1_x2_taylor_series((n as f64 + 0.5)*GRID_FINENESS)));

// x K1(x), using the memoized table of taylor expansions.
pub fn bessel_K1_x(x: f64) -> f64 {
  #![allow(non_snake_case)]
  if x >= GRID_MAX * GRID_MAX {
    return 0.0;
  } else if x < 0.0 {
    panic!("Bessel K cannot be evaluated for negative arguments.");
  }
  let i = (x.sqrt()/GRID_FINENESS - 0.5).round();
  let dx = x.sqrt() - (i + 0.5) * GRID_FINENESS;
  (*SERIES_TABLE)[i as usize].eval(dx)
}
