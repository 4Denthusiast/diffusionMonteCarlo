#[derive(Clone)]
struct Layer {
  latest: Option<(f64,f64)>,
  weight: f64,
  total: f64,
  square: f64,
}

#[derive(Clone)]
pub struct Variance {
  pub mean: f64,
  pub std_dev: f64,
  layers: Vec<Layer>,
}

impl Variance {
  pub fn empty() -> Variance {
    Variance {
      mean: 0.0,
      std_dev: f64::INFINITY,
      layers: vec![],
    }
  }
  
  pub fn add_data_point(&mut self, weight: f64, x: f64) {
    if x.is_nan() {
      return;
    }
    let mut weight = weight;
    let mut x = x;
    let mut i = 0;
    loop {
      if i >= self.layers.len() {
        self.layers.push(Layer{latest: None, weight: 0.0, total: 0.0, square: 0.0});
      }
      self.layers[i].weight += weight;
      self.layers[i].total  += x * weight;
      self.layers[i].square += x * x * weight;
      match self.layers[i].latest {
        None => {
          self.layers[i].latest = Some((x, weight));
          break;
        },
        Some((x1, weight1)) => {
          self.layers[i].latest = None;
          let weight0 = weight;
          weight = (weight0 + weight1) / 2.0;
          x = (x * weight0 + x1 * weight1) / (weight0 + weight1);
        }
      }
      i += 1;
    }
    
    // recompute mean and std_dev
    self.mean = self.layers[0].total / self.layers[0].weight;
    
    if self.layers.len() < 5 {
      self.std_dev = f64::INFINITY; // There isn't enough information for a good estimate, so use a conservative one instead.
      return;
    }
    let mut sum_c = 0.0;
    let mut sum_sc = 0.0;
    let mut sum_vc = 0.0;
    let mut sum_svc = 0.0;
    let mut sum_ssc = 0.0;
    let mut scale_log = 0.0;
    let mut n_samples : f64 = 0.0;
    let last_used_layer = (self.layers.len() / 2).min(self.layers.len()-5); // The variance typically follows roughly a power law in frequency, except that there's a single break where the power differs (I assume at the frequency corresponding to the lowest excitation energy of the system). The overall variance is linearly extrapolated on a log-log plot using the lower-frequency half of the layers, so that eventually the layers with frequencies above the break are ignored.
    for i in (last_used_layer..self.layers.len()).rev() {
      let layer = &self.layers[i];
      n_samples *= 2.0;
      if layer.latest.is_some() {
        n_samples += 1.0;
      }
      let variance = (layer.square - layer.total*layer.total/layer.weight)/layer.weight * n_samples/(n_samples-1.0);
      let confidence = n_samples - 1.0;
      if confidence > 0.0 {
        sum_c += confidence;
        sum_sc += scale_log * confidence;
        sum_vc += variance.log2() * confidence;
        sum_svc += scale_log * variance.log2() * confidence;
        sum_ssc += scale_log * scale_log * confidence;
      }
      scale_log += 1.0;
    }
    let average_variance = sum_vc / sum_c;
    let average_scale = sum_sc / sum_c;
    let final_scale = scale_log - n_samples.log2();
    let gradient = (sum_svc - sum_vc * sum_sc / sum_c) / (sum_ssc - sum_sc * sum_sc / sum_c);
    let intercept = average_variance + (final_scale - average_scale) * gradient;
    self.std_dev = (intercept/2.0).exp2();
  }
  
  pub fn time_at_bound(&self, b: f64) -> f64 {
    for i in 0..self.layers.len() {
      let layer = &self.layers[i];
      let layer_variance = (layer.square - layer.total*layer.total / layer.weight) / (layer.weight - 1.0).max(0.0);
      //println!("layer {}, variance {}", i, layer_variance);
      if !(layer_variance < b * b / ((1 << (2*i)) as f64)) {
        return (1 << i) as f64;
      }
    }
    return (1 << self.layers.len()) as f64;
  }
}
