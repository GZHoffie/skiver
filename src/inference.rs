use std::collections::HashMap;


use crate::types::*;
use crate::kvmer::KVmerStats;
use crate::constants::*;
use crate::cmdline::AnalyzeArgs;
use crate::huber::*;

use log::info;
use rand::Rng;

/**
 * Info of the estimated error rate and error spectrum
 * Each field is a tuple of (estimate, (5th_percentile, 95th_percentile))
 * where the confidence interval is estimated using bootstrap
 */
pub struct ReadPositionCalibration {
    pub index: u32,
    pub from_start: bool,
    pub num_correct: u64,
    pub num_error: u64,
}

pub struct QscoreCalibration {
    pub qscore: u8,
    pub num_correct: u64,
    pub num_error: u64,
    pub error_rate: f64,
    pub ci_lower: f64,
    pub ci_upper: f64,
}

pub struct ErrorSpectrum {
    // estimated Weibull parameters
    pub estimated_lambda: (f32, (f32, f32)),
    pub estimated_beta: (f32, (f32, f32)),

    // estimated error rates
    pub per_base_error_rate: (f32, (f32, f32)),
    pub effective_error_rate: (f32, (f32, f32)),

    // coverage information
    pub key_coverage: (f32, (f32, f32)),
    pub estimated_coverage: (f32, (f32, f32)),

    // error spectrum
    pub snp_rate: HashMap<(EditOperation, u8, u8), u32>,

    pub bidirectional: bool,
}


pub enum RatioEstimationMethod {
    Slope,
    LinearFit,
    RatioMean,
    SumRatio,
}

pub struct ErrorAnalyzer {
    /* 
    pub k: u8,

    pub bidirectional: bool,
    pub exclude_outliers: bool,
    pub outlier_threshold: f32,
    pub ratio_method: RatioEstimationMethod,

    // bootstrap parameters
    pub num_experiments: u32,
    pub bootstrap_sample_rate: f32,
    */
    pub args: AnalyzeArgs,

    pub ratio_method: RatioEstimationMethod,
}



impl ErrorAnalyzer {
    pub fn new(args: AnalyzeArgs) -> Self {
        let method = match args.estimation_method.as_str() {
            "slope" => RatioEstimationMethod::Slope,
            "linear_fit" => RatioEstimationMethod::LinearFit,
            "ratio_mean" => RatioEstimationMethod::RatioMean,
            "sum_ratio" => RatioEstimationMethod::SumRatio,
            _ => {
                panic!("Unknown ratio estimation method: {}. Supported methods are: slope, linear_fit, ratio_mean, sum_ratio.", args.estimation_method);
            }
        };


        Self {
            args,
            ratio_method: method,
        }
    }


    /*
    ========================
    Util functions for estimating the mean and variance of error rates
    ========================
    */

    /**
     * Perform linear regression with model y = k * x
     * return the slope k
     */
    fn slope(x: &Vec<u32>, y: &Vec<u32>, indices: &Vec<usize>) -> f32 {
        let n = indices.len() as f32;

        if n <= 1. {
            return 0.;
        }

        let sum_xy: f32 = indices.iter().map(|&i| x[i] * y[i]).sum::<u32>() as f32;
        let sum_x2: f32 = indices.iter().map(|&i| x[i] * x[i]).sum::<u32>() as f32;
        if sum_x2 == 0.0 {
            return 0.;
        }
        let k = sum_xy / sum_x2;

        k
    }

    /**
     * Perform linear regression with model y = k * x + b
     * return the slope k and intercept b
     */
    fn linear_fit(x: &Vec<u32>, y: &Vec<u32>, indices: &Vec<usize>) -> (f32, f32) {
        let n = indices.len() as f32;

        if n <= 2. {
            return (0., 0.);
        }
        let sum_x: f32 = indices.iter().map(|&i| x[i]).sum::<u32>() as f32;
        let sum_y: f32 = indices.iter().map(|&i| y[i]).sum::<u32>() as f32;
        let sum_xy: f32 = indices.iter().map(|&i| x[i] * y[i]).sum::<u32>() as f32;
        let sum_x2: f32 = indices.iter().map(|&i| x[i] * x[i]).sum::<u32>() as f32;

        let denom = n * sum_x2 - sum_x * sum_x;
        if denom == 0.0 {
            return (0., 0.);
        }

        let k = (n * sum_xy - sum_x * sum_y) / denom;
        let b = (sum_y - k * sum_x) / n;

        (k, b)
    }

    #[allow(dead_code)]
    fn linear_fit_f32(x: &Vec<f32>, y: &Vec<f32>) -> (f32, f32) {
        let n = x.len() as f32;

        if n <= 2. {
            return (0., 0.);
        }
        let sum_x: f32 = x.iter().sum::<f32>();
        let sum_y: f32 = y.iter().sum::<f32>();
        let sum_xy: f32 = x.iter().zip(y.iter()).map(|(&xi, &yi)| xi * yi).sum::<f32>();
        let sum_x2: f32 = x.iter().map(|&xi| xi * xi).sum::<f32>();

        let denom = n * sum_x2 - sum_x * sum_x;
        if denom == 0.0 {
            return (0., 0.);
        }

        let k = (n * sum_xy - sum_x * sum_y) / denom;
        let b = (sum_y - k * sum_x) / n;

        (k, b)
    }

    #[allow(dead_code)]
    fn ridge_fit_f32(x: &Vec<f32>, y: &Vec<f32>, lambda: f32) -> (f32, f32) {
        let n = x.len();
        // Means
        let sum_x = x.iter().sum::<f32>();
        let sum_y = y.iter().sum::<f32>();

        //println!("sum_x: {}, sum_y: {}", sum_x, sum_y);

        let mean_x = sum_x / n as f32;
        let mean_y = sum_y / n as f32;

        //println!("mean_x: {}, mean_y: {}", mean_x, mean_y);

        // Centered sums
        let mut sxx = 0.0f32;
        let mut sxy = 0.0f32;
        for i in 0..n {
            let dx = x[i] - mean_x;
            let dy = y[i] - mean_y;
            sxx += dx * dx;
            sxy += dx * dy;
        }

        //println!("sxx: {}, sxy: {}", sxx, sxy);

        // Ridge on slope only
        let denom = sxx + lambda;
        let k = if denom != 0.0 { sxy / denom } else { 0.0 };
        let b = mean_y - k * mean_x;

        //println!("k: {}, b: {}", k, b);

        (k, b)
    }

    fn linear_fit_huber_f32(x: &Vec<f32>, y: &Vec<f32>) -> (f32, f32) {
        let (slope, intercept) = huber_ridge_fit_1d(x, y, 0.1, 0.5, 100, 1e-6);
        
        (slope, intercept)
    }

    /**
     * Calculate the mean of the ratios y/x
     */
    fn ratio_mean(x: &Vec<u32>, y: &Vec<u32>, indices: &Vec<usize>) -> f32 {
        let n: f32 = indices.len() as f32;

        if n <= 1. {
            return 0.;
        }

        let ratio = indices.iter()
            .filter_map(|&i| if x[i] != 0 { Some(y[i] as f32 / x[i] as f32) } else { None })
            .collect::<Vec<f32>>();

        if ratio.is_empty() {
            return 0.;
        }

        let k: f32 = ratio.iter().sum::<f32>() / ratio.len() as f32;

        k
    }

    /**
     * Calculate the sum of ratios sum(y) / sum(x)
     */
    fn sum_ratio(x: &Vec<u32>, y: &Vec<u32>, indices: &Vec<usize>) -> f32 {
        let n: f32 = indices.len() as f32;

        if n == 0. {
            return 0.;
        }

        let sum_x: f32 = indices.iter().map(|&i| x[i]).sum::<u32>() as f32;
        if sum_x == 0.0 {
            return 0.;
        }
        let sum_y: f32 = indices.iter().map(|&i| y[i]).sum::<u32>() as f32;
        //println!("sum_y: {}, sum_x: {}", sum_y, sum_x);
        sum_y / sum_x
    }

    /**
     * The function that uses different methods to calculate the ratio
     */
    fn calculate_ratio(&self, x: &Vec<u32>, y: &Vec<u32>, indices: &Vec<usize>) -> f32 {
        match self.ratio_method {
            RatioEstimationMethod::Slope => Self::slope(x, y, indices),
            RatioEstimationMethod::LinearFit => {
                let (k, _) = Self::linear_fit(x, y, indices);
                k
            },
            RatioEstimationMethod::RatioMean => Self::ratio_mean(x, y, indices),
            RatioEstimationMethod::SumRatio => Self::sum_ratio(x, y, indices),
        }
    }

    fn sum_indices(&self, x: &Vec<u32>, indices: &Vec<usize>) -> u32 {
        indices.iter().map(|&i| x[i]).sum()
    }

    fn random_subsample_with_replacement(x: &Vec<usize>, n: usize) -> Vec<usize> {
        let mut rng = rand::rng();
        (0..n)
            .map(|_| x[rng.random_range(0..x.len())])
            .collect()
    }

    #[allow(dead_code)]
    fn random_subsample_without_replacement(x: &Vec<usize>, n: usize) -> Vec<usize> {
        use rand::seq::SliceRandom;
        let mut rng = rand::rng();
        let mut x_clone = x.clone();
        x_clone.shuffle(&mut rng);
        x_clone.truncate(n);
        x_clone
    }


    /*
    ========================
    Util functions for estimating the parameters of beta distribution
    ========================
    */
    /*
    fn residual_hazard_ratio_beta_distribution(params: &[f64], data: &[f64]) -> Array1<f64> {
        let n = data.len() / 2;
        let x = &data[0..n];
        let y = &data[n..];

        let mut res = Array1::zeros(n);
        for i in 0..n {
            res[i] = y[i] - (params[0] / (params[0] + params[1] + x[i]));
        }
        res
    }

    fn residual_hazard_ratio_beta_distribution_fixed_kappa(params: &[f64], data: &[f64]) -> Array1<f64> {
        let n = data.len() / 2;
        let x = &data[0..n];
        let y = &data[n..];

        let mut res = Array1::zeros(n);
        for i in 0..n {
            res[i] = y[i] - (params[0] / (MIN_KAPPA as f64 + x[i]));
        }
        res
    }

    fn fit_hazard_ratio_beta_distribution(&self, hazard_ratios: &Vec<f32>) -> (f32, f32) {
        let n = hazard_ratios.len();
        let mut vec_data: Vec<f64> = Vec::with_capacity(n * 2);
        for i in 1..=n {
            vec_data.push(i as f64 + self.args.k as f64);
        }
        for &hr in hazard_ratios.iter() {
            vec_data.push(hr as f64);
        }

        // convert the Vec<f64> to an Array1<f64> as required by robust_least_squares
        let data = Array1::from_vec(vec_data);

        let initial_params = array![1.0f64, 1.0f64];
        let result = least_squares(
            &Self::residual_hazard_ratio_beta_distribution,
            &initial_params,
            Method::LevenbergMarquardt,
            None::<fn(&[f64], &[f64]) -> scirs2_core::ndarray::Array2<f64>>, 
            &data, None
        ).expect("robust_least_squares failed");

        if result.x[1] <= MIN_KAPPA as f64 {
            // refit with fixed kappa
            let initial_params_fixed = array![1.0f64];
            //warn!("Estimated beta parameter is too small ({}), probably due to high error rate and low coverage. Refitting with fixed alpha + beta = {}.", result.x[1], MIN_KAPPA);
            //warn!("Consider increasing the coverage or using bidirectional kmers to improve the estimation.");
            let result_fixed = least_squares(
                &Self::residual_hazard_ratio_beta_distribution_fixed_kappa,
                &initial_params_fixed,
                Method::LevenbergMarquardt,
                None::<fn(&[f64], &[f64]) -> scirs2_core::ndarray::Array2<f64>>, 
                &data, None
            ).expect("robust_least_squares failed");

            return (result_fixed.x[0] as f32, MIN_KAPPA - result_fixed.x[0] as f32);
        }

        //println!("data: {:?}", data);

        (result.x[0] as f32, result.x[1] as f32)
    }
    */

    /*
    ========================
    Util functions for estimating the parameters of Weibull distribution
    ========================
    */

    /// Assume hazard ratio follows Weibull distribution: hazard ratio = a * (i + k)^b
    #[allow(dead_code)]
    fn fit_hazard_ratio_weibull_distribution_power_law(&self, hazard_ratios: &Vec<f32>) -> (f32, f32) {
        // Fit hazard ratio = a * (i + k)^b, or log(hazard ratio) = log(a) + b * log(i + k)
        let x = hazard_ratios.iter().enumerate().
            map(|(i, _)| (i as f32 + self.args.k as f32).ln())
            .collect::<Vec<f32>>();
        let y = hazard_ratios.iter()
            .map(|&hr| if hr > 0.0 { hr.ln() } else { 0.0 })
            .collect::<Vec<f32>>();
        //let (b, log_a) = Self::linear_fit_f32(&x, &y);
        let (b, log_a) = Self::ridge_fit_f32(&x, &y, 1.);
        
        
        let a = log_a.exp();

        (a, b)
    }


    /// Assume hazard ratio follows discrete Weibull distribution 
    /// h(t) = 1 - exp(-lambda * ((t+1)^beta - t^beta))
    /// By approximation, we can fit log(-log(1 - hazard ratio)) \approx log(lambda) + beta * log(t)
    fn fit_hazard_ratio_weibull_distribution_cloglog(&self, hazard_ratios: &Vec<f32>) -> (f32, f32) {
        // Fit hazard ratio = a * (i + k)^b, or log(hazard ratio) = log(a) + b * log(i + k)
        let x = hazard_ratios.iter().enumerate().
            map(|(i, _)| (i as f32 + self.args.k as f32).ln())
            .collect::<Vec<f32>>();
        // complementary log-log, clip hazard ratios to avoid log(0)
        let y = hazard_ratios.iter()
            .map(|&hr|
                (-(- hr.clamp(EPSILON, 1.0 - EPSILON)).ln_1p()).ln())
            .collect::<Vec<f32>>();
        //let (b, log_a) = Self::linear_fit_f32(&x, &y);
        //let (slope, intercept) = Self::ridge_fit_f32(&x, &y, 1.);
        let (slope, intercept) = Self::linear_fit_huber_f32(&x, &y);
        
        
        let beta = slope + 1.;
        let lambda = intercept.exp() / beta;

        (lambda, beta)
    }

    fn fit_hazard_ratio_constant(&self, hazard_ratios: &Vec<f32>) -> f32 {
        let n = hazard_ratios.len();
        if n == 0 {
            return 0.;
        }
        let mean = hazard_ratios.iter().sum::<f32>() / n as f32;
        mean
    }

    fn fit_hazard_ratio(&self, hazard_ratios: &Vec<f32>) -> (f32, f32) {
        match self.args.hazard_model.as_str() {
            "weibull" => self.fit_hazard_ratio_weibull_distribution_cloglog(hazard_ratios),
            "constant" => (self.fit_hazard_ratio_constant(hazard_ratios), 1.0),
            _ => {
                panic!("Unknown hazard model: {}. Supported models are: weibull, constant.", self.args.hazard_model);
            }
        }
    }

    /**
     * Estimate 1/E[T], where T~DiscreteWeibull(lambda, beta)
     * E[T] = sum_{t=1}^{\infty} P(T >= t) = sum_{t=1}^{\infty} exp(-lambda * t^beta)
     * approximate the sum until the terms are small enough
     */
    fn estimate_effective_error_rate(&self, lambda: f32, beta: f32) -> f32 {
        let mut expected_t = 0.0;
        let epsilon: f32 = 1e-6;
        let max_iterations: usize = 10000;
        for t in 1..max_iterations {
            let survival_prob = (- lambda * (t as f32).powf(beta)).exp();
            if survival_prob < epsilon {
                break;
            }
            expected_t += survival_prob;
        }
        1.0 / expected_t
    }


    

    /// Approximate erfc(x) using a rational polynomial (Abramowitz & Stegun 7.1.26, max error 1.5e-7).
    fn erfc_approx(x: f64) -> f64 {
        if x < 0.0 {
            return 2.0 - Self::erfc_approx(-x);
        }
        let t = 1.0 / (1.0 + 0.3275911 * x);
        let poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
        poly * (-x * x).exp()
    }

    fn normal_cdf(z: f64) -> f64 {
        0.5 * Self::erfc_approx(-z / std::f64::consts::SQRT_2)
    }

    /// P(X <= k_obs) where X ~ Binomial(n, p), using normal approximation with continuity correction.
    fn binomial_cdf_lower(n: u32, k_obs: u32, p: f64) -> f64 {
        if n == 0 { return 1.0; }
        if p <= 0.0 { return 1.0; }
        if p >= 1.0 { return if k_obs >= n { 1.0 } else { 0.0 }; }
        let mean = n as f64 * p;
        let var = n as f64 * p * (1.0 - p);
        if var < 1e-10 { return if k_obs as f64 >= mean { 1.0 } else { 0.0 }; }
        let z = (k_obs as f64 + 0.5 - mean) / var.sqrt();
        Self::normal_cdf(z)
    }

    /**
     * Identify outliers iteratively using the Weibull hazard model.
     *
     * Algorithm:
     *  1. Estimate lambda and beta from all active keys.
     *  2. For each v and each active key, compute h(t) = 1 - exp(-lambda*(t^beta - (t-1)^beta)).
     *     Treat the count at t as Binomial(n, 1-h(t)) where n is the count at t-1.
     *     If P(X <= observed) < outlier_threshold (default 1e-9), mark the key as an outlier.
     *  3. Re-estimate lambda and beta from the remaining keys.
     *  4. Repeat until lambda and beta change by less than 1e-4.
     */
    pub fn find_hazard_ratio_outliers(&self, stats: &KVmerStats) -> Vec<usize> {
        let n_keys = stats.error_summary.consensus_counts.len();
        let mut active = vec![true; n_keys];

        let max_iter = 10;
        let convergence_tol = 1e-5_f32;
        let p_threshold = self.args.outlier_threshold as f64;

        let mut prev_lambda = f32::INFINITY;
        let mut prev_beta = f32::INFINITY;

        let v_max = stats.v - self.args.ignore_last_hazard_ratios as u8;

        for iter in 0..max_iter {
            let indices: Vec<usize> = (0..n_keys).filter(|&i| active[i]).collect();
            if indices.is_empty() { break; }

            let (lambda, beta, _, _, _) = self.estimate_hazard_ratio(stats, &indices);

            if (lambda - prev_lambda).abs() < convergence_tol && (beta - prev_beta).abs() < convergence_tol {
                info!("Iterative outlier removal converged after {} iteration(s).", iter);
                break;
            }
            prev_lambda = lambda;
            prev_beta = beta;

            for v in 1..=v_max {
                // t is the time coordinate used when fitting: t = (v-1) + k
                let t = (v - 1) as f64 + self.args.k as f64;
                let t_beta = t.powf(beta as f64);
                let t1_beta = if t > 1.0 { (t - 1.0).powf(beta as f64) } else { 0.0 };
                let h_t = 1.0 - (-(lambda as f64) * (t_beta - t1_beta)).exp();
                let p_survival = (1.0 - h_t).clamp(0.0, 1.0);

                for i in 0..n_keys {
                    if !active[i] { continue; }

                    let n = if v == 1 {
                        stats.error_summary.total_counts[i]
                    } else {
                        stats.error_summary.consensus_up_to_v_counts[(v - 2) as usize][i]
                    };
                    let k_obs = stats.error_summary.consensus_up_to_v_counts[(v - 1) as usize][i];

                    if n == 0 { continue; }

                    if Self::binomial_cdf_lower(n, k_obs, p_survival) < p_threshold {
                        active[i] = false;
                    }
                }
            }
        }

        let indices: Vec<usize> = (0..n_keys).filter(|&i| active[i]).collect();
        info!("Identified {} inliers out of {} data points based on iterative Binomial outlier removal ({}%).",
            indices.len(), n_keys, (indices.len() as f32 / n_keys as f32) * 100.0);
        indices
    }


    /**
     * Estimate the error rate and 5-95% confidence interval using bootstrap
     * for all the error types
     */
    pub fn estimate_error_rate(&self, stats: &KVmerStats, indices: &Vec<usize>) -> HashMap<(EditOperation, u8, u8), u32> {
        // initialize the error count arrays
        let mut error_counts: HashMap<(EditOperation, u8, u8), u32> = HashMap::new();

        indices.iter().for_each(|&i| {
            for (ni, count_map) in stats.error_spectrum.error_counts[i].iter() {
                let count = error_counts.entry((ni.op, ni.prev_base, ni.next_base)).or_insert(0);
                *count += *count_map;
            }
        });
        
        // calculate the mean for each error type using the full error_counts vector
        let mut estimates: HashMap<(EditOperation, u8, u8), u32> = HashMap::new();
        for op in ALL_OPERATIONS.iter() {
            for prev_base in 0..4 {
                for next_base in 0..4 {
                    let count = error_counts.get(&(*op, prev_base, next_base)).unwrap_or(&0);
                    estimates.insert((*op, prev_base, next_base), *count);
                }
            }
        }
        
        /*
        // bootstrap to estimate the 5-95% confidence interval
        let mut bootstrap_estimates: HashMap<EditOperation, Vec<f32>> = HashMap::new();
        for op in operations.iter() {
            bootstrap_estimates.insert(*op, Vec::new());
        }

        for _ in 0..self.num_experiments {
            let indices_sample = Self::random_subsample_with_replacement(indices, indices.len() as usize);
            
            let mut error_count: HashMap<EditOperation, u32> = HashMap::new();
            for op in operations.iter() {
                let sum = self.sum_indices(&error_counts[op], &indices_sample);
                error_count.insert(*op, sum);
            }
            // normalize the error counts so that they sum to 1
            let total_count: u32 = error_count.values().sum();
            for op in operations.iter() {
                let rate = if total_count > 0 {
                    error_count[op] as f32 / total_count as f32
                } else {
                    0.0
                };
                bootstrap_estimates.get_mut(op).unwrap().push(rate);
            }
        }

        // calculate the 5-95% confidence interval
        let mut result: Vec<(f32, (f32, f32))> = Vec::new();
        for op in operations.iter() {
            let mut estimates_op = bootstrap_estimates[op].clone();
            estimates_op.sort_by(f32::total_cmp);
            let n = estimates_op.len();
            let mean = estimates[op];
            let lower = estimates_op[(n as f32 * 0.05) as usize];
            let upper = estimates_op[(n as f32 * 0.95) as usize];
            result.push((mean, (lower, upper)));
        }
        */

        estimates
    }

    /**
     * Returns ((lower_lambda, upper_lambda), (lower_beta, upper_beta), hazard_ratio_list, (lower_error_rate, upper_error_rate))
     */
    pub fn estimate_hazard_ratio_confidence_interval(&self, stats: &KVmerStats, indices: &Vec<usize>) -> ((f32, f32), (f32, f32), Vec<(f32, f32)>, (f32, f32)) {
        let mut x: &Vec<u32>;
        let mut y: &Vec<u32>;

        // record the estimated a and b
        let mut lambda_list: Vec<f32> = Vec::new();
        let mut beta_list: Vec<f32> = Vec::new();
        let mut error_rate_list: Vec<f32> = Vec::new();

        // record hazard ratios for each v
        let mut hazard_ratio_list: Vec<Vec<f32>> = Vec::new();
        for _v in 1..=(stats.v - self.args.ignore_last_hazard_ratios as u8) {
            hazard_ratio_list.push(Vec::new());
        }

        for _ in 0..self.args.num_experiments {
            let indices_sample = Self::random_subsample_with_replacement(indices, indices.len() as usize);

            let mut hazard_ratios: Vec<f32> = Vec::new();

            for v in 1..=(stats.v - self.args.ignore_last_hazard_ratios as u8) {
                if v - 1 == 0 {
                    x = &stats.error_summary.total_counts;
                    y = &stats.error_summary.consensus_up_to_v_counts[0];
                } else {
                    x = &stats.error_summary.consensus_up_to_v_counts[(v - 1 - 1) as usize];
                    y = &stats.error_summary.consensus_up_to_v_counts[(v - 1) as usize];
                }

                let h = self.calculate_ratio(x, y, &indices_sample);
                hazard_ratios.push(1. - h);
                hazard_ratio_list[(v - 1) as usize].push(1. - h);
            }
            // estimate the parameters of the beta distribution
            //let (alpha, beta) = self.fit_hazard_ratio_beta_distribution(&hazard_ratios, (indices.len() as f32 * self.bootstrap_sample_rate) as usize);
            let (lambda, beta) = self.fit_hazard_ratio(&hazard_ratios);
            lambda_list.push(lambda);
            beta_list.push(beta);
            error_rate_list.push(self.estimate_effective_error_rate(lambda, beta));
        }

        lambda_list.sort_by(f32::total_cmp);
        let lower_lambda = lambda_list[(self.args.num_experiments as f32 * 0.05) as usize];
        let upper_lambda = lambda_list[(self.args.num_experiments as f32 * 0.95) as usize];

        beta_list.sort_by(f32::total_cmp);
        let lower_beta = beta_list[(self.args.num_experiments as f32 * 0.05) as usize];
        let upper_beta = beta_list[(self.args.num_experiments as f32 * 0.95) as usize];

        let mut hazard_ratio_range_list: Vec<(f32, f32)> = Vec::new();
        for v in 0..hazard_ratio_list.len() {
            hazard_ratio_list[v].sort_by(f32::total_cmp);
            let h_lower = hazard_ratio_list[v][(self.args.num_experiments as f32 * 0.05) as usize];
            let h_upper = hazard_ratio_list[v][(self.args.num_experiments as f32 * 0.95) as usize];
            hazard_ratio_range_list.push((h_lower, h_upper));
        }

        error_rate_list.sort_by(f32::total_cmp);
        let lower_error_rate = error_rate_list[(self.args.num_experiments as f32 * 0.05) as usize];
        let upper_error_rate = error_rate_list[(self.args.num_experiments as f32 * 0.95) as usize];

        ((lower_lambda, upper_lambda), (lower_beta, upper_beta), hazard_ratio_range_list, (lower_error_rate, upper_error_rate))
    }


    // returns (estimated_lambda, estimated_beta, hazard_ratios, x_sum, y_sum)
    pub fn estimate_hazard_ratio(&self, stats: &KVmerStats, indices: &Vec<usize>) -> (f32, f32, Vec<f32>, Vec<u32>, Vec<u32>) {
        let mut x: &Vec<u32>;
        let mut y: &Vec<u32>;

        let mut hazard_ratios: Vec<f32> = Vec::new();
        let mut x_sum: Vec<u32> = Vec::new();
        let mut y_sum: Vec<u32> = Vec::new();

        for v in 1..=(stats.v - self.args.ignore_last_hazard_ratios as u8) {
            if v - 1 == 0 {
                x = &stats.error_summary.total_counts;
                y = &stats.error_summary.consensus_up_to_v_counts[0];
            } else {
                x = &stats.error_summary.consensus_up_to_v_counts[(v - 1 - 1) as usize];
                y = &stats.error_summary.consensus_up_to_v_counts[(v - 1) as usize];
            }

            let h = self.calculate_ratio(x, y, indices);
            hazard_ratios.push(1. - h);
            x_sum.push(self.sum_indices(x, indices));
            y_sum.push(self.sum_indices(y, indices));
        }

        let (lambda, beta) = self.fit_hazard_ratio(&hazard_ratios);
        //println!("Weibull parameters: alpha = {}, beta = {}", a, b);
        (lambda, beta, hazard_ratios, x_sum, y_sum) 
    }

    pub fn key_coverage(&self, stats: &KVmerStats, indices: &Vec<usize>) -> (f32, (f32, f32)) {
        let mut coverages: Vec<u32> = indices.iter().map(|&i| stats.error_summary.total_counts[i]).collect();
        coverages.sort_unstable();
        let n = coverages.len();
        if n == 0 {
            return (0., (0., 0.));
        }

        let median_coverage = if n % 2 == 0 {
            (coverages[n / 2 - 1] + coverages[n / 2]) as f32 / 2.0
        } else {
            coverages[n / 2] as f32
        };
        let coverage_ci_lower = coverages[(n as f32 * 0.05) as usize] as f32;
        let coverage_ci_upper = coverages[(n as f32 * 0.95) as usize] as f32;
        
        (median_coverage, (coverage_ci_lower, coverage_ci_upper))
    }

    pub fn estimate_true_coverage(&self, estimated_lambda: f32, estimated_beta: f32, key_coverage: (f32, (f32, f32))) -> (f32, (f32, f32)) {
        // estimate survival rate at k
        let survival_rate: f32 = (- estimated_lambda * ((self.args.k as f32).powf(estimated_beta))).exp();
        if survival_rate <= 0.0 || survival_rate > 1.0 {
            return (0., (0., 0.));
        }
        let mut estimated_coverage = key_coverage.0 / survival_rate;
        let mut estimated_coverage_ci_lower = (key_coverage.1).0 / survival_rate;
        let mut estimated_coverage_ci_upper = (key_coverage.1).1 / survival_rate;

        if self.args.forward_only {
            estimated_coverage *= 2.0;
            estimated_coverage_ci_lower *= 2.0;
            estimated_coverage_ci_upper *= 2.0;
        }

        (estimated_coverage, (estimated_coverage_ci_lower, estimated_coverage_ci_upper))
    }


    /// For each observed Phred score, compute the empirical error rate using only
    /// inlier keys (filtered by `find_hazard_ratio_outliers`), then bootstrap over
    /// those key indices to estimate the 5th–95th percentile confidence interval.
    ///
    /// Returns a vector of `(qscore, num_correct, num_error, error_rate, ci_lower, ci_upper)`.
    pub fn calibrate_qscores(&self, stats: &KVmerStats) -> Vec<QscoreCalibration> {
        // Filter outlier keys
        let indices = if !self.args.use_all {
            self.find_hazard_ratio_outliers(stats)
        } else {
            (0..stats.error_summary.consensus_counts.len()).collect()
        };

        // Aggregate qscore counts from inlier keys
        let mut qscore_correct: HashMap<u8, u64> = HashMap::new();
        let mut qscore_error: HashMap<u8, u64> = HashMap::new();
        for &i in &indices {
            for (&q, &c) in &stats.phred_summary.correct_per_key[i] {
                *qscore_correct.entry(q).or_insert(0) += c;
            }
            for (&q, &e) in &stats.phred_summary.error_per_key[i] {
                *qscore_error.entry(q).or_insert(0) += e;
            }
        }

        let mut qscores: Vec<u8> = qscore_correct.keys()
            .chain(qscore_error.keys())
            .cloned()
            .collect();
        qscores.sort_unstable();
        qscores.dedup();

        // Bootstrap: resample key indices with replacement, recompute error rates
        let mut bootstrap_rates: HashMap<u8, Vec<f64>> = qscores.iter().map(|&q| (q, Vec::new())).collect();
        for _ in 0..self.args.num_experiments {
            let sample = Self::random_subsample_with_replacement(&indices, indices.len());
            let mut c_sample: HashMap<u8, u64> = HashMap::new();
            let mut e_sample: HashMap<u8, u64> = HashMap::new();
            for &i in &sample {
                for (&q, &c) in &stats.phred_summary.correct_per_key[i] {
                    *c_sample.entry(q).or_insert(0) += c;
                }
                for (&q, &e) in &stats.phred_summary.error_per_key[i] {
                    *e_sample.entry(q).or_insert(0) += e;
                }
            }
            for &q in &qscores {
                let c = *c_sample.get(&q).unwrap_or(&0);
                let e = *e_sample.get(&q).unwrap_or(&0);
                let total = c + e;
                let rate = if total > 0 { e as f64 / total as f64 } else { 0.0 };
                bootstrap_rates.get_mut(&q).unwrap().push(rate);
            }
        }

        // Build result with point estimates and CI
        let mut result = Vec::new();
        for &q in &qscores {
            let correct = *qscore_correct.get(&q).unwrap_or(&0);
            let error   = *qscore_error.get(&q).unwrap_or(&0);
            let total   = correct + error;
            let error_rate = if total > 0 { error as f64 / total as f64 } else { 0.0 };

            let mut rates = bootstrap_rates[&q].clone();
            rates.sort_by(f64::total_cmp);
            let n = rates.len();
            let lower = if n > 0 { rates[(n as f64 * 0.05) as usize] } else { 0.0 };
            let upper = if n > 0 { rates[((n as f64 * 0.95) as usize).min(n - 1)] } else { 0.0 };

            result.push(QscoreCalibration { qscore: q, num_correct: correct, num_error: error, error_rate, ci_lower: lower, ci_upper: upper });
        }
        result
    }

    /// Like `calibrate_qscores`, but aggregates correct/error counts by read position
    /// (from start and from end) across inlier keys.  No bootstrap CI is computed
    /// since the output does not include confidence intervals.
    pub fn calibrate_read_positions(&self, stats: &KVmerStats) -> Vec<ReadPositionCalibration> {
        let indices = if !self.args.use_all {
            self.find_hazard_ratio_outliers(stats)
        } else {
            (0..stats.error_summary.consensus_counts.len()).collect()
        };

        let mut correct_from_start: HashMap<u32, u64> = HashMap::new();
        let mut correct_from_end: HashMap<u32, u64> = HashMap::new();
        let mut error_from_start: HashMap<u32, u64> = HashMap::new();
        let mut error_from_end: HashMap<u32, u64> = HashMap::new();

        for &i in &indices {
            for (&pos, &c) in &stats.read_position_summary.correct_from_start_per_key[i] {
                *correct_from_start.entry(pos).or_insert(0) += c;
            }
            for (&pos, &c) in &stats.read_position_summary.correct_from_end_per_key[i] {
                *correct_from_end.entry(pos).or_insert(0) += c;
            }
            for (&pos, &e) in &stats.read_position_summary.error_from_start_per_key[i] {
                *error_from_start.entry(pos).or_insert(0) += e;
            }
            for (&pos, &e) in &stats.read_position_summary.error_from_end_per_key[i] {
                *error_from_end.entry(pos).or_insert(0) += e;
            }
        }

        let mut result = Vec::new();

        let mut start_positions: Vec<u32> = correct_from_start.keys().chain(error_from_start.keys()).copied().collect();
        start_positions.sort_unstable();
        start_positions.dedup();
        for pos in start_positions {
            result.push(ReadPositionCalibration {
                index: pos,
                from_start: true,
                num_correct: *correct_from_start.get(&pos).unwrap_or(&0),
                num_error: *error_from_start.get(&pos).unwrap_or(&0),
            });
        }

        let mut end_positions: Vec<u32> = correct_from_end.keys().chain(error_from_end.keys()).copied().collect();
        end_positions.sort_unstable();
        end_positions.dedup();
        for pos in end_positions {
            result.push(ReadPositionCalibration {
                index: pos,
                from_start: false,
                num_correct: *correct_from_end.get(&pos).unwrap_or(&0),
                num_error: *error_from_end.get(&pos).unwrap_or(&0),
            });
        }

        result
    }

    pub fn analyze(&self, stats: &KVmerStats) -> ErrorSpectrum {
        // exclude the hazard ratio outliers
        let indices = if !self.args.use_all {
            self.find_hazard_ratio_outliers(stats)
        } else {
            (0..stats.error_summary.consensus_counts.len()).collect()
        };

        // estimate SNP rates
        let error_rates = self.estimate_error_rate(stats, &indices);

        // estimate hazard ratio parameters
        let (lambda, beta, hazard_ratio, x_sum, y_sum) = self.estimate_hazard_ratio(stats, &indices);
        let (lambda_ci, beta_ci, hazard_ratio_ci, error_rate_ci) = self.estimate_hazard_ratio_confidence_interval(stats, &indices);
        let effective_error_rate = self.estimate_effective_error_rate(lambda, beta);
        let per_base_error_rate = 1.0 - (-lambda).exp();
        let per_base_error_rate_ci = (1.0 - (-(lambda_ci.0)).exp(), 1.0 - (-(lambda_ci.1)).exp());

        if let Some(prefix) = &self.args.output_prefix {
            use std::fs;
            use std::fs::File;
            use std::io::{BufWriter, Write};

            let hazard_ratio_output = format!("{}.hazard_rate.csv", prefix);
            let file = File::create(&hazard_ratio_output).expect("Could not create hazard ratio output file.");
            let mut writer = BufWriter::new(file);

            writeln!(writer, "t,num_candidates,num_survival,hazard_ratio,5th_percentile,95th_percentile").expect("Could not write to hazard ratio output file.");
            for v in 0..hazard_ratio.len() {
                writeln!(writer, "{},{},{},{:.6},{:.6},{:.6}",
                    v + 1 + self.args.k as usize,
                    x_sum[v],
                    y_sum[v],
                    hazard_ratio[v],
                    hazard_ratio_ci[v].0,
                    hazard_ratio_ci[v].1
                ).expect("Could not write to hazard ratio output file.");
            }

            fs::write(format!("{}.kvmer.csv", prefix), stats.error_summary.to_csv(Some(&indices))).unwrap();
            fs::write(format!("{}.summary_error_spectrum.csv", prefix), stats.error_spectrum.to_csv(Some(&indices))).unwrap();
            fs::write(format!("{}.summary_error_spectrum_dependence_on_t.csv", prefix), stats.error_spectrum.to_dependence_on_t_csv(Some(&indices), self.args.k as usize, self.args.ignore_last_hazard_ratios)).unwrap();
            fs::write(format!("{}.summary_phred.csv", prefix), stats.phred_summary.to_csv(Some(&indices))).unwrap();
            fs::write(format!("{}.summary_read_position.csv", prefix), stats.read_position_summary.to_csv(Some(&indices))).unwrap();
        }

        // estimate key coverage
        let key_coverage = self.key_coverage(stats, &indices);
        let estimated_coverage = self.estimate_true_coverage(lambda, beta, key_coverage);

        ErrorSpectrum {
            estimated_lambda: (lambda, lambda_ci),
            estimated_beta: (beta, beta_ci),
            per_base_error_rate: (per_base_error_rate, per_base_error_rate_ci),
            effective_error_rate: (effective_error_rate, error_rate_ci),

            key_coverage: key_coverage,
            estimated_coverage: estimated_coverage,

            snp_rate: error_rates,
            bidirectional: !self.args.forward_only,
        }
    }
}


/**
 * Format the spectrum into a line in a csv file
 */
pub fn spectrum_to_str(spectrum: &ErrorSpectrum, bidirectional: bool) -> String {
    let mut result = String::new();

    if bidirectional != spectrum.bidirectional {
        panic!("The bidirectional flag does not match the spectrum data.");
    }

    // per-base and effective error rate
    result.push_str(&format!("{:.6},{:.6}~{:.6},", spectrum.per_base_error_rate.0, (spectrum.per_base_error_rate.1).0, (spectrum.per_base_error_rate.1).1));
    result.push_str(&format!("{:.6},{:.6}~{:.6},", spectrum.effective_error_rate.0, (spectrum.effective_error_rate.1).0, (spectrum.effective_error_rate.1).1));

    // hazard ratio parameters a and b
    result.push_str(&format!("{:.6},{:.6}~{:.6},", spectrum.estimated_lambda.0, (spectrum.estimated_lambda.1).0, (spectrum.estimated_lambda.1).1));
    result.push_str(&format!("{:.6},{:.6}~{:.6},", spectrum.estimated_beta.0, (spectrum.estimated_beta.1).0, (spectrum.estimated_beta.1).1));

    // key coverage and estimated true coverage
    result.push_str(&format!("{:.6},{:.6}~{:.6},", spectrum.key_coverage.0, (spectrum.key_coverage.1).0, (spectrum.key_coverage.1).1));
    result.push_str(&format!("{:.6},{:.6}~{:.6},", spectrum.estimated_coverage.0, (spectrum.estimated_coverage.1).0, (spectrum.estimated_coverage.1).1));

    // remove the last comma
    result.pop();

    result
}


pub fn header_str(bidirectional: bool) -> String {
    let mut result = String::new();

    result.push_str("per_base_error_rate,per_base_error_rate_5-95th_percentile,");
    result.push_str("effective_error_rate,effective_error_rate_5-95th_percentile,");
    result.push_str("lambda,lambda_5-95th_percentile,");
    result.push_str("beta,beta_5-95th_percentile,");

    result.push_str("key_median_coverage,key_coverage_5-95th_percentile,");
    result.push_str("true_median_coverage,true_coverage_5-95th_percentile");

    result
}