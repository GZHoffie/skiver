/// Fits a Huber regression model (1D: single feature + intercept) using IRLS,
/// matching sklearn's HuberRegressor with its default parameters:
///   epsilon=1.35, alpha=0.0001, max_iter=100, tol=1e-5, fit_intercept=True
///
/// The objective minimised is (following sklearn's formulation):
///   L(slope, intercept, sigma)
///     = 0.5 * Σ_inliers  r_i² / sigma
///     + 2·epsilon · Σ_outliers |r_i|
///     − epsilon² · sigma · n_outliers
///     + n · sigma
///     + alpha · slope²
///
/// where r_i = y_i − slope·x_i − intercept and the inlier set is
/// { i : |r_i| ≤ epsilon · sigma }.
///
/// Returns (slope, intercept).
pub fn huber_regression_1d(
    x: &[f32],
    y: &[f32],
    epsilon: f32,    // sklearn default: 1.35
    alpha: f32,      // L2 penalty on slope (not intercept); sklearn default: 0.0001
    max_iter: usize, // sklearn default: 100
    tol: f32,        // sklearn default: 1e-5
) -> (f32, f32) {
    let n = x.len();
    if n < 2 || n != y.len() || epsilon < 1.0 || alpha < 0.0 || max_iter == 0 || tol <= 0.0 {
        return (0.0, 0.0);
    }

    // Solve weighted ridge normal equations (ridge on slope only):
    //   (sxx + ridge) · slope + sx · intercept = sxy
    //   sx · slope  + sw · intercept            = sy
    // Accumulate in f64 for numerical stability.
    fn solve_weighted(x: &[f32], y: &[f32], w: &[f32], ridge: f64) -> Option<(f32, f32)> {
        let mut sw = 0.0f64;
        let mut sx = 0.0f64;
        let mut sxx = 0.0f64;
        let mut sy = 0.0f64;
        let mut sxy = 0.0f64;
        for i in 0..x.len() {
            let wi = w[i] as f64;
            let xi = x[i] as f64;
            let yi = y[i] as f64;
            sw += wi;
            sx += wi * xi;
            sxx += wi * xi * xi;
            sy += wi * yi;
            sxy += wi * xi * yi;
        }
        let det = (sxx + ridge) * sw - sx * sx;
        if !det.is_finite() || det.abs() < 1e-12 {
            return None;
        }
        let slope = (sxy * sw - sy * sx) / det;
        let intercept = ((sxx + ridge) * sy - sx * sxy) / det;
        if slope.is_finite() && intercept.is_finite() {
            Some((slope as f32, intercept as f32))
        } else {
            None
        }
    }

    // Initialise with OLS (faster convergence than sklearn's w=0,b=0 start,
    // but converges to the same convex minimum).
    let unit_w = vec![1.0f32; n];
    let (mut slope, mut intercept) = match solve_weighted(x, y, &unit_w, 0.0) {
        Some(v) => v,
        None => return (0.0, 0.0),
    };
    // sklearn initialises sigma=1.0; we do the same.
    let mut sigma = 1.0f32;

    // The sklearn gradient of alpha·slope² is 2·alpha·slope, so the ridge term
    // added to the (slope,slope) entry of the normal equations is 2·alpha.
    let ridge = 2.0 * alpha as f64;
    let mut weights = vec![0.0f32; n];

    for _ in 0..max_iter {
        let mut inlier_sq_sum = 0.0f32;
        let mut n_outliers = 0usize;

        for i in 0..n {
            let r = y[i] - slope * x[i] - intercept;
            let ar = r.abs();
            if ar <= epsilon * sigma {
                inlier_sq_sum += r * r;
                weights[i] = 1.0 / sigma;
            } else {
                n_outliers += 1;
                weights[i] = (2.0 * epsilon) / ar;
            }
        }

        // Closed-form sigma update: d_L/d_sigma = 0
        //   sigma² = 0.5 · Σ_inliers r² / (n − epsilon² · n_outliers)
        let denom = n as f32 - epsilon * epsilon * n_outliers as f32;
        let new_sigma = if denom > 1e-10 {
            (0.5 * inlier_sq_sum / denom).sqrt().max(1e-6)
        } else {
            sigma // degenerate case: keep current sigma
        };

        let (new_slope, new_intercept) = match solve_weighted(x, y, &weights, ridge) {
            Some(v) => v,
            None => return (slope, intercept),
        };

        let d_slope = (new_slope - slope).abs();
        let d_intercept = (new_intercept - intercept).abs();
        let d_sigma = (new_sigma - sigma).abs();

        slope = new_slope;
        intercept = new_intercept;
        sigma = new_sigma;

        if d_slope.max(d_intercept).max(d_sigma) < tol {
            break;
        }
    }

    (slope, intercept)
}
