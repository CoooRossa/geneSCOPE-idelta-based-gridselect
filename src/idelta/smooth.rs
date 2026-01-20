pub fn gaussian_smooth_logx(widths_um: &[u32], y: &[f64], sigma: f64) -> Vec<f64> {
    assert_eq!(widths_um.len(), y.len());
    let n = widths_um.len();
    if n == 0 {
        return Vec::new();
    }

    let x: Vec<f64> = widths_um.iter().map(|&w| (w as f64).ln()).collect();

    let mut out = vec![f64::NAN; n];
    for i in 0..n {
        if !y[i].is_finite() {
            continue;
        }
        let xi = x[i];
        let mut num = 0.0;
        let mut den = 0.0;
        for j in 0..n {
            let yj = y[j];
            if !yj.is_finite() {
                continue;
            }
            let dx = xi - x[j];
            let w = (-0.5 * (dx * dx) / (sigma * sigma)).exp();
            num += w * yj;
            den += w;
        }
        if den > 0.0 {
            out[i] = num / den;
        }
    }
    out
}
