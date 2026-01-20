pub fn knee_max_distance_chord(widths_um: &[u32], y: &[f64]) -> Option<u32> {
    if widths_um.len() != y.len() || widths_um.len() < 3 {
        return None;
    }

    // Use log(width) as x.
    let mut x: Vec<f64> = widths_um.iter().map(|&w| (w as f64).ln()).collect();
    let mut yv: Vec<f64> = y.to_vec();

    // Filter to finite points.
    let mut idx: Vec<usize> = (0..x.len()).filter(|&i| yv[i].is_finite()).collect();
    if idx.len() < 3 {
        return None;
    }
    idx.sort_unstable_by_key(|&i| widths_um[i]);
    x = idx.iter().map(|&i| x[i]).collect();
    yv = idx.iter().map(|&i| yv[i]).collect();
    let w_sorted: Vec<u32> = idx.iter().map(|&i| widths_um[i]).collect();

    // Normalize x and y to [0,1].
    let (x0, x1) = (x[0], x[x.len() - 1]);
    let (y0, y1) = {
        let mut miny = f64::INFINITY;
        let mut maxy = f64::NEG_INFINITY;
        for &v in &yv {
            miny = miny.min(v);
            maxy = maxy.max(v);
        }
        (miny, maxy)
    };
    if (x1 - x0).abs() < 1e-12 || (y1 - y0).abs() < 1e-12 {
        return None;
    }
    let xn: Vec<f64> = x.iter().map(|&v| (v - x0) / (x1 - x0)).collect();
    let yn: Vec<f64> = yv.iter().map(|&v| (v - y0) / (y1 - y0)).collect();

    // Chord from first to last point.
    let (ax, ay) = (xn[0], yn[0]);
    let (bx, by) = (xn[xn.len() - 1], yn[yn.len() - 1]);

    let dx = bx - ax;
    let dy = by - ay;
    let denom = (dx * dx + dy * dy).sqrt();
    if denom <= 0.0 {
        return None;
    }

    let mut best_i = 0usize;
    let mut best_d = -1.0;
    for i in 1..(xn.len() - 1) {
        // Perpendicular distance from point to line AB.
        let px = xn[i];
        let py = yn[i];
        let dist = ((dy * px - dx * py + bx * ay - by * ax).abs()) / denom;
        if dist > best_d {
            best_d = dist;
            best_i = i;
        }
    }

    Some(w_sorted[best_i])
}

fn peak_width_um(widths_um: &[u32], y: &[f64]) -> Option<u32> {
    if widths_um.len() != y.len() || widths_um.is_empty() {
        return None;
    }
    let mut best_y = f64::NEG_INFINITY;
    let mut best_w: Option<u32> = None;
    for (&w, &v) in widths_um.iter().zip(y.iter()) {
        if !v.is_finite() {
            continue;
        }
        if v > best_y {
            best_y = v;
            best_w = Some(w);
        } else if v == best_y {
            // If the curve plateaus, treat the "peak" as the *last* max-width point.
            if best_w.map(|bw| w > bw).unwrap_or(true) {
                best_w = Some(w);
            }
        }
    }
    best_w
}

/// Knee detection with a descending-segment-only search domain:
/// Evaluate candidates only for widths w >= max(knee_search_min_um, w_peak + descending_delta_um).
pub fn knee_max_distance_chord_descending(
    widths_um: &[u32],
    y: &[f64],
    knee_search_min_um: u32,
    descending_delta_um: u32,
) -> Option<u32> {
    if widths_um.len() != y.len() || widths_um.len() < 3 {
        return None;
    }
    let Some(w_peak) = peak_width_um(widths_um, y) else {
        return None;
    };
    let min_w = knee_search_min_um.max(w_peak.saturating_add(descending_delta_um));

    let mut w_f: Vec<u32> = Vec::new();
    let mut y_f: Vec<f64> = Vec::new();
    for (&w, &v) in widths_um.iter().zip(y.iter()) {
        if w >= min_w {
            w_f.push(w);
            y_f.push(v);
        }
    }
    knee_max_distance_chord(&w_f, &y_f)
}
