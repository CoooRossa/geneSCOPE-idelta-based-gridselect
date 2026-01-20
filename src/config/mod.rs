use anyhow::{anyhow, Context, Result};
use std::fs;
use std::path::Path;

pub fn read_widths_um<P: AsRef<Path>>(path: P) -> Result<Vec<u32>> {
    let raw = fs::read_to_string(&path)
        .with_context(|| format!("failed to read widths file: {}", path.as_ref().display()))?;

    let mut widths: Vec<u32> = Vec::new();
    for (line_idx, line) in raw.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let w: f64 = line
            .parse()
            .map_err(|e| anyhow!("invalid width at line {}: {} ({})", line_idx + 1, line, e))?;
        if !w.is_finite() || w <= 0.0 {
            return Err(anyhow!("invalid width at line {}: {}", line_idx + 1, line));
        }
        let w_rounded = w.round();
        if (w - w_rounded).abs() > 1e-9 {
            return Err(anyhow!(
                "width must be an integer number of µm; got {} at line {}",
                line,
                line_idx + 1
            ));
        }
        widths.push(w_rounded as u32);
    }

    widths.sort_unstable();
    widths.dedup();
    if widths.is_empty() {
        return Err(anyhow!("width set is empty: {}", path.as_ref().display()));
    }
    Ok(widths)
}
