use anyhow::{anyhow, Context, Result};
use std::path::Path;

pub mod config;
pub mod idelta;
pub mod io;
pub mod roi;

pub const COORD_SCALE_PER_UM: i64 = 100; // 0.01 µm units

#[derive(Debug, Clone)]
pub struct RecommendGridConfig {
    pub molecules_path: String,
    pub roi_path: String,
    pub out_dir: String,
    pub w0_um: u32,
    pub tile_um: u32,
    pub prescreen_fraction: f64,
    pub prescreen_reps: u32,
    pub seed: u64,
    pub knee_window_um: u32,
    pub knee_search_min_um: u32,
    pub descending_delta_um: u32,
    pub threads: usize,
    pub widths_prescreen_path: Option<String>,
    pub widths_anchor_path: Option<String>,
    pub b_min: u32,
    pub n_min: u32,
    pub min_points: usize,
    pub min_informative_frac_ge2: f64,
}

pub fn recommend_grid(cfg: RecommendGridConfig) -> Result<()> {
    if cfg.w0_um == 0 {
        return Err(anyhow!("--w0-um must be >= 1"));
    }
    if cfg.threads == 0 {
        return Err(anyhow!("--threads must be >= 1"));
    }
    if !(0.0 < cfg.prescreen_fraction && cfg.prescreen_fraction <= 1.0) {
        return Err(anyhow!("--prescreen-fraction must be in (0, 1]"));
    }
    if cfg.prescreen_reps == 0 {
        return Err(anyhow!("--prescreen-reps must be >= 1"));
    }
    if cfg.b_min == 0 {
        return Err(anyhow!("--b-min must be >= 1"));
    }
    if cfg.n_min == 0 {
        return Err(anyhow!("--n-min must be >= 1"));
    }
    if cfg.min_points < 3 {
        return Err(anyhow!("--min-points must be >= 3"));
    }
    if cfg.knee_search_min_um == 0 {
        return Err(anyhow!("--knee-search-min-um must be >= 1"));
    }
    if !(0.0 < cfg.min_informative_frac_ge2 && cfg.min_informative_frac_ge2 <= 1.0) {
        return Err(anyhow!("--min-informative-frac-ge2 must be in (0, 1]"));
    }

    let (widths_prescreen, widths_prescreen_version) = match cfg.widths_prescreen_path.as_deref()
    {
        Some(path) => (
            config::read_widths_um(path)
                .with_context(|| format!("failed reading --widths-prescreen: {}", path))?,
            format!(
                "file:{}",
                Path::new(path)
                    .file_name()
                    .and_then(|s| s.to_str())
                    .unwrap_or("widths_prescreen")
            ),
        ),
        None => (
            config::parse_widths_um_str(
                config::DEFAULT_WIDTHS_PRESCREEN_V2,
                "embedded: widths_prescreen_v2",
            )
            .context("failed parsing embedded prescreen widths")?,
            "embedded:v2".to_string(),
        ),
    };
    let (widths_anchor, widths_anchor_version) = match cfg.widths_anchor_path.as_deref() {
        Some(path) => (
            config::read_widths_um(path)
                .with_context(|| format!("failed reading --widths-anchor: {}", path))?,
            format!(
                "file:{}",
                Path::new(path)
                    .file_name()
                    .and_then(|s| s.to_str())
                    .unwrap_or("widths_anchor")
            ),
        ),
        None => (
            config::parse_widths_um_str(
                config::DEFAULT_WIDTHS_ANCHOR_V2,
                "embedded: widths_anchor_v2",
            )
            .context("failed parsing embedded anchor widths")?,
            "embedded:v2".to_string(),
        ),
    };

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(cfg.threads)
        .build()
        .context("failed creating rayon thread pool")?;

    pool.install(|| {
        crate::idelta::pipeline::recommend_grid(crate::idelta::pipeline::PipelineConfig {
            molecules_path: cfg.molecules_path,
            roi_path: cfg.roi_path,
            out_dir: cfg.out_dir,
            w0_um: cfg.w0_um,
            tile_um: cfg.tile_um,
            prescreen_fraction: cfg.prescreen_fraction,
            prescreen_reps: cfg.prescreen_reps,
            seed: cfg.seed,
            knee_window_um: cfg.knee_window_um,
            knee_search_min_um: cfg.knee_search_min_um,
            descending_delta_um: cfg.descending_delta_um,
            threads: cfg.threads,
            widths_prescreen_um: widths_prescreen,
            widths_anchor_um: widths_anchor,
            widths_prescreen_version,
            widths_anchor_version,
            b_min: cfg.b_min,
            n_min: cfg.n_min,
            min_points: cfg.min_points,
            min_informative_frac_ge2: cfg.min_informative_frac_ge2,
        })
    })
}
