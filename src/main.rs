use anyhow::Result;
use clap::{Parser, Subcommand};

#[derive(Parser, Debug)]
#[command(
    name = "idelta-gridselect",
    about = "Dedicated Morisita Iδ grid-size selector (recommend best grid width)",
    version
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Recommend an optimal grid size (µm) with a minimal diagnostics bundle.
    RecommendGrid(RecommendGridArgs),
}

#[derive(Parser, Debug)]
struct RecommendGridArgs {
    /// Molecules table (CSV/TSV/Parquet). Required columns: gene, x_um, y_um.
    #[arg(long)]
    molecules: String,

    /// ROI polygon path (GeoJSON / WKT / vertices CSV/TSV).
    #[arg(long)]
    roi: String,

    /// Output directory (created if missing).
    #[arg(long)]
    out: String,

    /// Base grid size (µm). Must evenly divide all evaluated widths.
    #[arg(long, default_value_t = 1)]
    w0_um: u32,

    /// Tile size (µm) for prescreen stratification (must be a multiple of w0).
    #[arg(long, default_value_t = 250)]
    tile_um: u32,

    /// Fraction of valid tiles sampled per prescreen repetition.
    #[arg(long, default_value_t = 0.25)]
    prescreen_fraction: f64,

    /// Number of prescreen repetitions.
    #[arg(long, default_value_t = 25)]
    prescreen_reps: u32,

    /// RNG seed for deterministic reproducibility.
    #[arg(long, default_value_t = 1)]
    seed: u64,

    /// Knee window half-width (µm) around the global prescreen knee.
    #[arg(long, default_value_t = 20)]
    knee_window_um: u32,

    /// Minimum knee-search width (µm). Knee search ignores widths < this value.
    #[arg(long, default_value_t = 5)]
    knee_search_min_um: u32,

    /// Descending-segment constraint (µm): knee search starts at w >= w_peak + descending_delta_um.
    #[arg(long, default_value_t = 2)]
    descending_delta_um: u32,

    /// Number of threads (rayon worker threads).
    #[arg(long, default_value_t = 1)]
    threads: usize,

    /// Optional prescreen widths file (one width per line, µm); defaults to embedded `widths_prescreen_v2`.
    #[arg(long)]
    widths_prescreen: Option<String>,

    /// Optional anchor widths file (one width per line, µm); defaults to embedded `widths_anchor_v2`.
    #[arg(long)]
    widths_anchor: Option<String>,

    /// Minimum number of bins with n_i >= 2 required for a non-degenerate point.
    #[arg(long, default_value_t = 50)]
    b_min: u32,

    /// If sum_n_n1 == 0 and N_total >= N_min, mark point as degenerate_01.
    #[arg(long, default_value_t = 100)]
    n_min: u32,

    /// Minimum number of non-degenerate widths required to estimate a gene knee.
    #[arg(long, default_value_t = 6)]
    min_points: usize,

    /// Sanity check: at the recommended width, require this fraction of informative genes to have n_bins_ge2 >= B_min.
    #[arg(long, default_value_t = 0.6)]
    min_informative_frac_ge2: f64,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Commands::RecommendGrid(args) => {
            let cfg = idelta_gridselect::RecommendGridConfig {
                molecules_path: args.molecules,
                roi_path: args.roi,
                out_dir: args.out,
                w0_um: args.w0_um,
                tile_um: args.tile_um,
                prescreen_fraction: args.prescreen_fraction,
                prescreen_reps: args.prescreen_reps,
                seed: args.seed,
                knee_window_um: args.knee_window_um,
                knee_search_min_um: args.knee_search_min_um,
                descending_delta_um: args.descending_delta_um,
                threads: args.threads,
                widths_prescreen_path: args.widths_prescreen,
                widths_anchor_path: args.widths_anchor,
                b_min: args.b_min,
                n_min: args.n_min,
                min_points: args.min_points,
                min_informative_frac_ge2: args.min_informative_frac_ge2,
            };
            idelta_gridselect::recommend_grid(cfg)
        }
    }
}
