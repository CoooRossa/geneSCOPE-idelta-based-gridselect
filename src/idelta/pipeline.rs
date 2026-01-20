use crate::idelta::knee::knee_max_distance_chord_descending;
use crate::idelta::smooth::gaussian_smooth_logx;
use crate::io::molecules::RawMolecule;
use crate::roi::raster::{GridSpec, RoiMask};
use crate::roi::tiles::{list_valid_tiles, TileId};
use anyhow::{anyhow, Context, Result};
use rand::prelude::*;
use rayon::prelude::*;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs;
use std::path::{Path, PathBuf};

#[derive(Debug, Clone)]
pub struct PipelineConfig {
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
    pub widths_prescreen_um: Vec<u32>,
    pub widths_anchor_um: Vec<u32>,
    pub widths_prescreen_version: String,
    pub widths_anchor_version: String,
    pub b_min: u32,
    pub n_min: u32,
    pub min_points: usize,
    pub min_informative_frac_ge2: f64,
}

#[derive(Debug, Clone)]
struct MoleculeCell {
    gene_id: usize,
    cell_i: u32,
    cell_j: u32,
    tile: TileId,
}

#[derive(Debug, Clone, Copy)]
struct BinDiag {
    sum_n_n1: u64,
    n_bins_ge2: u32,
    max_bin_count: u32,
}

#[derive(Debug, Clone)]
struct CurveRow {
    gene: String,
    width_um: u32,
    n_total: u64,
    q_bins: u64,
    n_bins_ge2: u32,
    frac_bins_ge2: f64,
    max_bin_count: u32,
    sum_n_n1: u64,
    idelta_raw: Option<f64>,
    idelta_smooth: Option<f64>,
    status: &'static str,
}

#[derive(Debug, Clone)]
struct GeneSummary {
    gene: String,
    knee_um: Option<u32>,
    gene_status: &'static str,
}

pub fn recommend_grid(cfg: PipelineConfig) -> Result<()> {
    fs::create_dir_all(&cfg.out_dir)
        .with_context(|| format!("failed creating out dir: {}", cfg.out_dir))?;

    if cfg.tile_um % cfg.w0_um != 0 {
        return Err(anyhow!(
            "--tile-um ({}) must be a multiple of --w0-um ({})",
            cfg.tile_um,
            cfg.w0_um
        ));
    }
    for &w in cfg
        .widths_prescreen_um
        .iter()
        .chain(cfg.widths_anchor_um.iter())
    {
        if w % cfg.w0_um != 0 {
            return Err(anyhow!(
                "width {}µm is not a multiple of w0 {}µm",
                w,
                cfg.w0_um
            ));
        }
    }

    // Step 1: ROI parse + w0 mask.
    let poly = crate::io::roi::read_roi_polygon(&cfg.roi_path)?;
    let grid = GridSpec::from_polygon_bbox(&poly, cfg.w0_um)?;
    let roi = RoiMask::build(&poly, grid);

    // Step 1: read + filter molecules.
    let raw = crate::io::molecules::read_molecules(&cfg.molecules_path)?;
    let (gene_names, molecules) = filter_and_index_molecules(&raw, &poly, &roi, cfg.tile_um)?;
    let gene_order = sorted_gene_order(&gene_names);

    if gene_names.is_empty() || molecules.is_empty() {
        return Err(anyhow!("no molecules remained after ROI filtering"));
    }

    // Step 2: prescreen global knee (deterministic).
    let q_prescreen: Vec<u64> = cfg
        .widths_prescreen_um
        .iter()
        .copied()
        .map(|w| compute_q_bins(&roi, w))
        .collect::<Result<Vec<u64>>>()?;
    let prescreen_global_knee_um = prescreen_global_knee(
        &cfg,
        &roi,
        &gene_names,
        &molecules,
        &q_prescreen,
    )?
    .unwrap_or_else(|| cfg.widths_prescreen_um[cfg.widths_prescreen_um.len() / 2]);

    // Step 3: build per-gene width sets:
    // W_final(g) = W_anchor ∪ W_window_dense(g), where W_window_dense(g) is all integer widths
    // within ±knee_window_um of the gene's *anchor-only* knee estimate.
    let min_w = cfg
        .widths_anchor_um
        .iter()
        .copied()
        .chain(cfg.widths_prescreen_um.iter().copied())
        .min()
        .unwrap_or(cfg.w0_um);
    let max_w = cfg
        .widths_anchor_um
        .iter()
        .copied()
        .chain(cfg.widths_prescreen_um.iter().copied())
        .max()
        .unwrap_or(min_w);

    let gene_count = gene_names.len();
    let mut gene_cells: Vec<Vec<(u32, u32)>> = vec![Vec::new(); gene_count];
    for m in &molecules {
        gene_cells[m.gene_id].push((m.cell_i, m.cell_j));
    }

    // Precompute Q(w) for anchor widths.
    let widths_anchor_um = cfg.widths_anchor_um.clone();
    let widths_anchor_k: Vec<u32> = widths_anchor_um.iter().map(|&w| w / cfg.w0_um).collect();
    let q_anchor: Vec<u64> = widths_anchor_um
        .iter()
        .copied()
        .map(|w| compute_q_bins(&roi, w))
        .collect::<Result<Vec<u64>>>()?;

    // Anchor-only per-gene knees (used only to define dense windows).
    let knee_anchor_by_gene: Vec<Option<u32>> = (0..gene_count)
        .into_par_iter()
        .map(|g| {
            let cells = &gene_cells[g];
            let (_valid_points, knee) = gene_knee_on_widths(
                &cfg,
                cells,
                &widths_anchor_um,
                &widths_anchor_k,
                &q_anchor,
            );
            knee
        })
        .collect();

    // Final width sets per gene: anchor ∪ dense window around the anchor-only knee.
    let widths_final_by_gene: Vec<Vec<u32>> = (0..gene_count)
        .into_par_iter()
        .map(|g| {
            let mut set: BTreeSet<u32> = cfg.widths_anchor_um.iter().copied().collect();
            if let Some(knee) = knee_anchor_by_gene[g] {
                let lo = knee.saturating_sub(cfg.knee_window_um).max(min_w);
                let hi = (knee + cfg.knee_window_um).min(max_w);
                for w in lo..=hi {
                    if w % cfg.w0_um == 0 {
                        set.insert(w);
                    }
                }
            }
            set.into_iter().collect()
        })
        .collect();

    // Q(w) for all widths used by any gene.
    let mut widths_union: BTreeSet<u32> = BTreeSet::new();
    for ws in &widths_final_by_gene {
        widths_union.extend(ws.iter().copied());
    }
    let mut q_by_width: HashMap<u32, u64> = HashMap::new();
    for &w in &widths_union {
        q_by_width.insert(w, compute_q_bins(&roi, w)?);
    }

    // Step 4: full-ROI final curves + gene-wise knees (per-gene width sets).
    let rows_by_gene: Vec<Vec<CurveRow>> = (0..gene_count)
        .into_par_iter()
        .map(|g| {
            let gene = gene_names[g].clone();
            let cells = &gene_cells[g];
            let widths_um = &widths_final_by_gene[g];
            let widths_k: Vec<u32> = widths_um.iter().map(|&w| w / cfg.w0_um).collect();
            let q_vec: Vec<u64> = widths_um
                .iter()
                .map(|w| q_by_width.get(w).copied().unwrap_or(0))
                .collect();
            compute_rows_for_gene(&cfg, gene, cells, widths_um, &widths_k, &q_vec)
        })
        .collect();

    let gene_summaries: Vec<GeneSummary> = (0..gene_count)
        .into_par_iter()
        .map(|g| {
            let gene = gene_names[g].clone();
            let rows = &rows_by_gene[g];
            let valid_points = rows
                .iter()
                .filter(|r| r.status == "ok" && r.idelta_smooth.is_some())
                .count();

            if valid_points < cfg.min_points {
                return GeneSummary {
                    gene,
                    knee_um: None,
                    gene_status: "excluded_degenerate",
                };
            }

            let widths: Vec<u32> = rows.iter().map(|r| r.width_um).collect();
            let y: Vec<f64> = rows
                .iter()
                .map(|r| r.idelta_smooth.unwrap_or(f64::NAN))
                .collect();
            let knee = knee_max_distance_chord_descending(
                &widths,
                &y,
                cfg.knee_search_min_um,
                cfg.descending_delta_um,
            );
            let gene_status = if knee.is_some() { "informative" } else { "no_knee" };
            GeneSummary {
                gene,
                knee_um: knee,
                gene_status,
            }
        })
        .collect();

    let informative_gene_ids: Vec<usize> = gene_order
        .iter()
        .copied()
        .filter(|&g| gene_summaries[g].gene_status == "informative")
        .collect();

    // Mean curve (diagnostic; computed across informative genes where available).
    let (avg_widths_um, mean_curve) =
        compute_mean_curve_sparse(&rows_by_gene, &informative_gene_ids);
    let knee_of_mean_curve_um = knee_max_distance_chord_descending(
        &avg_widths_um,
        &mean_curve,
        cfg.knee_search_min_um,
        cfg.descending_delta_um,
    );

    let (recommended_grid_um, rec_status, rec_note, rec_frac_ge2) = recommend_and_summarize(
        &cfg,
        &rows_by_gene,
        &gene_summaries,
        &informative_gene_ids,
        &gene_cells,
    );

    // Write outputs (stable, minimal bundle).
    write_final_curves_gz(
        Path::new(&cfg.out_dir).join("final_curves.tsv.gz"),
        &rows_by_gene,
        &gene_order,
    )?;
    write_avg_curve(
        Path::new(&cfg.out_dir).join("avg_curve.tsv"),
        &avg_widths_um,
        &mean_curve,
    )?;
    write_knee_summary(
        Path::new(&cfg.out_dir).join("knee_summary.tsv"),
        &gene_summaries,
        &gene_order,
        KneeSummaryMeta {
            prescreen_global_knee_um,
            knee_of_mean_curve_um,
            recommended_grid_um,
            recommendation_status: rec_status,
            recommendation_note: rec_note,
            recommendation_frac_ge2: rec_frac_ge2,
            n_genes_total: gene_names.len(),
            n_genes_informative: informative_gene_ids.len(),
            n_genes_excluded_degenerate: gene_names.len() - informative_gene_ids.len(),
            b_min: cfg.b_min,
            n_min: cfg.n_min,
            min_points: cfg.min_points,
            min_informative_frac_ge2: cfg.min_informative_frac_ge2,
            knee_search_min_um: cfg.knee_search_min_um,
            descending_delta_um: cfg.descending_delta_um,
            widths_anchor_version: cfg.widths_anchor_version.clone(),
        },
    )?;

    if let Some(w) = recommended_grid_um {
        println!("{w}");
    } else {
        println!("NA");
    }

    Ok(())
}

fn sorted_gene_order(gene_names: &[String]) -> Vec<usize> {
    let mut idx: Vec<usize> = (0..gene_names.len()).collect();
    idx.sort_unstable_by(|&a, &b| gene_names[a].cmp(&gene_names[b]));
    idx
}

fn filter_and_index_molecules(
    raw: &[RawMolecule],
    poly: &crate::roi::polygon::Polygon,
    roi: &RoiMask,
    tile_um: u32,
) -> Result<(Vec<String>, Vec<MoleculeCell>)> {
    let tile_100 = tile_um as i64 * crate::COORD_SCALE_PER_UM;
    let tile_k = (tile_100 / roi.grid.w0_100) as u32;
    if tile_k == 0 {
        return Err(anyhow!("tile size too small vs w0"));
    }

    let (min_x, max_x, min_y, max_y) = poly.bbox_100();

    let mut gene_to_id: HashMap<String, usize> = HashMap::new();
    let mut gene_names: Vec<String> = Vec::new();
    let mut out: Vec<MoleculeCell> = Vec::new();

    for m in raw {
        // bbox prefilter
        if m.x_100 < min_x || m.x_100 > max_x || m.y_100 < min_y || m.y_100 > max_y {
            continue;
        }
        if !poly.contains_point(m.x_100, m.y_100) {
            continue;
        }
        let Some((cell_i, cell_j)) = roi.grid.point_to_cell(m.x_100, m.y_100) else {
            continue;
        };

        let gene_id = match gene_to_id.get(&m.gene) {
            Some(&id) => id,
            None => {
                let id = gene_names.len();
                gene_names.push(m.gene.clone());
                gene_to_id.insert(m.gene.clone(), id);
                id
            }
        };

        let tile = TileId {
            tx: cell_i / tile_k,
            ty: cell_j / tile_k,
        };
        out.push(MoleculeCell {
            gene_id,
            cell_i,
            cell_j,
            tile,
        });
    }

    Ok((gene_names, out))
}

fn compute_q_bins(roi: &RoiMask, width_um: u32) -> Result<u64> {
    let w_100 = width_um as i64 * crate::COORD_SCALE_PER_UM;
    if w_100 % roi.grid.w0_100 != 0 {
        return Err(anyhow!(
            "width {}µm is not a multiple of w0 {}µm",
            width_um,
            roi.grid.w0_100 / crate::COORD_SCALE_PER_UM
        ));
    }
    let k = (w_100 / roi.grid.w0_100) as usize;
    if k == 0 {
        return Err(anyhow!("invalid width"));
    }

    let nx_bins = (roi.grid.nx + k - 1) / k;
    let ny_bins = (roi.grid.ny + k - 1) / k;

    let mut q = 0u64;
    for by in 0..ny_bins {
        for bx in 0..nx_bins {
            let i0 = bx * k;
            let j0 = by * k;
            let i1 = ((bx + 1) * k).min(roi.grid.nx);
            let j1 = ((by + 1) * k).min(roi.grid.ny);
            let s = roi.sum_in_block(i0, j0, i1, j1);
            if s > 0 {
                q += 1;
            }
        }
    }
    Ok(q)
}

fn bin_diag_for_width(cells: &[(u32, u32)], k: u32) -> BinDiag {
    if k == 0 || cells.is_empty() {
        return BinDiag {
            sum_n_n1: 0,
            n_bins_ge2: 0,
            max_bin_count: 0,
        };
    }

    // Accumulate per-bin counts; update diagnostics incrementally.
    let mut counts: HashMap<u64, u32> = HashMap::new();
    let mut sum_n_n1 = 0u64;
    let mut n_bins_ge2 = 0u32;
    let mut max_bin_count = 0u32;

    for &(i, j) in cells {
        let bi = i / k;
        let bj = j / k;
        let key = ((bi as u64) << 32) | (bj as u64);

        let entry = counts.entry(key).or_insert(0);
        let old = *entry;
        let new = old + 1;
        *entry = new;

        // For this bin, contribution c(c-1) increases by 2*old when count increments.
        sum_n_n1 += 2 * (old as u64);
        if new == 2 {
            n_bins_ge2 += 1;
        }
        max_bin_count = max_bin_count.max(new);
    }

    BinDiag {
        sum_n_n1,
        n_bins_ge2,
        max_bin_count,
    }
}

fn idelta_from_counts(q: u64, sum_n_n1: u64, n: u64) -> Option<f64> {
    if n < 2 || q == 0 {
        return None;
    }
    let denom = (n as f64) * ((n - 1) as f64);
    Some((q as f64) * (sum_n_n1 as f64) / denom)
}

fn compute_rows_for_gene(
    cfg: &PipelineConfig,
    gene: String,
    cells: &[(u32, u32)],
    widths_um: &[u32],
    widths_k: &[u32],
    q_by_index: &[u64],
) -> Vec<CurveRow> {
    assert_eq!(widths_um.len(), widths_k.len());
    assert_eq!(widths_um.len(), q_by_index.len());

    let n = cells.len() as u64;
    let w_len = widths_um.len();

    let mut raw: Vec<f64> = vec![f64::NAN; w_len];
    let mut status: Vec<&'static str> = vec!["insufficient_n"; w_len];
    let mut sum_n_n1: Vec<u64> = vec![0u64; w_len];
    let mut n_bins_ge2: Vec<u32> = vec![0u32; w_len];
    let mut max_bin_count: Vec<u32> = vec![0u32; w_len];

    for (i, &k) in widths_k.iter().enumerate() {
        let q = q_by_index[i];
        let diag = bin_diag_for_width(cells, k);
        sum_n_n1[i] = diag.sum_n_n1;
        n_bins_ge2[i] = diag.n_bins_ge2;
        max_bin_count[i] = diag.max_bin_count;

        let st = if n < 2 {
            "insufficient_n"
        } else if q == 0 {
            "invalid_q"
        } else if diag.n_bins_ge2 < cfg.b_min {
            "degenerate_01"
        } else if diag.sum_n_n1 == 0 && n >= cfg.n_min as u64 {
            "degenerate_01"
        } else {
            "ok"
        };
        status[i] = st;

        if st == "ok" {
            raw[i] = idelta_from_counts(q, diag.sum_n_n1, n).unwrap_or(f64::NAN);
        }
    }

    let smooth = gaussian_smooth_logx(widths_um, &raw, 0.35);

    let mut rows: Vec<CurveRow> = Vec::with_capacity(w_len);
    for i in 0..w_len {
        let q = q_by_index[i];
        let frac = if q > 0 {
            (n_bins_ge2[i] as f64) / (q as f64)
        } else {
            f64::NAN
        };
        rows.push(CurveRow {
            gene: gene.clone(),
            width_um: widths_um[i],
            n_total: n,
            q_bins: q,
            n_bins_ge2: n_bins_ge2[i],
            frac_bins_ge2: frac,
            max_bin_count: max_bin_count[i],
            sum_n_n1: sum_n_n1[i],
            idelta_raw: if raw[i].is_finite() { Some(raw[i]) } else { None },
            idelta_smooth: if smooth[i].is_finite() { Some(smooth[i]) } else { None },
            status: status[i],
        });
    }
    rows
}

fn gene_knee_on_widths(
    cfg: &PipelineConfig,
    cells: &[(u32, u32)],
    widths_um: &[u32],
    widths_k: &[u32],
    q_by_index: &[u64],
) -> (usize, Option<u32>) {
    assert_eq!(widths_um.len(), widths_k.len());
    assert_eq!(widths_um.len(), q_by_index.len());

    let n = cells.len() as u64;
    let w_len = widths_um.len();
    let mut raw: Vec<f64> = vec![f64::NAN; w_len];
    let mut status: Vec<&'static str> = vec!["insufficient_n"; w_len];

    for (i, &k) in widths_k.iter().enumerate() {
        let q = q_by_index[i];
        let diag = bin_diag_for_width(cells, k);

        let st = if n < 2 {
            "insufficient_n"
        } else if q == 0 {
            "invalid_q"
        } else if diag.n_bins_ge2 < cfg.b_min {
            "degenerate_01"
        } else if diag.sum_n_n1 == 0 && n >= cfg.n_min as u64 {
            "degenerate_01"
        } else {
            "ok"
        };
        status[i] = st;
        if st == "ok" {
            raw[i] = idelta_from_counts(q, diag.sum_n_n1, n).unwrap_or(f64::NAN);
        }
    }

    let smooth = gaussian_smooth_logx(widths_um, &raw, 0.35);
    let valid_points = status
        .iter()
        .zip(smooth.iter())
        .filter(|(&st, &v)| st == "ok" && v.is_finite())
        .count();

    if valid_points < cfg.min_points {
        return (valid_points, None);
    }

    let knee = knee_max_distance_chord_descending(
        widths_um,
        &smooth,
        cfg.knee_search_min_um,
        cfg.descending_delta_um,
    );
    (valid_points, knee)
}

fn compute_mean_curve_sparse(
    rows_by_gene: &[Vec<CurveRow>],
    informative_gene_ids: &[usize],
) -> (Vec<u32>, Vec<f64>) {
    let mut sum_cnt: BTreeMap<u32, (f64, u64)> = BTreeMap::new();
    for &g in informative_gene_ids {
        for r in &rows_by_gene[g] {
            if r.status != "ok" {
                continue;
            }
            let Some(v) = r.idelta_smooth else {
                continue;
            };
            let entry = sum_cnt.entry(r.width_um).or_insert((0.0, 0));
            entry.0 += v;
            entry.1 += 1;
        }
    }

    let widths: Vec<u32> = sum_cnt.keys().copied().collect();
    let mean: Vec<f64> = sum_cnt
        .values()
        .map(|(s, c)| if *c > 0 { *s / (*c as f64) } else { f64::NAN })
        .collect();
    (widths, mean)
}

fn prescreen_global_knee(
    cfg: &PipelineConfig,
    roi: &RoiMask,
    gene_names: &[String],
    molecules: &[MoleculeCell],
    q_bins_by_w: &[u64],
) -> Result<Option<u32>> {
    let tile_100 = cfg.tile_um as i64 * crate::COORD_SCALE_PER_UM;
    let tile_k = (tile_100 / roi.grid.w0_100) as usize;
    if tile_k == 0 {
        return Err(anyhow!("tile size too small vs w0"));
    }

    let valid_tiles = list_valid_tiles(roi, tile_k);
    if valid_tiles.is_empty() {
        return Err(anyhow!("ROI had no valid tiles at tile_um={}", cfg.tile_um));
    }
    let tile_to_index: HashMap<u64, usize> = valid_tiles
        .iter()
        .enumerate()
        .map(|(i, t)| (t.to_u64(), i))
        .collect();

    let mut tile_mols: Vec<Vec<(usize, (u32, u32))>> = vec![Vec::new(); valid_tiles.len()];
    for m in molecules {
        let tid = m.tile.to_u64();
        let Some(&ti) = tile_to_index.get(&tid) else {
            continue;
        };
        tile_mols[ti].push((m.gene_id, (m.cell_i, m.cell_j)));
    }

    let gene_count = gene_names.len();
    let wv = &cfg.widths_prescreen_um;
    let w_len = wv.len();
    if q_bins_by_w.len() != w_len {
        return Err(anyhow!("internal error: Q(w) length mismatch"));
    }

    let sample_n = ((cfg.prescreen_fraction * valid_tiles.len() as f64).ceil() as usize)
        .clamp(1, valid_tiles.len());

    let mut sum_values: Vec<Vec<u64>> = vec![Vec::new(); gene_count * w_len];
    let mut n_values: Vec<Vec<u64>> = vec![Vec::new(); gene_count];

    let mut gene_cells: Vec<Vec<(u32, u32)>> = vec![Vec::new(); gene_count];
    let mut gene_n: Vec<u64> = vec![0; gene_count];

    for rep in 0..cfg.prescreen_reps {
        for v in gene_cells.iter_mut() {
            v.clear();
        }
        gene_n.fill(0);

        let mut idx: Vec<usize> = (0..valid_tiles.len()).collect();
        let mut rng = StdRng::seed_from_u64(
            cfg.seed ^ (rep as u64).wrapping_mul(0x9E3779B97F4A7C15),
        );
        idx.shuffle(&mut rng);
        idx.truncate(sample_n);

        for &ti in &idx {
            for &(gene_id, cell) in &tile_mols[ti] {
                gene_cells[gene_id].push(cell);
                gene_n[gene_id] += 1;
            }
        }

        let active: Vec<usize> = (0..gene_count).filter(|&g| gene_n[g] >= 2).collect();
        let rep_sums: Vec<Vec<u64>> = active
            .par_iter()
            .map(|&g| {
                let cells = &gene_cells[g];
                wv.iter()
                    .map(|&w| {
                        let k = (w / cfg.w0_um) as u32;
                        bin_diag_for_width(cells, k).sum_n_n1
                    })
                    .collect::<Vec<u64>>()
            })
            .collect();

        for (ai, &g) in active.iter().enumerate() {
            n_values[g].push(gene_n[g]);
            for wi in 0..w_len {
                sum_values[g * w_len + wi].push(rep_sums[ai][wi]);
            }
        }
    }

    let mut knees: Vec<u32> = Vec::new();
    for g in 0..gene_count {
        let n_med = median_u64(&mut n_values[g]);
        if n_med < 2 {
            continue;
        }

        let mut raw: Vec<f64> = vec![f64::NAN; w_len];
        for wi in 0..w_len {
            let sum_med = median_u64(&mut sum_values[g * w_len + wi]);
            raw[wi] = idelta_from_counts(q_bins_by_w[wi], sum_med, n_med).unwrap_or(f64::NAN);
        }
        let smooth = gaussian_smooth_logx(wv, &raw, 0.35);
        if let Some(knee) = knee_max_distance_chord_descending(
            wv,
            &smooth,
            cfg.knee_search_min_um,
            cfg.descending_delta_um,
        ) {
            knees.push(knee);
        }
    }

    knees.sort_unstable();
    Ok(knees.get(knees.len() / 2).copied())
}

fn recommend_and_summarize(
    cfg: &PipelineConfig,
    rows_by_gene: &[Vec<CurveRow>],
    gene_summaries: &[GeneSummary],
    informative_gene_ids: &[usize],
    gene_cells: &[Vec<(u32, u32)>],
) -> (Option<u32>, String, String, Option<f64>) {
    if informative_gene_ids.is_empty() {
        let max_w = cfg.widths_anchor_um.iter().copied().max();
        return (
            max_w,
            "no_informative_genes_fallback_max_width".to_string(),
            "no informative genes after degenerate-point filtering; using max anchor width"
                .to_string(),
            None,
        );
    }

    // Gene-wise knees across informative genes.
    let mut knees: Vec<u32> = informative_gene_ids
        .iter()
        .filter_map(|&g| gene_summaries.get(g).and_then(|s| s.knee_um))
        .collect();
    knees.sort_unstable();

    let mut recommended_grid_um = knees.get(knees.len() / 2).copied();
    let mut rec_status = "median_gene_knees".to_string();
    let mut rec_note =
        "recommended_grid_um is the median of informative gene knees (overall_knee_um is the mean)"
            .to_string();

    // Sanity check: require >= min_informative_frac_ge2 of informative genes to have n_bins_ge2 >= b_min at the recommendation.
    let frac_ge2_at = |w: u32| -> Option<f64> {
        let mut pass = 0u64;
        for &g in informative_gene_ids {
            let bins_ge2 = if let Some(r) = find_row_by_width(&rows_by_gene[g], w) {
                r.n_bins_ge2
            } else {
                // Not evaluated for this gene: compute on-the-fly for the sanity check only.
                let k = w / cfg.w0_um;
                bin_diag_for_width(&gene_cells[g], k).n_bins_ge2
            };
            if bins_ge2 >= cfg.b_min {
                pass += 1;
            }
        }
        Some(pass as f64 / informative_gene_ids.len() as f64)
    };

    let mut rec_frac_ge2 = recommended_grid_um.and_then(frac_ge2_at);
    if let (Some(rec), Some(frac)) = (recommended_grid_um, rec_frac_ge2) {
        if frac < cfg.min_informative_frac_ge2 {
            // Deterministic conservative fallback: scan anchor widths to find the smallest larger width that passes.
            let mut candidates = cfg.widths_anchor_um.clone();
            candidates.sort_unstable();
            candidates.dedup();

            let mut found = None;
            for &w2 in &candidates {
                if w2 <= rec {
                    continue;
                }
                if let Some(f2) = frac_ge2_at(w2) {
                    if f2 >= cfg.min_informative_frac_ge2 {
                        found = Some((w2, f2));
                        break;
                    }
                }
            }
            match found {
                Some((w2, f2)) => {
                    rec_status = "fallback_increase_width".to_string();
                    rec_note = format!(
                        "sanity_check failed at {}µm (frac_ge2={:.3}); moved to {}µm (frac_ge2={:.3}); policy=median_gene_knees",
                        rec, frac, w2, f2
                    );
                    recommended_grid_um = Some(w2);
                    rec_frac_ge2 = Some(f2);
                }
                None => {
                    let max_w = candidates.last().copied().unwrap_or(rec);
                    let f2 = frac_ge2_at(max_w);
                    rec_status = "fallback_max_width".to_string();
                    rec_note = format!(
                        "sanity_check failed at {}µm (frac_ge2={:.3}); no larger anchor width passed; using max anchor width {}µm; policy=median_gene_knees",
                        rec, frac, max_w
                    );
                    recommended_grid_um = Some(max_w);
                    rec_frac_ge2 = f2;
                }
            }
        }
    }

    (recommended_grid_um, rec_status, rec_note, rec_frac_ge2)
}

fn find_row_by_width<'a>(rows: &'a [CurveRow], w: u32) -> Option<&'a CurveRow> {
    let mut lo = 0usize;
    let mut hi = rows.len();
    while lo < hi {
        let mid = (lo + hi) / 2;
        let mw = rows[mid].width_um;
        if mw < w {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    if lo < rows.len() && rows[lo].width_um == w {
        Some(&rows[lo])
    } else {
        None
    }
}

fn median_u64(v: &mut Vec<u64>) -> u64 {
    if v.is_empty() {
        return 0;
    }
    v.sort_unstable();
    v[v.len() / 2]
}

fn opt_f64(v: Option<f64>) -> String {
    match v {
        Some(x) => format!("{:.10}", x),
        None => "NA".to_string(),
    }
}

fn final_curves_header() -> &'static str {
    "gene\twidth_um\tN_total\tQ_bins\tn_bins_ge2\tfrac_bins_ge2\tmax_bin_count\tsum_n_n1\tidelta_raw\tidelta_smooth\tstatus"
}

fn write_final_curves_gz(
    path: PathBuf,
    rows_by_gene: &[Vec<CurveRow>],
    gene_order: &[usize],
) -> Result<()> {
    let rows = gene_order.iter().flat_map(|&g| {
        rows_by_gene[g].iter().map(|r| {
            format!(
                "{}\t{}\t{}\t{}\t{}\t{:.10}\t{}\t{}\t{}\t{}\t{}",
                r.gene,
                r.width_um,
                r.n_total,
                r.q_bins,
                r.n_bins_ge2,
                r.frac_bins_ge2,
                r.max_bin_count,
                r.sum_n_n1,
                opt_f64(r.idelta_raw),
                opt_f64(r.idelta_smooth),
                r.status
            )
        })
    });
    crate::io::tsv::write_tsv_gz(path, final_curves_header(), rows)
}

fn write_avg_curve(path: PathBuf, widths_um: &[u32], mean_curve: &[f64]) -> Result<()> {
    let header = "width_um\tidelta_mean";
    let rows = widths_um.iter().zip(mean_curve.iter()).map(|(&w, &v)| {
        let s = if v.is_finite() {
            format!("{:.10}", v)
        } else {
            "NA".to_string()
        };
        format!("{w}\t{s}")
    });
    crate::io::tsv::write_tsv(path, header, rows)
}

struct KneeSummaryMeta {
    prescreen_global_knee_um: u32,
    knee_of_mean_curve_um: Option<u32>,
    recommended_grid_um: Option<u32>,
    recommendation_status: String,
    recommendation_note: String,
    recommendation_frac_ge2: Option<f64>,
    n_genes_total: usize,
    n_genes_informative: usize,
    n_genes_excluded_degenerate: usize,
    b_min: u32,
    n_min: u32,
    min_points: usize,
    min_informative_frac_ge2: f64,
    knee_search_min_um: u32,
    descending_delta_um: u32,
    widths_anchor_version: String,
}

fn write_knee_summary(
    path: PathBuf,
    gene_summaries: &[GeneSummary],
    gene_order: &[usize],
    meta: KneeSummaryMeta,
) -> Result<()> {
    let header_cols = [
        // Per-gene fields (populated for gene rows).
        "gene",
        "knee_um",
        "gene_status",
        // Global summary fields (populated only on gene='__SUMMARY__').
        "prescreen_global_knee_um",
        "overall_knee_um",
        "recommended_grid_um",
        "recommendation_frac_ge2",
        "recommendation_status",
        "recommendation_note",
        "knee_mean",
        "knee_sd",
        "knee_min",
        "knee_max",
        "n_genes_total",
        "n_genes_informative",
        "n_genes_excluded_degenerate",
        "b_min",
        "n_min",
        "min_points",
        "min_informative_frac_ge2",
        // Added global fields (summary row only).
        "overall_knee_median_um",
        "overall_knee_trimmed_mean_um",
        "overall_knee_sd_um",
        "overall_knee_iqr_um",
        "knee_of_mean_curve_um",
        "knee_search_min_um",
        "descending_delta_um",
        "widths_anchor_version",
    ];
    let header = header_cols.join("\t");

    let mut knees: Vec<u32> = gene_summaries
        .iter()
        .filter(|g| g.gene_status == "informative")
        .filter_map(|g| g.knee_um)
        .collect();
    knees.sort_unstable();
    let (knee_mean, knee_sd, knee_min, knee_max) = knee_stats(&knees);
    let overall_knee_um = knee_mean;
    let overall_knee_median_um = knees.get(knees.len() / 2).copied();
    let overall_knee_trimmed_mean_um = trimmed_mean_u32(&knees, 0.10);
    let overall_knee_sd_um = knee_sd;
    let overall_knee_iqr_um = iqr_u32(&knees);

    let summary_row = {
        let mut fields: Vec<String> = vec![String::from("NA"); header_cols.len()];
        fields[0] = "__SUMMARY__".to_string();
        fields[1] = "NA".to_string();
        fields[2] = "summary".to_string();
        fields[3] = meta.prescreen_global_knee_um.to_string();
        fields[4] = overall_knee_um.map(|v| format!("{:.10}", v)).unwrap_or_else(|| "NA".to_string());
        fields[5] = meta
            .recommended_grid_um
            .map(|v| v.to_string())
            .unwrap_or_else(|| "NA".to_string());
        fields[6] = meta
            .recommendation_frac_ge2
            .map(|v| format!("{:.6}", v))
            .unwrap_or_else(|| "NA".to_string());
        fields[7] = meta.recommendation_status.clone();
        fields[8] = meta.recommendation_note.replace('\t', " ");
        fields[9] = opt_f64(knee_mean);
        fields[10] = opt_f64(knee_sd);
        fields[11] = knee_min.map(|v| v.to_string()).unwrap_or_else(|| "NA".to_string());
        fields[12] = knee_max.map(|v| v.to_string()).unwrap_or_else(|| "NA".to_string());
        fields[13] = meta.n_genes_total.to_string();
        fields[14] = meta.n_genes_informative.to_string();
        fields[15] = meta.n_genes_excluded_degenerate.to_string();
        fields[16] = meta.b_min.to_string();
        fields[17] = meta.n_min.to_string();
        fields[18] = meta.min_points.to_string();
        fields[19] = format!("{:.6}", meta.min_informative_frac_ge2);
        fields[20] = overall_knee_median_um.map(|v| v.to_string()).unwrap_or_else(|| "NA".to_string());
        fields[21] = overall_knee_trimmed_mean_um.map(|v| format!("{:.10}", v)).unwrap_or_else(|| "NA".to_string());
        fields[22] = overall_knee_sd_um.map(|v| format!("{:.10}", v)).unwrap_or_else(|| "NA".to_string());
        fields[23] = overall_knee_iqr_um.map(|v| format!("{:.10}", v)).unwrap_or_else(|| "NA".to_string());
        fields[24] = meta
            .knee_of_mean_curve_um
            .map(|v| v.to_string())
            .unwrap_or_else(|| "NA".to_string());
        fields[25] = meta.knee_search_min_um.to_string();
        fields[26] = meta.descending_delta_um.to_string();
        fields[27] = meta.widths_anchor_version.clone();
        fields.join("\t")
    };

    let mut rows: Vec<String> = Vec::new();
    rows.push(summary_row);
    rows.extend(gene_order.iter().map(|&g| {
        let gs = &gene_summaries[g];
        let mut fields: Vec<String> = vec![String::from("NA"); header_cols.len()];
        fields[0] = gs.gene.clone();
        fields[1] = gs.knee_um.map(|v| v.to_string()).unwrap_or_else(|| "NA".to_string());
        fields[2] = gs.gene_status.to_string();
        fields.join("\t")
    }));

    crate::io::tsv::write_tsv(path, &header, rows)
}

fn knee_stats(knees: &[u32]) -> (Option<f64>, Option<f64>, Option<u32>, Option<u32>) {
    if knees.is_empty() {
        return (None, None, None, None);
    }
    let n = knees.len() as f64;
    let mean = knees.iter().map(|&v| v as f64).sum::<f64>() / n;
    let var = if knees.len() >= 2 {
        knees
            .iter()
            .map(|&v| {
                let d = v as f64 - mean;
                d * d
            })
            .sum::<f64>()
            / (n - 1.0)
    } else {
        0.0
    };
    let sd = var.sqrt();
    (Some(mean), Some(sd), knees.first().copied(), knees.last().copied())
}

fn trimmed_mean_u32(knees_sorted: &[u32], trim: f64) -> Option<f64> {
    if knees_sorted.is_empty() {
        return None;
    }
    let n = knees_sorted.len();
    let mut k = (trim.max(0.0).min(0.49) * n as f64).floor() as usize;
    if 2 * k >= n {
        k = n / 2;
    }
    let slice = &knees_sorted[k..(n - k)];
    if slice.is_empty() {
        return None;
    }
    Some(slice.iter().map(|&v| v as f64).sum::<f64>() / slice.len() as f64)
}

fn iqr_u32(knees_sorted: &[u32]) -> Option<f64> {
    if knees_sorted.is_empty() {
        return None;
    }
    let n = knees_sorted.len();
    let q1 = knees_sorted[n / 4] as f64;
    let q3 = knees_sorted[(3 * n) / 4] as f64;
    Some(q3 - q1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;

    #[test]
    fn exactness_naive_vs_w0_coarsen_sum_n_n1() {
        // Toy: 2 genes, small coordinates, origin at 0.
        let w0_um = 1u32;
        let x0 = 0.0f64;
        let y0 = 0.0f64;
        let widths = vec![1u32, 2u32, 3u32];

        let mols = vec![
            ("A", 0.10, 0.10),
            ("A", 0.20, 0.20),
            ("A", 1.10, 1.10),
            ("A", 1.20, 1.20),
            ("B", 0.10, 0.10),
            ("B", 2.10, 0.10),
        ];

        let mut cells_by_gene: HashMap<&str, Vec<(u32, u32)>> = HashMap::new();
        let mut coords_by_gene: HashMap<&str, Vec<(f64, f64)>> = HashMap::new();
        for (g, x, y) in mols {
            coords_by_gene.entry(g).or_default().push((x, y));
            let ci = ((x - x0) / w0_um as f64).floor() as u32;
            let cj = ((y - y0) / w0_um as f64).floor() as u32;
            cells_by_gene.entry(g).or_default().push((ci, cj));
        }

        for g in ["A", "B"] {
            let coords = &coords_by_gene[g];
            let cells = &cells_by_gene[g];
            for &w in &widths {
                // Naive: rebin coords directly at width w.
                let mut naive: HashMap<(u32, u32), u32> = HashMap::new();
                for &(x, y) in coords {
                    let bi = ((x - x0) / w as f64).floor() as u32;
                    let bj = ((y - y0) / w as f64).floor() as u32;
                    *naive.entry((bi, bj)).or_insert(0) += 1;
                }
                let naive_sum: u64 =
                    naive.values().map(|&c| (c as u64) * (c as u64 - 1)).sum();

                // w0-based: coarsen from w0 cells via integer division by k.
                let k = (w / w0_um) as u32;
                let coarsen_sum = bin_diag_for_width(cells, k).sum_n_n1;
                assert_eq!(naive_sum, coarsen_sum, "gene={g} w={w}");
            }
        }
    }

    #[test]
    fn q_consistency_rectangle() {
        // Rectangle ROI: (0,0)-(3,2) in µm, w0=1.
        use crate::roi::polygon::{Point, Polygon};
        let poly = Polygon::new_closed(vec![
            Point::from_um_f64(0.0, 0.0),
            Point::from_um_f64(3.0, 0.0),
            Point::from_um_f64(3.0, 2.0),
            Point::from_um_f64(0.0, 2.0),
        ])
        .unwrap();
        let grid = GridSpec::from_polygon_bbox(&poly, 1).unwrap();
        let roi = RoiMask::build(&poly, grid);

        // Direct enumeration for rectangle:
        // Q(w) = number of w-bins whose rectangle overlaps ROI with positive area.
        fn q_rect(w: u32) -> u64 {
            let w = w as i64;
            let roi_x0 = 0i64;
            let roi_x1 = 3i64;
            let roi_y0 = 0i64;
            let roi_y1 = 2i64;
            let nx = ((roi_x1 - roi_x0) + w - 1) / w;
            let ny = ((roi_y1 - roi_y0) + w - 1) / w;
            let mut q = 0u64;
            for by in 0..ny {
                for bx in 0..nx {
                    let x0 = bx * w;
                    let x1 = (bx + 1) * w;
                    let y0 = by * w;
                    let y1 = (by + 1) * w;
                    let ox = (x0 < roi_x1) && (x1 > roi_x0);
                    let oy = (y0 < roi_y1) && (y1 > roi_y0);
                    if ox && oy {
                        q += 1;
                    }
                }
            }
            q
        }

        for &w in &[1u32, 2u32] {
            let q_pool = compute_q_bins(&roi, w).unwrap();
            assert_eq!(q_rect(w), q_pool, "w={w}");
        }
    }
}
