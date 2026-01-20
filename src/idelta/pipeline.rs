use crate::idelta::knee::knee_max_distance_chord;
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
    pub threads: usize,
    pub widths_prescreen_um: Vec<u32>,
    pub widths_anchor_um: Vec<u32>,
    pub widths_prescreen_path: String,
    pub widths_anchor_path: String,
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

    // Step 3: global final width set for all genes (stable bundle).
    let widths_final_um = build_global_final_widths(&cfg, prescreen_global_knee_um)?;
    let widths_final_k: Vec<u32> = widths_final_um.iter().map(|&w| w / cfg.w0_um).collect();
    let width_to_idx: HashMap<u32, usize> =
        widths_final_um.iter().enumerate().map(|(i, &w)| (w, i)).collect();

    // Q(w) for evaluated widths.
    let mut q_by_width: BTreeMap<u32, u64> = BTreeMap::new();
    for &w in &widths_final_um {
        q_by_width.insert(w, compute_q_bins(&roi, w)?);
    }
    let q_vec: Vec<u64> = widths_final_um
        .iter()
        .map(|w| q_by_width.get(w).copied().unwrap_or(0))
        .collect();

    // Step 4: full-ROI final curves on the global width set.
    let (rows_by_gene, gene_summaries, informative_gene_ids) = compute_final_curves_and_knees(
        &cfg,
        &gene_names,
        &molecules,
        &widths_final_um,
        &widths_final_k,
        &q_vec,
    )?;

    // avg_curve.tsv + overall knee of mean curve.
    let mean_curve: Vec<f64> =
        compute_mean_curve(&rows_by_gene, &informative_gene_ids, widths_final_um.len());
    let overall_knee_um = knee_max_distance_chord(&widths_final_um, &mean_curve);

    // Recommend grid size (with sanity check + deterministic fallback).
    let (recommended_grid_um, rec_status, rec_note, rec_frac_ge2) =
        recommend_with_sanity_check(
            &cfg,
            &rows_by_gene,
            &gene_summaries,
            &informative_gene_ids,
            &width_to_idx,
            &widths_final_um,
            overall_knee_um,
        );

    // Write outputs (stable, minimal bundle).
    write_final_curves_gz(
        Path::new(&cfg.out_dir).join("final_curves.tsv.gz"),
        &rows_by_gene,
        &gene_order,
    )?;
    write_avg_curve(
        Path::new(&cfg.out_dir).join("avg_curve.tsv"),
        &widths_final_um,
        &mean_curve,
    )?;
    write_knee_summary(
        Path::new(&cfg.out_dir).join("knee_summary.tsv"),
        &gene_summaries,
        &gene_order,
        KneeSummaryMeta {
            prescreen_global_knee_um,
            overall_knee_um,
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
        if let Some(knee) = knee_max_distance_chord(wv, &smooth) {
            knees.push(knee);
        }
    }

    knees.sort_unstable();
    Ok(knees.get(knees.len() / 2).copied())
}

fn build_global_final_widths(cfg: &PipelineConfig, global_knee_um: u32) -> Result<Vec<u32>> {
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
        .unwrap_or(global_knee_um);

    let lo = global_knee_um.saturating_sub(cfg.knee_window_um).max(min_w);
    let hi = (global_knee_um + cfg.knee_window_um).min(max_w);

    let mut set: BTreeSet<u32> = BTreeSet::new();
    set.extend(cfg.widths_anchor_um.iter().copied());

    for w in lo..=hi {
        if w % cfg.w0_um == 0 {
            set.insert(w);
        }
    }

    let mut widths: Vec<u32> = set.into_iter().collect();
    widths.sort_unstable();
    if widths.is_empty() {
        return Err(anyhow!("final width set is empty"));
    }
    Ok(widths)
}

fn compute_final_curves_and_knees(
    cfg: &PipelineConfig,
    gene_names: &[String],
    molecules: &[MoleculeCell],
    widths_um: &[u32],
    widths_k: &[u32],
    q_by_index: &[u64],
) -> Result<(Vec<Vec<CurveRow>>, Vec<GeneSummary>, Vec<usize>)> {
    let gene_count = gene_names.len();
    if widths_um.len() != widths_k.len() || widths_um.len() != q_by_index.len() {
        return Err(anyhow!("internal error: width/Q length mismatch"));
    }

    // Collect full-ROI cells per gene.
    let mut gene_cells: Vec<Vec<(u32, u32)>> = vec![Vec::new(); gene_count];
    for m in molecules {
        gene_cells[m.gene_id].push((m.cell_i, m.cell_j));
    }

    let rows_by_gene: Vec<Vec<CurveRow>> = (0..gene_count)
        .into_par_iter()
        .map(|g| {
            let gene = gene_names[g].clone();
            let cells = &gene_cells[g];
            let n = cells.len() as u64;

            let mut raw: Vec<f64> = vec![f64::NAN; widths_um.len()];
            let mut status: Vec<&'static str> = vec!["insufficient_n"; widths_um.len()];
            let mut sum_n_n1: Vec<u64> = vec![0u64; widths_um.len()];
            let mut n_bins_ge2: Vec<u32> = vec![0u32; widths_um.len()];
            let mut max_bin_count: Vec<u32> = vec![0u32; widths_um.len()];

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

            let mut rows: Vec<CurveRow> = Vec::with_capacity(widths_um.len());
            for i in 0..widths_um.len() {
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
            let knee = knee_max_distance_chord(&widths, &y);
            let gene_status = if knee.is_some() { "informative" } else { "no_knee" };
            GeneSummary {
                gene,
                knee_um: knee,
                gene_status,
            }
        })
        .collect();

    let gene_order = sorted_gene_order(gene_names);
    let informative_gene_ids: Vec<usize> = gene_order
        .into_iter()
        .filter(|&g| gene_summaries[g].gene_status == "informative")
        .collect();

    Ok((rows_by_gene, gene_summaries, informative_gene_ids))
}

fn compute_mean_curve(
    rows_by_gene: &[Vec<CurveRow>],
    informative_gene_ids: &[usize],
    width_len: usize,
) -> Vec<f64> {
    let mut sum = vec![0.0f64; width_len];
    let mut cnt = vec![0u64; width_len];

    // Deterministic sum order: informative_gene_ids is already sorted by gene name.
    for &g in informative_gene_ids {
        for (i, r) in rows_by_gene[g].iter().enumerate() {
            if let Some(v) = r.idelta_smooth {
                sum[i] += v;
                cnt[i] += 1;
            }
        }
    }

    (0..width_len)
        .map(|i| if cnt[i] > 0 { sum[i] / (cnt[i] as f64) } else { f64::NAN })
        .collect()
}

fn recommend_with_sanity_check(
    cfg: &PipelineConfig,
    rows_by_gene: &[Vec<CurveRow>],
    gene_summaries: &[GeneSummary],
    informative_gene_ids: &[usize],
    width_to_idx: &HashMap<u32, usize>,
    widths_um: &[u32],
    overall_knee_um: Option<u32>,
) -> (Option<u32>, String, String, Option<f64>) {
    if informative_gene_ids.is_empty() {
        let max_w = widths_um.last().copied();
        return (
            max_w,
            "no_informative_genes_fallback_max_width".to_string(),
            "no informative genes after degenerate-point filtering; using max evaluated width"
                .to_string(),
            None,
        );
    }

    let mut recommended = overall_knee_um;
    if recommended.is_none() {
        // Fallback: median of informative gene knees.
        let mut knees: Vec<u32> = informative_gene_ids
            .iter()
            .filter_map(|&g| gene_summaries.get(g).and_then(|s| s.knee_um))
            .collect();
        knees.sort_unstable();
        recommended = knees.get(knees.len() / 2).copied();
    }

    let Some(mut rec) = recommended else {
        return (None, "no_recommendation".to_string(), "".to_string(), None);
    };

    let mut status = "ok".to_string();
    let mut note = String::new();

    let frac_at = |w: u32| -> Option<f64> {
        let idx = *width_to_idx.get(&w)?;
        let mut pass = 0u64;
        for &g in informative_gene_ids {
            let r = &rows_by_gene[g][idx];
            if r.n_bins_ge2 >= cfg.b_min {
                pass += 1;
            }
        }
        Some(pass as f64 / informative_gene_ids.len() as f64)
    };

    let mut frac = frac_at(rec);
    if let Some(f) = frac {
        if f < cfg.min_informative_frac_ge2 {
            // Deterministic conservative fallback: choose the smallest larger evaluated width that passes.
            let mut found = None;
            for &w in widths_um {
                if w <= rec {
                    continue;
                }
                if let Some(fw) = frac_at(w) {
                    if fw >= cfg.min_informative_frac_ge2 {
                        found = Some((w, fw));
                        break;
                    }
                }
            }
            match found {
                Some((w2, f2)) => {
                    status = "fallback_increase_width".to_string();
                    note = format!(
                        "sanity_check failed at {}µm (frac_ge2={:.3}); moved to {}µm (frac_ge2={:.3})",
                        rec, f, w2, f2
                    );
                    rec = w2;
                    frac = Some(f2);
                }
                None => {
                    status = "fallback_max_width".to_string();
                    note = format!(
                        "sanity_check failed at {}µm (frac_ge2={:.3}); no larger evaluated width passed; using max width {}µm",
                        rec,
                        f,
                        widths_um.last().copied().unwrap_or(rec)
                    );
                    rec = widths_um.last().copied().unwrap_or(rec);
                    frac = frac_at(rec);
                }
            }
        }
    }

    (Some(rec), status, note, frac)
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
    overall_knee_um: Option<u32>,
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
}

fn write_knee_summary(
    path: PathBuf,
    gene_summaries: &[GeneSummary],
    gene_order: &[usize],
    meta: KneeSummaryMeta,
) -> Result<()> {
    let header = [
        "gene",
        "knee_um",
        "gene_status",
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
    ]
    .join("\t");

    let mut knees: Vec<u32> = gene_summaries
        .iter()
        .filter(|g| g.gene_status == "informative")
        .filter_map(|g| g.knee_um)
        .collect();
    knees.sort_unstable();
    let (knee_mean, knee_sd, knee_min, knee_max) = knee_stats(&knees);

    let summary_row = {
        format!(
            "__SUMMARY__\tNA\tsummary\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}",
            meta.prescreen_global_knee_um,
            meta.overall_knee_um.map(|v| v.to_string()).unwrap_or_else(|| "NA".to_string()),
            meta.recommended_grid_um.map(|v| v.to_string()).unwrap_or_else(|| "NA".to_string()),
            meta.recommendation_frac_ge2
                .map(|v| format!("{:.6}", v))
                .unwrap_or_else(|| "NA".to_string()),
            meta.recommendation_status,
            meta.recommendation_note.replace('\t', " "),
            opt_f64(knee_mean),
            opt_f64(knee_sd),
            knee_min.map(|v| v.to_string()).unwrap_or_else(|| "NA".to_string()),
            knee_max.map(|v| v.to_string()).unwrap_or_else(|| "NA".to_string()),
            meta.n_genes_total,
            meta.n_genes_informative,
            meta.n_genes_excluded_degenerate,
            meta.b_min,
            meta.n_min,
            meta.min_points,
            meta.min_informative_frac_ge2,
        )
    };

    let mut rows: Vec<String> = Vec::new();
    rows.push(summary_row);
    rows.extend(gene_order.iter().map(|&g| {
        let gs = &gene_summaries[g];
        format!(
            "{}\t{}\t{}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA",
            gs.gene,
            gs.knee_um.map(|v| v.to_string()).unwrap_or_else(|| "NA".to_string()),
            gs.gene_status
        )
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
