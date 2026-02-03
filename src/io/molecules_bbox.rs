use anyhow::{anyhow, Context, Result};
use flate2::read::GzDecoder;
use std::fs::{self, File};
use std::io::{self, Read};
use std::path::Path;

#[derive(Debug, Clone, Copy)]
struct BBox {
    min_x: f64,
    max_x: f64,
    min_y: f64,
    max_y: f64,
}

impl BBox {
    fn contains(&self, x: f64, y: f64) -> bool {
        x >= self.min_x && x <= self.max_x && y >= self.min_y && y <= self.max_y
    }
}

pub const DEFAULT_EXCLUDE_PREFIXES: [&str; 8] = [
    "Unassigned",
    "NegControl",
    "Background",
    "DeprecatedCodeword",
    "SystemControl",
    "Negative",
    "BlankCodeword",
    "Blank",
];

#[derive(Debug, Clone)]
pub struct MakeMoleculesBBoxConfig {
    pub roi_csv: String,
    pub transcripts_path: String,
    pub out_tsv: String,
    pub gene_col: String,
    pub x_col: String,
    pub y_col: String,
    pub exclude_prefixes: Vec<String>,
}

fn read_roi_bbox(roi_csv: &Path) -> Result<BBox> {
    let raw = fs::read_to_string(roi_csv)
        .with_context(|| format!("failed reading roi_csv: {}", roi_csv.display()))?;
    let mut cleaned: Vec<&str> = Vec::new();
    for ln in raw.lines() {
        let t = ln.trim();
        if t.is_empty() {
            continue;
        }
        if t.trim_start().starts_with('#') {
            continue;
        }
        cleaned.push(ln);
    }
    if cleaned.len() < 2 {
        return Err(anyhow!(
            "ROI file had no data rows after stripping comments: {}",
            roi_csv.display()
        ));
    }

    let header = cleaned[0];
    let delim = if header.contains('\t') && !header.contains(',') {
        b'\t'
    } else {
        b','
    };

    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(delim)
        .has_headers(true)
        .from_reader(io::Cursor::new(cleaned.join("\n")));
    let headers = rdr.headers().context("failed reading ROI header")?.clone();

    let headers_low: Vec<String> = headers.iter().map(|h| h.trim().to_ascii_lowercase()).collect();
    let idx_x = headers_low
        .iter()
        .position(|h| h == "x_um" || h == "x")
        .ok_or_else(|| anyhow!("ROI header must include x_um (or x)"))?;
    let idx_y = headers_low
        .iter()
        .position(|h| h == "y_um" || h == "y")
        .ok_or_else(|| anyhow!("ROI header must include y_um (or y)"))?;

    let mut min_x = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    let mut n = 0usize;
    for rec in rdr.records() {
        let rec = rec.context("failed reading ROI row")?;
        let x: f64 = rec
            .get(idx_x)
            .unwrap_or("")
            .parse()
            .with_context(|| format!("invalid ROI x value: {}", rec.get(idx_x).unwrap_or("")))?;
        let y: f64 = rec
            .get(idx_y)
            .unwrap_or("")
            .parse()
            .with_context(|| format!("invalid ROI y value: {}", rec.get(idx_y).unwrap_or("")))?;
        if !x.is_finite() || !y.is_finite() {
            return Err(anyhow!("ROI vertices contain non-finite x/y values"));
        }
        min_x = min_x.min(x);
        max_x = max_x.max(x);
        min_y = min_y.min(y);
        max_y = max_y.max(y);
        n += 1;
    }
    if n < 3 || !(min_x < max_x && min_y < max_y) {
        return Err(anyhow!(
            "ROI bbox invalid (need >=3 vertices and non-degenerate bounds): {}",
            roi_csv.display()
        ));
    }
    Ok(BBox {
        min_x,
        max_x,
        min_y,
        max_y,
    })
}

fn open_maybe_gz(path: &Path) -> Result<Box<dyn Read>> {
    let file = File::open(path).with_context(|| format!("failed to open: {}", path.display()))?;
    let path_low = path.to_string_lossy().to_ascii_lowercase();
    if path_low.ends_with(".gz") {
        Ok(Box::new(GzDecoder::new(file)))
    } else {
        Ok(Box::new(file))
    }
}

fn infer_transcripts_kind_and_delim(path: &Path) -> Result<(TranscriptsKind, u8)> {
    let p = path.to_string_lossy().to_ascii_lowercase();
    if p.ends_with(".parquet") || p.ends_with(".pq") {
        return Ok((TranscriptsKind::Parquet, b','));
    }
    if p.ends_with(".csv") || p.ends_with(".csv.gz") {
        return Ok((TranscriptsKind::Delimited, b','));
    }
    if p.ends_with(".tsv") || p.ends_with(".tsv.gz") || p.ends_with(".txt") || p.ends_with(".txt.gz") {
        return Ok((TranscriptsKind::Delimited, b'\t'));
    }
    Err(anyhow!(
        "unsupported transcripts extension: {} (expected .csv(.gz), .tsv(.gz), .parquet)",
        path.display()
    ))
}

#[derive(Debug, Clone, Copy)]
enum TranscriptsKind {
    Delimited,
    Parquet,
}

fn should_exclude_gene(gene: &str, prefixes: &[String]) -> bool {
    prefixes.iter().any(|p| !p.is_empty() && gene.starts_with(p))
}

fn make_from_delimited(
    bbox: BBox,
    transcripts: &Path,
    out_tsv: &Path,
    delim: u8,
    gene_col: &str,
    x_col: &str,
    y_col: &str,
    exclude_prefixes: &[String],
) -> Result<(u64, u64, u64)> {
    let reader = open_maybe_gz(transcripts)?;
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(delim)
        .has_headers(true)
        .from_reader(reader);

    let headers = rdr
        .headers()
        .with_context(|| format!("failed reading transcripts header: {}", transcripts.display()))?
        .clone();
    let idx_gene = headers
        .iter()
        .position(|h| h == gene_col)
        .ok_or_else(|| anyhow!("missing transcripts column '{}'", gene_col))?;
    let idx_x = headers
        .iter()
        .position(|h| h == x_col)
        .ok_or_else(|| anyhow!("missing transcripts column '{}'", x_col))?;
    let idx_y = headers
        .iter()
        .position(|h| h == y_col)
        .ok_or_else(|| anyhow!("missing transcripts column '{}'", y_col))?;

    if let Some(parent) = out_tsv.parent() {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed creating out dir: {}", parent.display()))?;
    }
    let fout = File::create(out_tsv)
        .with_context(|| format!("failed creating out_tsv: {}", out_tsv.display()))?;
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(fout);
    wtr.write_record(["gene", "x_um", "y_um"])
        .context("failed writing molecules header")?;

    let mut n_in: u64 = 0;
    let mut n_out: u64 = 0;
    let mut n_excl_prefix: u64 = 0;

    for rec in rdr.records() {
        let rec = rec.context("failed reading transcripts row")?;
        n_in += 1;
        let gene = rec.get(idx_gene).unwrap_or("").to_string();
        if gene.is_empty() {
            continue;
        }
        if should_exclude_gene(&gene, exclude_prefixes) {
            n_excl_prefix += 1;
            continue;
        }
        let x: f64 = match rec.get(idx_x).unwrap_or("").parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let y: f64 = match rec.get(idx_y).unwrap_or("").parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        if !x.is_finite() || !y.is_finite() {
            continue;
        }
        if !bbox.contains(x, y) {
            continue;
        }
        let x_s = x.to_string();
        let y_s = y.to_string();
        wtr.write_record([gene.as_str(), x_s.as_str(), y_s.as_str()])
            .context("failed writing molecules row")?;
        n_out += 1;
    }

    wtr.flush().context("failed flushing molecules TSV")?;
    Ok((n_in, n_out, n_excl_prefix))
}

fn row_get_string(row: &parquet::record::Row, idx: usize) -> Result<String> {
    use parquet::record::RowAccessor;
    if let Ok(s) = row.get_string(idx) {
        return Ok(s.clone());
    }
    if let Ok(b) = row.get_bytes(idx) {
        return Ok(b.as_utf8()?.to_string());
    }
    Err(anyhow!("cannot decode parquet string column at index {}", idx))
}

fn row_get_f64(row: &parquet::record::Row, idx: usize) -> Result<f64> {
    use parquet::record::RowAccessor;
    if let Ok(v) = row.get_double(idx) {
        return Ok(v);
    }
    if let Ok(v) = row.get_float(idx) {
        return Ok(v as f64);
    }
    if let Ok(v) = row.get_int(idx) {
        return Ok(v as f64);
    }
    if let Ok(v) = row.get_long(idx) {
        return Ok(v as f64);
    }
    if let Ok(v) = row.get_uint(idx) {
        return Ok(v as f64);
    }
    if let Ok(v) = row.get_ulong(idx) {
        return Ok(v as f64);
    }
    Err(anyhow!("cannot decode parquet numeric column at index {}", idx))
}

fn make_from_parquet(
    bbox: BBox,
    transcripts: &Path,
    out_tsv: &Path,
    gene_col: &str,
    x_col: &str,
    y_col: &str,
    exclude_prefixes: &[String],
) -> Result<(u64, u64, u64)> {
    use parquet::file::reader::{FileReader, SerializedFileReader};

    let file = File::open(transcripts)
        .with_context(|| format!("failed to open: {}", transcripts.display()))?;
    let reader = SerializedFileReader::new(file)
        .with_context(|| format!("failed creating parquet reader: {}", transcripts.display()))?;
    let mut iter = reader
        .get_row_iter(None)
        .context("failed creating parquet row iterator")?;

    if let Some(parent) = out_tsv.parent() {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed creating out dir: {}", parent.display()))?;
    }
    let fout = File::create(out_tsv)
        .with_context(|| format!("failed creating out_tsv: {}", out_tsv.display()))?;
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(fout);
    wtr.write_record(["gene", "x_um", "y_um"])
        .context("failed writing molecules header")?;

    let gene_col_l = gene_col.to_ascii_lowercase();
    let x_col_l = x_col.to_ascii_lowercase();
    let y_col_l = y_col.to_ascii_lowercase();

    let mut idx_gene: Option<usize> = None;
    let mut idx_x: Option<usize> = None;
    let mut idx_y: Option<usize> = None;

    let mut n_in: u64 = 0;
    let mut n_out: u64 = 0;
    let mut n_excl_prefix: u64 = 0;

    while let Some(row) = iter.next() {
        let row = row.context("failed reading parquet row")?;
        n_in += 1;

        if idx_gene.is_none() || idx_x.is_none() || idx_y.is_none() {
            for (i, (name, _field)) in row.get_column_iter().enumerate() {
                let name_l = name.to_ascii_lowercase();
                if idx_gene.is_none() && name_l == gene_col_l {
                    idx_gene = Some(i);
                } else if idx_x.is_none() && name_l == x_col_l {
                    idx_x = Some(i);
                } else if idx_y.is_none() && name_l == y_col_l {
                    idx_y = Some(i);
                }
            }
        }

        let idx_gene = idx_gene.ok_or_else(|| anyhow!("missing transcripts column '{}'", gene_col))?;
        let idx_x = idx_x.ok_or_else(|| anyhow!("missing transcripts column '{}'", x_col))?;
        let idx_y = idx_y.ok_or_else(|| anyhow!("missing transcripts column '{}'", y_col))?;

        let gene = row_get_string(&row, idx_gene)?;
        if gene.is_empty() {
            continue;
        }
        if should_exclude_gene(&gene, exclude_prefixes) {
            n_excl_prefix += 1;
            continue;
        }

        let x = row_get_f64(&row, idx_x)?;
        let y = row_get_f64(&row, idx_y)?;
        if !x.is_finite() || !y.is_finite() {
            continue;
        }
        if !bbox.contains(x, y) {
            continue;
        }

        let x_s = x.to_string();
        let y_s = y.to_string();
        wtr.write_record([gene.as_str(), x_s.as_str(), y_s.as_str()])
            .context("failed writing molecules row")?;
        n_out += 1;
    }

    wtr.flush().context("failed flushing molecules TSV")?;
    Ok((n_in, n_out, n_excl_prefix))
}

pub fn make_molecules_bbox(cfg: MakeMoleculesBBoxConfig) -> Result<()> {
    let roi_csv = Path::new(&cfg.roi_csv);
    let transcripts = Path::new(&cfg.transcripts_path);
    let out_tsv = Path::new(&cfg.out_tsv);

    let bbox = read_roi_bbox(roi_csv)?;
    eprintln!(
        "[make-molecules-bbox] ROI bbox: x=[{}, {}] y=[{}, {}]",
        bbox.min_x, bbox.max_x, bbox.min_y, bbox.max_y
    );

    if !cfg.exclude_prefixes.is_empty() {
        eprintln!(
            "[make-molecules-bbox] exclude_prefixes={}",
            cfg.exclude_prefixes.join(",")
        );
    }

    let (kind, delim) = infer_transcripts_kind_and_delim(transcripts)?;
    let (n_in, n_out, n_excl_prefix) = match kind {
        TranscriptsKind::Delimited => make_from_delimited(
            bbox,
            transcripts,
            out_tsv,
            delim,
            &cfg.gene_col,
            &cfg.x_col,
            &cfg.y_col,
            &cfg.exclude_prefixes,
        )?,
        TranscriptsKind::Parquet => make_from_parquet(
            bbox,
            transcripts,
            out_tsv,
            &cfg.gene_col,
            &cfg.x_col,
            &cfg.y_col,
            &cfg.exclude_prefixes,
        )?,
    };

    eprintln!(
        "[make-molecules-bbox] wrote {} / {} transcripts rows to {} (excluded_by_prefix={})",
        n_out,
        n_in,
        out_tsv.display(),
        n_excl_prefix
    );
    Ok(())
}
