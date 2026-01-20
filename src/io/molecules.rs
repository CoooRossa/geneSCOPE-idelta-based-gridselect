use crate::COORD_SCALE_PER_UM;
use anyhow::{anyhow, Context, Result};
use std::fs::File;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct RawMolecule {
    pub gene: String,
    pub x_100: i64,
    pub y_100: i64,
}

pub fn read_molecules<P: AsRef<Path>>(path: P) -> Result<Vec<RawMolecule>> {
    let path = path.as_ref();
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();

    match ext.as_str() {
        "csv" => read_molecules_delimited(path, b','),
        "tsv" | "txt" => read_molecules_delimited(path, b'\t'),
        "parquet" | "pq" => read_molecules_parquet(path),
        _ => Err(anyhow!(
            "unsupported molecules extension '{}': {}",
            ext,
            path.display()
        )),
    }
}

fn read_molecules_delimited(path: &Path, delimiter: u8) -> Result<Vec<RawMolecule>> {
    let file = File::open(path).with_context(|| format!("failed to open: {}", path.display()))?;
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(delimiter)
        .has_headers(true)
        .from_reader(file);

    let headers = rdr
        .headers()
        .context("failed reading header row")?
        .clone();

    let idx_gene = headers
        .iter()
        .position(|h| h == "gene")
        .ok_or_else(|| anyhow!("missing required column 'gene'"))?;
    let idx_x = headers
        .iter()
        .position(|h| h == "x_um" || h == "x")
        .ok_or_else(|| anyhow!("missing required column 'x_um' (or 'x')"))?;
    let idx_y = headers
        .iter()
        .position(|h| h == "y_um" || h == "y")
        .ok_or_else(|| anyhow!("missing required column 'y_um' (or 'y')"))?;

    let mut out: Vec<RawMolecule> = Vec::new();
    for (row_idx, rec) in rdr.records().enumerate() {
        let rec = rec.with_context(|| format!("failed reading record {}", row_idx + 1))?;
        let gene = rec
            .get(idx_gene)
            .ok_or_else(|| anyhow!("missing gene at record {}", row_idx + 1))?
            .to_string();
        let x_um: f64 = rec
            .get(idx_x)
            .ok_or_else(|| anyhow!("missing x at record {}", row_idx + 1))?
            .parse()
            .with_context(|| format!("invalid x at record {}", row_idx + 1))?;
        let y_um: f64 = rec
            .get(idx_y)
            .ok_or_else(|| anyhow!("missing y at record {}", row_idx + 1))?
            .parse()
            .with_context(|| format!("invalid y at record {}", row_idx + 1))?;

        out.push(RawMolecule {
            gene,
            x_100: (x_um * COORD_SCALE_PER_UM as f64).round() as i64,
            y_100: (y_um * COORD_SCALE_PER_UM as f64).round() as i64,
        });
    }
    Ok(out)
}

fn read_molecules_parquet(path: &Path) -> Result<Vec<RawMolecule>> {
    use parquet::file::reader::{FileReader, SerializedFileReader};

    let file = File::open(path).with_context(|| format!("failed to open: {}", path.display()))?;
    let reader =
        SerializedFileReader::new(file).context("failed creating parquet reader")?;

    let mut iter = reader
        .get_row_iter(None)
        .context("failed creating parquet row iterator")?;

    let mut idx_gene: Option<usize> = None;
    let mut idx_x: Option<usize> = None;
    let mut idx_y: Option<usize> = None;

    let mut out: Vec<RawMolecule> = Vec::new();
    while let Some(row) = iter.next() {
        let row = row.context("failed reading parquet row")?;

        if idx_gene.is_none() || idx_x.is_none() || idx_y.is_none() {
            for (i, (name, _field)) in row.get_column_iter().enumerate() {
                match name.as_str() {
                    "gene" | "gene_name" if idx_gene.is_none() => idx_gene = Some(i),
                    "x_um" | "x" if idx_x.is_none() => idx_x = Some(i),
                    "y_um" | "y" if idx_y.is_none() => idx_y = Some(i),
                    _ => {}
                }
            }
        }

        let idx_gene = idx_gene.ok_or_else(|| anyhow!("missing required column 'gene'"))?;
        let idx_x = idx_x.ok_or_else(|| anyhow!("missing required column 'x_um' (or 'x')"))?;
        let idx_y = idx_y.ok_or_else(|| anyhow!("missing required column 'y_um' (or 'y')"))?;

        let gene = row_get_string(&row, idx_gene)?;
        let x_um = row_get_f64(&row, idx_x)?;
        let y_um = row_get_f64(&row, idx_y)?;

        out.push(RawMolecule {
            gene,
            x_100: (x_um * COORD_SCALE_PER_UM as f64).round() as i64,
            y_100: (y_um * COORD_SCALE_PER_UM as f64).round() as i64,
        });
    }
    Ok(out)
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
