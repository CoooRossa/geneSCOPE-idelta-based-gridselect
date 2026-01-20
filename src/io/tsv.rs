use anyhow::{Context, Result};
use flate2::write::GzEncoder;
use flate2::Compression;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

pub fn write_tsv<P: AsRef<Path>>(path: P, header: &str, rows: impl IntoIterator<Item = String>) -> Result<()> {
    let file = File::create(&path).with_context(|| format!("failed to create: {}", path.as_ref().display()))?;
    let mut w = BufWriter::new(file);
    writeln!(w, "{header}")?;
    for row in rows {
        writeln!(w, "{row}")?;
    }
    Ok(())
}

pub fn write_tsv_gz<P: AsRef<Path>>(path: P, header: &str, rows: impl IntoIterator<Item = String>) -> Result<()> {
    let file = File::create(&path).with_context(|| format!("failed to create: {}", path.as_ref().display()))?;
    let gz = GzEncoder::new(file, Compression::default());
    let mut w = BufWriter::new(gz);
    writeln!(w, "{header}")?;
    for row in rows {
        writeln!(w, "{row}")?;
    }
    Ok(())
}
