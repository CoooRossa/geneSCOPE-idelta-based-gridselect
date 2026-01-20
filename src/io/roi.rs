use crate::roi::polygon::{Point, Polygon};
use anyhow::{anyhow, Context, Result};
use serde::Deserialize;
use std::fs;
use std::path::Path;

pub fn read_roi_polygon<P: AsRef<Path>>(path: P) -> Result<Polygon> {
    let path = path.as_ref();
    let raw = fs::read_to_string(path).with_context(|| format!("failed to read: {}", path.display()))?;

    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();

    match ext.as_str() {
        "geojson" | "json" => parse_geojson_polygon(&raw),
        "wkt" | "txt" => parse_wkt_polygon(&raw),
        "csv" | "tsv" => parse_vertices_table(&raw, if ext == "tsv" { b'\t' } else { b',' }),
        _ => {
            // As a fallback, interpret file contents as WKT.
            parse_wkt_polygon(&raw)
        }
    }
}

#[derive(Debug, Deserialize)]
#[serde(untagged)]
enum GeoJsonAny {
    FeatureCollection(GeoJsonFeatureCollection),
    Feature(GeoJsonFeature),
    Geometry(GeoJsonGeometry),
}

#[derive(Debug, Deserialize)]
struct GeoJsonFeatureCollection {
    #[serde(rename = "type")]
    _kind: String,
    features: Vec<GeoJsonFeature>,
}

#[derive(Debug, Deserialize)]
struct GeoJsonFeature {
    #[serde(rename = "type")]
    _kind: String,
    geometry: GeoJsonGeometry,
}

#[derive(Debug, Deserialize)]
struct GeoJsonGeometry {
    #[serde(rename = "type")]
    kind: String,
    coordinates: serde_json::Value,
}

fn parse_geojson_polygon(raw: &str) -> Result<Polygon> {
    let parsed: GeoJsonAny =
        serde_json::from_str(raw).context("failed to parse GeoJSON")?;

    let geom = match parsed {
        GeoJsonAny::Geometry(g) => g,
        GeoJsonAny::Feature(f) => f.geometry,
        GeoJsonAny::FeatureCollection(fc) => fc
            .features
            .into_iter()
            .next()
            .ok_or_else(|| anyhow!("GeoJSON FeatureCollection had no features"))?
            .geometry,
    };

    if geom.kind != "Polygon" {
        return Err(anyhow!("GeoJSON geometry type must be Polygon; got {}", geom.kind));
    }

    // Expect coordinates: [ [ [x,y], ... ] ] (outer ring first)
    let rings: Vec<Vec<[f64; 2]>> = serde_json::from_value(geom.coordinates)
        .context("failed to parse Polygon.coordinates")?;
    let outer = rings
        .into_iter()
        .next()
        .ok_or_else(|| anyhow!("GeoJSON Polygon had no rings"))?;

    let mut pts: Vec<Point> = Vec::with_capacity(outer.len());
    for xy in outer {
        pts.push(Point::from_um_f64(xy[0], xy[1]));
    }
    Ok(Polygon::new_closed(pts)?)
}

fn parse_wkt_polygon(raw: &str) -> Result<Polygon> {
    // Minimal WKT POLYGON parser: POLYGON((x y, x y, ...))
    let raw = raw.trim();
    let raw_upper = raw.to_ascii_uppercase();
    let prefix = "POLYGON((";
    let suffix = "))";
    let start = raw_upper
        .find(prefix)
        .ok_or_else(|| anyhow!("WKT must start with 'POLYGON(('"))?;
    if start != 0 {
        return Err(anyhow!("WKT must start with 'POLYGON(('"));
    }
    if !raw_upper.ends_with(suffix) {
        return Err(anyhow!("WKT must end with '))'"));
    }
    let inner = &raw[prefix.len()..raw.len() - suffix.len()];

    let mut pts: Vec<Point> = Vec::new();
    for (idx, pair) in inner.split(',').enumerate() {
        let pair = pair.trim();
        if pair.is_empty() {
            continue;
        }
        let mut it = pair.split_whitespace();
        let x: f64 = it
            .next()
            .ok_or_else(|| anyhow!("missing x at vertex {}", idx + 1))?
            .parse()
            .with_context(|| format!("invalid x at vertex {}", idx + 1))?;
        let y: f64 = it
            .next()
            .ok_or_else(|| anyhow!("missing y at vertex {}", idx + 1))?
            .parse()
            .with_context(|| format!("invalid y at vertex {}", idx + 1))?;
        pts.push(Point::from_um_f64(x, y));
    }
    Ok(Polygon::new_closed(pts)?)
}

fn parse_vertices_table(raw: &str, delimiter: u8) -> Result<Polygon> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(delimiter)
        .has_headers(true)
        .from_reader(raw.as_bytes());

    let headers = rdr.headers().context("failed reading vertices header")?.clone();
    let idx_x = headers
        .iter()
        .position(|h| h == "x_um" || h == "x")
        .ok_or_else(|| anyhow!("vertices table missing 'x_um' (or 'x')"))?;
    let idx_y = headers
        .iter()
        .position(|h| h == "y_um" || h == "y")
        .ok_or_else(|| anyhow!("vertices table missing 'y_um' (or 'y')"))?;

    let mut pts: Vec<Point> = Vec::new();
    for (row_idx, rec) in rdr.records().enumerate() {
        let rec = rec.with_context(|| format!("failed reading vertex record {}", row_idx + 1))?;
        let x: f64 = rec
            .get(idx_x)
            .ok_or_else(|| anyhow!("missing x at vertex record {}", row_idx + 1))?
            .parse()
            .with_context(|| format!("invalid x at vertex record {}", row_idx + 1))?;
        let y: f64 = rec
            .get(idx_y)
            .ok_or_else(|| anyhow!("missing y at vertex record {}", row_idx + 1))?
            .parse()
            .with_context(|| format!("invalid y at vertex record {}", row_idx + 1))?;
        pts.push(Point::from_um_f64(x, y));
    }
    Ok(Polygon::new_closed(pts)?)
}
