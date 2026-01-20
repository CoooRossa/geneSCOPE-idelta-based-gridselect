use crate::COORD_SCALE_PER_UM;
use anyhow::{anyhow, Result};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Point {
    pub x_100: i64,
    pub y_100: i64,
}

impl Point {
    pub fn from_um_f64(x_um: f64, y_um: f64) -> Self {
        Self {
            x_100: (x_um * COORD_SCALE_PER_UM as f64).round() as i64,
            y_100: (y_um * COORD_SCALE_PER_UM as f64).round() as i64,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Polygon {
    pub vertices: Vec<Point>,
}

impl Polygon {
    pub fn new_closed(mut vertices: Vec<Point>) -> Result<Self> {
        if vertices.len() < 3 {
            return Err(anyhow!("polygon needs >= 3 vertices"));
        }
        if vertices.first() != vertices.last() {
            let first = *vertices.first().unwrap();
            vertices.push(first);
        }
        Ok(Self { vertices })
    }

    pub fn bbox_100(&self) -> (i64, i64, i64, i64) {
        let mut min_x = i64::MAX;
        let mut max_x = i64::MIN;
        let mut min_y = i64::MAX;
        let mut max_y = i64::MIN;
        for p in &self.vertices {
            min_x = min_x.min(p.x_100);
            max_x = max_x.max(p.x_100);
            min_y = min_y.min(p.y_100);
            max_y = max_y.max(p.y_100);
        }
        (min_x, max_x, min_y, max_y)
    }

    /// Even-odd point-in-polygon. Boundary counts as inside.
    pub fn contains_point(&self, x_100: i64, y_100: i64) -> bool {
        // Ray casting in +x direction (even-odd rule).
        let mut inside = false;
        let n = self.vertices.len();
        if n < 4 {
            return false;
        }

        for i in 0..(n - 1) {
            let a = self.vertices[i];
            let b = self.vertices[i + 1];

            // Check if point is exactly on segment.
            if point_on_segment(x_100, y_100, a, b) {
                return true;
            }

            let (x1, y1) = (a.x_100, a.y_100);
            let (x2, y2) = (b.x_100, b.y_100);

            // Only consider edges that straddle the horizontal line through (x, y).
            // Half-open rule avoids double-counting vertices.
            if (y1 <= y_100 && y2 > y_100) || (y2 <= y_100 && y1 > y_100) {
                let y = y_100 as f64;
                let x_int = x1 as f64 + (y - y1 as f64) * (x2 - x1) as f64 / (y2 - y1) as f64;
                if x_int > x_100 as f64 {
                    inside = !inside;
                }
            }
        }
        inside
    }

    /// Sorted x-intersections (0.01µm units) of polygon boundary with horizontal scanline at `y_100`.
    pub fn scanline_x_intersections(&self, y_100: i64) -> Vec<f64> {
        let n = self.vertices.len();
        let mut xs: Vec<f64> = Vec::new();
        if n < 4 {
            return xs;
        }

        let y = y_100 as f64;
        for i in 0..(n - 1) {
            let a = self.vertices[i];
            let b = self.vertices[i + 1];
            let (x1, y1) = (a.x_100 as f64, a.y_100 as f64);
            let (x2, y2) = (b.x_100 as f64, b.y_100 as f64);
            if (y2 - y1).abs() < f64::EPSILON {
                continue; // horizontal edge
            }
            // Half-open rule: include lower endpoint, exclude upper.
            if (y1 <= y && y2 > y) || (y2 <= y && y1 > y) {
                let t = (y - y1) / (y2 - y1);
                let x = x1 + t * (x2 - x1);
                xs.push(x);
            }
        }

        xs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        xs
    }
}

fn point_on_segment(px: i64, py: i64, a: Point, b: Point) -> bool {
    let (x1, y1) = (a.x_100, a.y_100);
    let (x2, y2) = (b.x_100, b.y_100);

    let cross = (px - x1) as i128 * (y2 - y1) as i128 - (py - y1) as i128 * (x2 - x1) as i128;
    if cross != 0 {
        return false;
    }

    let dot = (px - x1) as i128 * (px - x2) as i128 + (py - y1) as i128 * (py - y2) as i128;
    dot <= 0
}
