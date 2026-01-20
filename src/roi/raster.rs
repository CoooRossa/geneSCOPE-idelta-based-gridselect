use crate::roi::polygon::Polygon;
use anyhow::{anyhow, Result};

#[derive(Debug, Clone)]
pub struct GridSpec {
    pub x0_100: i64,
    pub y0_100: i64,
    pub w0_100: i64,
    pub nx: usize,
    pub ny: usize,
}

impl GridSpec {
    pub fn from_polygon_bbox(poly: &Polygon, w0_um: u32) -> Result<Self> {
        let w0_100 = w0_um as i64 * crate::COORD_SCALE_PER_UM;
        if w0_100 <= 0 {
            return Err(anyhow!("w0 must be >= 1"));
        }

        let (min_x, max_x, min_y, max_y) = poly.bbox_100();
        let x0_100 = min_x.div_euclid(w0_100) * w0_100;
        let y0_100 = min_y.div_euclid(w0_100) * w0_100;

        let span_x = max_x - x0_100;
        let span_y = max_y - y0_100;
        let mut nx = ((span_x + w0_100 - 1) / w0_100) as usize;
        let mut ny = ((span_y + w0_100 - 1) / w0_100) as usize;
        nx = nx.max(1);
        ny = ny.max(1);

        Ok(Self {
            x0_100,
            y0_100,
            w0_100,
            nx,
            ny,
        })
    }

    pub fn cell_center_100(&self, i: usize, j: usize) -> (i64, i64) {
        let x = self.x0_100 + (i as i64) * self.w0_100 + self.w0_100 / 2;
        let y = self.y0_100 + (j as i64) * self.w0_100 + self.w0_100 / 2;
        (x, y)
    }

    pub fn point_to_cell(&self, x_100: i64, y_100: i64) -> Option<(u32, u32)> {
        let dx = x_100 - self.x0_100;
        let dy = y_100 - self.y0_100;
        if dx < 0 || dy < 0 {
            return None;
        }
        let i = dx.div_euclid(self.w0_100) as i64;
        let j = dy.div_euclid(self.w0_100) as i64;
        if i < 0 || j < 0 {
            return None;
        }
        let (i, j) = (i as usize, j as usize);
        if i >= self.nx || j >= self.ny {
            return None;
        }
        Some((i as u32, j as u32))
    }
}

#[derive(Debug, Clone)]
pub struct RoiMask {
    pub grid: GridSpec,
    pub mask: Vec<u8>, // length nx*ny
    pub prefix_sum: Vec<u32>, // (nx+1)*(ny+1) summed-area table
}

impl RoiMask {
    pub fn idx(&self, i: usize, j: usize) -> usize {
        j * self.grid.nx + i
    }

    pub fn cell(&self, i: usize, j: usize) -> u8 {
        self.mask[self.idx(i, j)]
    }

    pub fn build(poly: &Polygon, grid: GridSpec) -> Self {
        // Overlap rule (conservative rasterization):
        // mark a w0 cell as in-ROI if the polygon interior intersects the cell.
        //
        // Practical implementation: for each cell-row, OR three scanline fills within the row
        // (y_low + 1, y_mid, y_high - 1 in 0.01µm units). This avoids per-cell PIP checks and
        // is fast enough for ~5000×5000 grids.
        let mut mask = vec![0u8; grid.nx * grid.ny];

        for j in 0..grid.ny {
            let row_off = j * grid.nx;
            let y_low = grid.y0_100 + (j as i64) * grid.w0_100;
            let y_high = y_low + grid.w0_100;
            let y_mid = y_low + grid.w0_100 / 2;

            // Ensure samples stay inside the open cell to avoid boundary-only intersections.
            let ys = [y_low + 1, y_mid, y_high - 1];

            let origin_x = grid.x0_100 as f64;
            let w0 = grid.w0_100 as f64;
            let center_offset = 0.5 * w0;

            for &y in &ys {
                let xs = poly.scanline_x_intersections(y);
                if xs.len() < 2 {
                    continue;
                }
                for pair in xs.chunks_exact(2) {
                    let (x_left, x_right) = (pair[0], pair[1]);
                    if !(x_right > x_left) {
                        continue;
                    }

                    let i_start = ((x_left - origin_x - center_offset) / w0).ceil() as i64;
                    let i_end = ((x_right - origin_x - center_offset) / w0).ceil() as i64;
                    let i0 = i_start.max(0) as usize;
                    let i1 = (i_end.max(0) as usize).min(grid.nx);
                    for i in i0..i1 {
                        mask[row_off + i] = 1;
                    }
                }
            }
        }

        let prefix_sum = build_prefix_sum_u8(&mask, grid.nx, grid.ny);

        Self {
            grid,
            mask,
            prefix_sum,
        }
    }

    pub fn sum_in_block(&self, i0: usize, j0: usize, i1: usize, j1: usize) -> u32 {
        // Sum over rectangle [i0, i1) x [j0, j1).
        let nx1 = self.grid.nx + 1;
        let a = self.prefix_sum[j0 * nx1 + i0];
        let b = self.prefix_sum[j0 * nx1 + i1];
        let c = self.prefix_sum[j1 * nx1 + i0];
        let d = self.prefix_sum[j1 * nx1 + i1];
        d + a - b - c
    }
}

fn build_prefix_sum_u8(mask: &[u8], nx: usize, ny: usize) -> Vec<u32> {
    let mut ps = vec![0u32; (nx + 1) * (ny + 1)];
    for j in 0..ny {
        let mut row_sum = 0u32;
        for i in 0..nx {
            row_sum += mask[j * nx + i] as u32;
            let above = ps[j * (nx + 1) + (i + 1)];
            ps[(j + 1) * (nx + 1) + (i + 1)] = above + row_sum;
        }
    }
    ps
}
