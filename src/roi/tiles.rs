use crate::roi::raster::RoiMask;
use std::collections::BTreeSet;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct TileId {
    pub tx: u32,
    pub ty: u32,
}

impl TileId {
    pub fn to_u64(self) -> u64 {
        ((self.tx as u64) << 32) | self.ty as u64
    }
}

pub fn list_valid_tiles(roi: &RoiMask, tile_k: usize) -> Vec<TileId> {
    let nx_tiles = (roi.grid.nx + tile_k - 1) / tile_k;
    let ny_tiles = (roi.grid.ny + tile_k - 1) / tile_k;

    let mut out: BTreeSet<TileId> = BTreeSet::new();
    for ty in 0..ny_tiles {
        for tx in 0..nx_tiles {
            let i0 = tx * tile_k;
            let j0 = ty * tile_k;
            let i1 = ((tx + 1) * tile_k).min(roi.grid.nx);
            let j1 = ((ty + 1) * tile_k).min(roi.grid.ny);
            let s = roi.sum_in_block(i0, j0, i1, j1);
            if s > 0 {
                out.insert(TileId {
                    tx: tx as u32,
                    ty: ty as u32,
                });
            }
        }
    }
    out.into_iter().collect()
}
