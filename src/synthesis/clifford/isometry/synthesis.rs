use crate::structures::CliffordCircuit;
use crate::structures::IsometryTableau;
use crate::structures::Metric;

use super::count::isometry_count_synthesis;

pub fn isometry_synthesis(
    isometry: &IsometryTableau,
    metric: &Metric,
    niter: usize,
) -> CliffordCircuit {
    match metric {
        Metric::COUNT => {
            return isometry_count_synthesis(isometry, niter);
        }
        _ => panic!("Isometry synthesis targeting depth is not available (yet!)"),
    }
}
