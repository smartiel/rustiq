use crate::structures::CliffordCircuit;
use crate::structures::IsometryTableau;
use crate::structures::Metric;

use super::count::isometry_count_synthesis;

pub fn isometry_synthesis(
    isometry: &IsometryTableau,
    metric: &Metric,
    niter: usize,
) -> CliffordCircuit {
    return isometry_count_synthesis(isometry, niter);
}
