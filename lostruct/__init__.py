"""Python implementation of Lostruct, which calculates local population structure of a genome"""

from .lostruct import (
    cov_pca,
    dist_sq_from_pcs,
    eigen_windows,
    get_gts,
    get_landmarks,
    get_pc_dists,
    get_samples,
    get_snps,
    l1_norm,
    parse_vcf,
    partition_all,
)
