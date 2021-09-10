from itertools import islice

import cyvcf2 as vcf
import numpy as np
import pytest
import sparse
from skbio.stats.ordination import pcoa

import lostruct as ls

# Benchmarks to prevent performance regressions...

# This codeblock taken from https://docs.python.org/3/library/itertools.html
def take(n, iterable):
    "Return first n items of the iterable as a list"
    return list(islice(iterable, n))


# ...

vcf_file = "test_data/chr1-filtered.vcf.gz"


def test_getgts(benchmark):
    record = next(ls.get_snps(vcf_file, "chr1"))
    benchmark(ls.get_gts, record)


error_tolerance = 0.00000001


def test_cov_pca(benchmark):
    windows, _ = ls.parse_vcf(vcf_file, "chr1", 99)
    benchmark(ls.cov_pca, windows[0].todense(), 5, 1)


def test_eigen_windows(benchmark):
    windows, _ = ls.parse_vcf(vcf_file, "chr1", 99)
    benchmark(ls.eigen_windows, windows[0], 5, 1)


def test_l1_norm(benchmark):
    windows, _ = ls.parse_vcf(vcf_file, "chr1", 99)
    _, _, eigenvals, _ = ls.eigen_windows(windows[0], 5, 1)
    benchmark(ls.l1_norm, eigenvals)


@pytest.mark.benchmark(group="PCoA", disable_gc=True, min_rounds=50, warmup=True)
def test_pcoa_default(benchmark):
    windows, _ = ls.parse_vcf(vcf_file, "chr1", 99)
    result = list()
    for x in take(50, windows):
        result.append(ls.eigen_windows(x, 10, 1))
    result = np.vstack(result)
    dists = ls.get_pc_dists(result)
    benchmark(pcoa, dists, number_of_dimensions=10)


@pytest.mark.benchmark(group="PCoA", disable_gc=True, min_rounds=50, warmup=True)
def test_pcoa_fsvd_method(benchmark):
    windows, _ = ls.parse_vcf(vcf_file, "chr1", 99)
    result = list()
    for x in take(50, windows):
        result.append(ls.eigen_windows(x, 10, 1))
    result = np.vstack(result)
    dists = ls.get_pc_dists(result)
    benchmark(pcoa, dists, method="fsvd", number_of_dimensions=10)


@pytest.mark.benchmark(
    group="Get PCs Dists", disable_gc=True, min_rounds=100, warmup=True
)
def test_get_pcs_dists(benchmark):
    windows, _ = ls.parse_vcf(vcf_file, "chr1", 99)
    result = list()
    for x in take(50, windows):
        result.append(ls.eigen_windows(x, 10, 1))
    result = np.vstack(result)
    benchmark(ls.get_pc_dists, result)


@pytest.mark.benchmark(
    group="Get PCs Dists", disable_gc=True, min_rounds=100, warmup=True
)
def test_get_pcs_dists_fastmath(benchmark):
    windows, _ = ls.parse_vcf(vcf_file, "chr1", 99)
    result = list()
    for x in take(20, windows):
        result.append(ls.eigen_windows(x, 10, 1))
    result = np.vstack(result)
    benchmark(ls.get_pc_dists, result, fastmath=True)


@pytest.mark.benchmark(
    group="Get PCs Dists", disable_gc=True, min_rounds=100, warmup=True
)
def test_get_pcs_dists_jax(benchmark):
    windows, _ = ls.parse_vcf(vcf_file, "chr1", 99)
    result = list()
    for x in take(20, windows):
        result.append(ls.eigen_windows(x, 10, 1))
    result = np.vstack(result)
    benchmark(ls.get_pc_dists, result, jax=True)
