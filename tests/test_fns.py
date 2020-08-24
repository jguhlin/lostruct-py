import lostructpy as ls
import cyvcf2 as vcf
import sparse
import numpy as np
from itertools import islice
from skbio.stats.ordination import pcoa

# TODO: Rewrite this to make early VCF parsing and PCA stuff one-time run

# This codeblock taken from https://docs.python.org/3/library/itertools.html
def take(n, iterable):
    "Return first n items of the iterable as a list"
    return list(islice(iterable, n))
#...

vcf_file = "chr1-filtered.vcf.gz"

def test_readvcf():
    assert(ls.get_landmarks(vcf_file)[0] == "chl_Mt")
    assert(ls.get_samples(vcf_file)[0] == "HM017-I")
    assert(len(ls.get_samples(vcf_file)) == 50)

def test_partitionall():
    assert(len(list(ls.partition_all(3, range(10)))[3]) == 1)

def test_getsnps():
    record = next(ls.get_snps(vcf_file, "chr1"))
    assert(isinstance(record, vcf.Variant))

def test_getgts():
    record = next(ls.get_snps(vcf_file, "chr1"))
    gts = ls.get_gts(record)
    assert(gts.shape == (50,))
    assert(gts[0] == 2.)
    assert(gts[-1] == 2.)

def test_parse_vcf():
    windows, positions = ls.parse_vcf(vcf_file, "chr1", 99)
    assert(len(windows) == 119)
    assert(len(positions) == len(windows))
    assert(positions[0][0] == 59864)
    assert(isinstance(windows[0], sparse.COO))

def test_cov_pca():
    windows, _ = ls.parse_vcf(vcf_file, "chr1", 99)
    covmat, total_variance, eigenvals, eigenvecs = ls.cov_pca(windows[0].todense(), 5, 1)
    assert(np.sum(covmat) == -0.4784263654778492)
    assert(np.sum(total_variance) == 0.9265612493057297)
    assert(np.sum(eigenvals) == 1.735862014813605)
    assert(np.sum(eigenvecs) == 0.13157175919284625)

def test_eigen_windows():
    windows, _ = ls.parse_vcf(vcf_file, "chr1", 99)
    covmat, total_variance, eigenvals, eigenvecs = ls.eigen_windows(windows[0], 5, 1)
    assert(np.sum(covmat) == -0.4784263654778492)
    assert(np.sum(total_variance) == 0.9265612493057297)
    assert(np.sum(eigenvals) == 1.735862014813605)
    assert(np.sum(eigenvecs) == 0.13157175919284625)

def test_l1_norm():
    windows, _ = ls.parse_vcf(vcf_file, "chr1", 99)
    _, _, eigenvals, _ = ls.eigen_windows(windows[0], 5, 1)
    assert(np.sum(ls.l1_norm(eigenvals)) == 5.0)

# Many of these tests are redundant, so this one won't be...
# It also tests dist_sq_from_pcs so we won't test that separately...
def test_get_pcs_dists():
    windows, _ = ls.parse_vcf(vcf_file, "chr1", 99)
    result = list()
    for x in take(4, windows):
        result.append(ls.eigen_windows(x, 10, 1))
    result = np.vstack(result)
    pc_dists = ls.get_pc_dists(result)

    assert(pc_dists[0][0] == 0.0)
    assert(pc_dists[0][3] == 0.30474948474286145)

def test_compare_to_rcode():
    windows, _ = ls.parse_vcf(vcf_file, "chr1", 95)
    covmat, total_variance, eigenvals, eigenvecs = ls.cov_pca(windows[0].todense(), 10, 1)

    results = np.loadtxt("lostruct-results/chr1.filtered.pca.csv", 
                            delimiter=",", 
                            skiprows=1)

    totalandvalsR = results[0][0:11]
    totalandvalsPy = np.concatenate(([total_variance], eigenvals)),
    # Comes out as 0.9999921929150888
    assert(np.corrcoef(totalandvalsR, totalandvalsPy)[0][1] >= 0.99999)

    # Squared here, because signs are often opposite between the two analyses.
    eigenvecsR = np.square(results[0][11:61])
    eigenvecsPy = np.square(eigenvecs[0])
    # Comes out as 0.9999921929150888
    assert(np.corrcoef(eigenvecsR, eigenvecsPy)[0][1] >= 0.99999)
    assert(covmat.shape == (50, 50))

    mds_coords = np.loadtxt("lostruct-results/mds_coords.csv", 
            delimiter=",", skiprows=1, usecols=[2])

    result = list()
    for x in windows:
        result.append(ls.eigen_windows(x, 10, 1))
    result = np.vstack(result)
    pc_dists = ls.get_pc_dists(result)
    mds = pcoa(pc_dists)
    # Comes out as 0.9971509982243156
    assert(np.corrcoef(mds.samples['PC1'], mds_coords)[0][1] >= 0.995)
