from itertools import islice
import numpy as np
import operator as op
import pandas as pd
from math import sqrt
import scipy.linalg
from sklearn.preprocessing import normalize
from sklearn.manifold import MDS
from skbio.stats.ordination import pcoa
from numba import jit
from cyvcf2 import VCF

# Need to cite:
# https://github.com/brentp/cyvcf2

# newer, cyvcf2 version
def get_landmarks(vcf_file):
    return VCF(vcf_file).seqnames

def get_samples(vcf_file):
    return VCF(vcf_file).samples

def partition_all(n, iterable):
    x = iter(iterable)
    while True:
        y = list(islice(x, n))
        if not y:
            break
        yield y

# Split by chromosome
def get_snps(file, landmark, min_af = 0.0, max_af = 1.0):
    vcf_reader = VCF(file, gts012=True)(landmark)
    for record in vcf_reader:
        yield record
        # Let users handle this...
        #if record.is_snp:
        #    yield record
        
def get_gts(x):
    gts = x.gt_types
    gts = gts.astype(np.float64)
    gts[gts == 3.] = np.nan
    return gts

def parse_vcf(vcf_file, landmark, window_size):
    snp_generator = partition_all(window_size, get_snps(vcf_file, landmark))
    windows = list()
    positions = list()
    for window in snp_generator:
        # Should mask it instead of missing it...
        if (len(window) < window_size):
            break
        gts = [get_gts(i) for i in window]
        pos = [y.POS for y in window]
        windows.append(gts)
        positions.append(pos)

    return windows, positions

# Effectively, this is lostruct-py

# Input files should be bgzip vcfs with tabix idx, bcf with idx should probably work too

# DONE! cov_pca matches the output from lostruct now!
def cov_pca(snps,w,k):
    n = len(snps)
    rowmeans = np.nanmean(snps, axis=1)
    rowmeans = np.reshape(rowmeans, (n, 1))
    subtracted = np.array(snps - rowmeans, dtype=np.float64)
    #covmat = np.cov(subtracted) # Not using weights, so ignoring that part from lostruct
    covmat = pd.DataFrame(subtracted).cov()
    #covmat = np.ma.cov(np.ma.array(subtracted, mask=np.isnan(subtracted)))
    #covmat = np.ma.cov(subtracted)
    sqrt_w = np.repeat(sqrt(w), covmat.shape[0])

    # This is the first returned argument for cov_pca in R
    first_return_arg = np.sum(np.power(covmat, 2).values.flatten())
    #first_return_arg = np.sum(np.power(covmat, 2).flatten())
    vals, vectors = np.linalg.eig(covmat)

    eigenvecs = list()
    for i in range(k):
        eigenvecs.append(vectors[:,i])

    eigenvals = vals[0:k]
    #eigenvals.astype(np.float64)
    
    return covmat, first_return_arg, eigenvals, np.asarray(eigenvecs, dtype=np.float64)

def eigen_windows(snps, k):
    return cov_pca(snps, 1, k)

# TODO: Implement L2 norm (and others?)
def l1_norm(eigenvals):
    z = np.vstack(eigenvals)
    norm = np.linalg.norm(z, ord=1, axis=1, keepdims=True)
    return z / norm

# k = number of primary components
# w = weight to apply
# norm = normalization to apply (L1 or L2)

# return covmat, first_return_arg, eigenvals, eigenvecs
def get_pc_dists(windows):
    n = len(windows)
    vals = l1_norm([x[2] for x in windows])
    comparison = np.zeros((n, n), dtype=np.float64)
    upper_triangle = zip(*np.triu_indices(n, k=1))
    
    vals = vals.astype(np.float64)

    for i,j in upper_triangle:
        comparison[i,j] = dist_sq_from_pcs(vals[i], vals[j], windows[i][3], windows[j][3])

    # Make symmetric
    comparison = comparison + comparison.T
    # Remove negatives...
    comparison[comparison < 0] = 0

    # Get square root
    return np.sqrt(comparison)

@jit(nopython=True, parallel=False)
def dist_sq_from_pcs(v1, v2, xvec, yvec):
    Xt = np.dot(xvec, yvec.T)
    return np.square(v1).sum() + np.square(v2).sum() - 2 * np.sum(v1*np.sum(Xt*(v2*Xt), axis=1))