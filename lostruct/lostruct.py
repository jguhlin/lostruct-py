from itertools import islice
import numpy as np
from math import sqrt
from numba import jit
from cyvcf2 import VCF
import sparse
from enum import Enum

class Window(Enum):
    SNP = 1
    BP = 2

def get_landmarks(vcf_file):
    """ 
    Get landmarks (chromosomes, scaffolds, contigs, etc) from a VCF/BCF file.
    """
    return VCF(vcf_file).seqnames

def get_samples(vcf_file):
    """
    Get the samples from a VCF/BCF file. This is the order the data will remain in as well.
    """
    return VCF(vcf_file).samples

def partition_all(n, iterable):
    """
    Utility function to partition an iterable to equal sized chunks, except for the last chunk which may be smaller than the chunk size.
    """
    x = iter(iterable)
    while True:
        y = list(islice(x, n))
        if not y:
            break
        yield y

# Split by chromosome
def get_snps(file, landmark):
    """
    Get SNPs from CyVCF2 for a specific landmark. gts012=True. VCF/BCF files should be filtered
    """
    vcf_reader = VCF(file, gts012=True)(landmark)
    for record in vcf_reader:
        yield record
        # Let users handle this...
        #if record.is_snp:
        #    yield record
        
def get_gts(x):
    """
    Convert to numeric format and convert missing to np.nan
    """
    gts = x.gt_types
    gts = gts.astype(np.float64)
    gts[gts == 3.] = np.nan
    return gts

def parse_vcf(vcf_file, landmark, window_size, window_type=Window.SNP):
    """
    Parse a VCF (or BCF) file. Must be properly indexed (see: bcftools index command)

    Parameters:
    vcf_file (String): Location of VCF/BCF file, relative or absolute path.

    landmark (String): Identifier or Landmark (chromosome, contig, scaffold) exactly as in VCF file.

    window_size (int): Size of marker windows in SNPs or BP

    window_type (Enum): Set to lostruct.Window.SNP for window size set as SNPs, and lostruct.Window.BP for window size as base pairs.*
        * Currently, only window sizes of SNPs are supported, and this option is ignored.

    Returns:
    (list) Windows of SNPs, as sparse matrices, with missing data set to np.nan
    (list) Lists of Positions of SNPs found in Windows

    """
    snp_generator = partition_all(window_size, get_snps(vcf_file, landmark))
    windows = list()
    positions = list()
    for window in snp_generator:
        # Should mask it instead of skipping it...
        if (len(window) < window_size):
            break
        gts = [get_gts(i) for i in window]
        pos = [y.POS for y in window]
        windows.append(sparse.COO(np.asarray(gts)))
        positions.append(pos)

    return windows, positions

# Input files should be bgzip vcfs with tabix idx, bcf with idx should probably work too

def cov_pca(snps,k,w=1):
    """
    Returns the covariance matrix, total variance, eigenvalues, eigenvectors for a given SNP window. Returns top k components.

    Please note this differs from the R function by also returning the covariance matrix.

    Parameters:
    snps (np.array): SNPs as output from parse_vcf or get_gts function
    
    k (int): Number of components to return

    w (np.array, list, or tuple): Weights for samples. Works, but untested at this time. Please see: https://github.com/jguhlin/lostruct-py/issues/7

    Returns:
    covmat: Covariance matrix
    total_variance: List of the total variance captured by each component
    eigenvalues
    eigenvectors
    """
    n = len(snps)
    rowmeans = np.nanmean(snps, axis=1)
    rowmeans = np.reshape(rowmeans, (n, 1))
    subtracted = np.array(snps - rowmeans, dtype=np.float64)
    #covmat = np.cov(subtracted) # Not using weights, so ignoring that part from lostruct
    #covmat = pd.DataFrame(subtracted).cov()
    covmat = np.ma.cov(np.ma.array(subtracted, mask=np.isnan(subtracted)), rowvar=False)
    #covmat = np.ma.cov(subtracted)
    if np.all(w != 1):
        #sqrt_w = np.repeat(np.sqrt(w), covmat.shape[0])
        sqrt_w = np.sqrt(w)
        covmat = np.multiply(covmat, sqrt_w.T)
        covmat = np.multiply(covmat, sqrt_w)

    # This is the first returned argument for cov_pca in R
    total_variance = np.sum(np.power(covmat, 2).flatten())
    #first_return_arg = np.sum(np.power(covmat, 2).flatten())
    vals, vectors = np.linalg.eig(covmat)

    vals = vals[:k].real.astype(np.float64)
    vectors = vectors[:,:k].T.real.astype(np.float64)
    
    return covmat, total_variance, vals, vectors

# I think this became a wrapper function for cov_pca...
def eigen_windows(snps, k, w):
    """
    See cov_pca. This is mostly a placeholder but does handle conversion of snps back to a dense matrix.
    """
    return cov_pca(snps.todense(), k, w)

def l1_norm(eigenvals):
    """
    Applies l1 norm to a set of eigenvalues.
    """
    z = np.vstack(eigenvals)
    norm = np.linalg.norm(z, ord=1, axis=1, keepdims=True)
    return z / norm

# k = number of primary components
# w = weight to apply
# norm = normalization to apply (L1 or L2)

def get_pc_dists(windows):
    """
    Calculate distances between window matrices.

    Works on only the upper triangle of the matrix, but clones the data into the lower half as well.
    """
    n = len(windows)
    vals = l1_norm([x[2] for x in windows])
    comparison = np.zeros((n, n), dtype=np.float64)
    upper_triangle = zip(*np.triu_indices(n, k=1))
    
    vals = vals.real.astype(np.float64)

    for i,j in upper_triangle:
        comparison[i,j] = dist_sq_from_pcs(vals[i], vals[j], windows[i][3], windows[j][3])

    # Make symmetric
    comparison = comparison + comparison.T
    # Remove negatives...
    comparison[comparison < 0] = 0

    # Get square root
    return np.sqrt(comparison)

@jit(nopython=True, parallel=True)
def dist_sq_from_pcs(v1, v2, xvec, yvec):
    """
    Given two matrices, calculates the distance between them.
    """
    Xt = np.dot(xvec, yvec.T)
    return np.square(v1).sum() + np.square(v2).sum() - 2 * np.sum(v1*np.sum(Xt*(v2*Xt), axis=1))
