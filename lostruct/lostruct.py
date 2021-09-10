"""Python implementation of Lostruct, which calculates local population structure of a genome"""

from enum import Enum
from itertools import islice

import numpy as np
import sparse
from cyvcf2 import VCF
from numba import jit
import jax.numpy as jnp
import jax
from multiprocessing import Pool
from jax.config import config
from jax import pmap, vmap

config.update("jax_enable_x64", True)


class Window(Enum):
    """Species whether to treat window sizes as per-SNP or per-BP... Currently not implemented"""

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
    Utility function to partition an iterable to equal sized chunks, except for the
    last chunk which may be smaller than the chunk size.
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
    Get SNPs from CyVCF2 for a specific landmark. gts012=True. VCF/BCF files should be
    filtered before passing into this function.
    """
    vcf_reader = VCF(file, gts012=True)(landmark)
    for record in vcf_reader:
        yield record
        # Let users handle this...
        # if record.is_snp:
        #    yield record


def get_gts(x):
    """
    Convert to numeric format and convert missing to np.nan
    """
    gts = x.gt_types
    gts = gts.astype(np.float64)
    gts[gts == 3.0] = np.nan
    return gts


# Input files should be bgzip vcfs with tabix idx, bcf with idx should probably work too
def parse_vcf(vcf_file, landmark, window_size, window_type=Window.SNP):
    """
    Parse a VCF (or BCF) file. Must be properly indexed (see: bcftools index command)

    Parameters:
    vcf_file (String): Location of VCF/BCF file, relative or absolute path.

    landmark (String): Identifier or Landmark (chromosome, contig, scaffold) exactly
        as in VCF file.

    window_size (int): Size of marker windows in SNPs or BP

    window_type (Enum): Set to lostruct.Window.SNP for window size set as SNPs, and
        lostruct.Window.BP for window size as base pairs.*
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
        if len(window) < window_size:
            break
        gts = [get_gts(i) for i in window]
        pos = [y.POS for y in window]
        windows.append(sparse.COO(np.asarray(gts)))
        positions.append(pos)

    return windows, positions


def cov_pca(snps, k, w=1):
    """
    Returns the covariance matrix, total variance, eigenvalues, eigenvectors for a given
    SNP window. Returns top k components.

    Please note this differs from the R function by also returning the covariance matrix.

    Parameters:
    snps (np.array): SNPs as output from parse_vcf or get_gts function

    k (int): Number of components to return

    w (np.array, list, or tuple): Weights for samples. Works, but untested at this time.
    Please see: https://github.com/jguhlin/lostruct-py/issues/7

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
    # covmat = np.cov(subtracted) # Not using weights, so ignoring that part from lostruct
    # covmat = pd.DataFrame(subtracted).cov()
    covmat = np.ma.cov(np.ma.array(subtracted, mask=np.isnan(subtracted)), rowvar=False)
    # covmat = np.ma.cov(subtracted)
    if np.all(w != 1):
        # sqrt_w = np.repeat(np.sqrt(w), covmat.shape[0])
        sqrt_w = np.sqrt(w)
        covmat = np.multiply(covmat, sqrt_w.T)
        covmat = np.multiply(covmat, sqrt_w)

    # This is the first returned argument for cov_pca in R
    total_variance = np.sum(np.power(covmat, 2).flatten())
    # first_return_arg = np.sum(np.power(covmat, 2).flatten())
    vals, vectors = np.linalg.eig(covmat)

    vals = vals[:k].real.astype(np.float64)
    vectors = vectors[:, :k].T.real.astype(np.float64)

    return covmat, total_variance, vals, vectors


# I think this became a wrapper function for cov_pca...
def eigen_windows(snps, k, w):
    """
    See cov_pca. This is mostly a placeholder but does handle conversion of snps back
    to a dense matrix.
    """
    return cov_pca(snps.todense(), k, w)


# @jit(parallel=True)
def l1_norm(eigenvals):
    """
    Applies l1 norm to a set of eigenvalues.
    """
    z = np.vstack(eigenvals)
    # Because of the norm fn below, we can't use numba
    # Possible to redo the data, but not sure the benefit would be worth it
    norm = np.linalg.norm(z, ord=1, axis=1, keepdims=True)
    return z / norm


@jax.jit
def l1_norm_jax(eigenvals):
    """
    Applies l1 norm to a set of eigenvalues using JAX
    """
    z = jnp.vstack(eigenvals)
    # Because of the norm fn below, we can't use numba
    # Possible to redo the data, but not sure the benefit would be worth it
    norm = jnp.linalg.norm(z, ord=1, axis=1, keepdims=True)
    return z / norm


@jit(nopython=True, parallel=True, fastmath=True)
def calc_dists_fastmath(vals, eigenvecs):
    """
    Calculate the distances of windows given a set of eigenvectors and normalized
    eigenvalues. Called from get_pc_dist() function.

    This is separated out to take advantage of Numba optimizations.
    """
    #comparison = np.zeros((n, n), dtype=np.float64)
    #upper_triangle = zip(*np.triu_indices(n, k=1))

    comparison = np.zeros((vals.shape[0], vals.shape[0]), dtype=np.float64)
    upper_triangle = np.stack(np.triu_indices_from(comparison, k=1), axis=1)

    vals = vals.real.astype(np.float64)

    for i, j in upper_triangle:
        comparison[i, j] = dist_sq_from_pcs_fastmath(
            vals[i], vals[j], eigenvecs[i], eigenvecs[j]
        )

    # Make symmetric
    comparison = comparison + comparison.T

    return comparison


@jit(nopython=True, parallel=True)
def calc_dists(vals, eigenvecs):
    """
    Calculate the distances of windows given a set of eigenvectors and normalized
    eigenvalues. Called from get_pc_dist() function.

    This is separated out to take advantage of Numba optimizations.
    """
    #comparison = np.zeros((n, n), dtype=np.float64)
    comparison = np.zeros((vals.shape[0], vals.shape[0]), dtype=np.float64)
    #upper_triangle = zip(*np.triu_indices(n, k=1))
    upper_triangle = np.stack(np.triu_indices_from(comparison, k=1), axis=1)

    vals = vals.real.astype(np.float64)

    for i,j in upper_triangle:
        comparison[i, j] = dist_sq_from_pcs(
            vals[i], vals[j], eigenvecs[i], eigenvecs[j]
        )

    # Make symmetric
    comparison = comparison + comparison.T

    return comparison


# k = number of primary components
# w = weight to apply
# norm = normalization to apply (L1 or L2)
def get_pc_dists(windows, fastmath=False, jax=False, w=1):
    """
    Calculate distances between window matrices.

    Works on only the upper triangle of the matrix, but clones the data into the lower
    half as well.
    """
    #n = len(windows)

    # Remove negatives... Can't be placed within Numba code
    # Get square root

    if fastmath:
        vals = np.asarray(l1_norm_jax(np.asarray([x[2] for x in windows])))
        vals = vals.real.astype(np.float64)
        comparison = calc_dists_fastmath(vals, np.asarray([x[3] for x in windows]))
        return np.sqrt(np.where(comparison > 0, comparison, 0))
    elif jax:
        vals = l1_norm_jax(jnp.asarray([x[2] for x in windows]))
        comparison = calc_dists_jax(vals, jnp.asarray([x[3] for x in windows]))
        return jnp.sqrt(jnp.where(comparison > 0, comparison, 0))
    else:
        vals = l1_norm(np.asarray([x[2] for x in windows]))
        vals = vals.real.astype(np.float64)
        comparison = calc_dists(vals, np.asarray([x[3] for x in windows]))
        return np.sqrt(np.where(comparison > 0, comparison, 0))


@jit(nopython=True, parallel=True, fastmath=True)
def dist_sq_from_pcs_fastmath(v1, v2, xvec, yvec):
    """
    Given two matrices, calculates the distance between them.
    """
    Xt = np.dot(xvec, yvec.T)
    return (
        np.square(v1).sum()
        + np.square(v2).sum()
        - 2 * np.sum(v1 * np.sum(Xt * (v2 * Xt), axis=1))
    )


@jit(nopython=True, parallel=True, fastmath=False)
def dist_sq_from_pcs(v1, v2, xvec, yvec):
    """
    Given two matrices, calculates the distance between them.
    """
    Xt = np.dot(xvec, yvec.T)
    return (
        np.square(v1).sum()
        + np.square(v2).sum()
        - 2 * np.sum(v1 * np.sum(Xt * (v2 * Xt), axis=1))
    )


# @jax.jit
def calc_dists_jax(vals, eigenvecs):
    """
    Calculate the distances of windows given a set of eigenvectors and normalized
    eigenvalues. Called from get_pc_dist(..., jax=True) function.
    """
    comparison = jnp.zeros((vals.shape[0], vals.shape[0]), dtype=jnp.float64)

    upper_triangle = jnp.triu_indices_from(comparison, k=1)
    # upper_triangle = jnp.stack(upper_triangle, axis=1)

    upper_triangle_vals = jnp.stack(
        (
            jnp.take(vals, upper_triangle[0], axis=0),
            jnp.take(vals, upper_triangle[1], axis=0),
        ),
        axis=1,
    )

    upper_triangle_eigenvecs = jnp.stack(
        (
            jnp.take(eigenvecs, upper_triangle[0], axis=0),
            jnp.take(eigenvecs, upper_triangle[1], axis=0),
        ),
        axis=1,
    )

    #    upper_triangle = jnp.stack((upper_triangle_vals, upper_triangle_eigenvecs), axis=1)

    #    def do_calc(i, j):
    #        return dist_sq_from_pcs_jax(vals[i], vals[j], eigenvecs[i], eigenvecs[j])

    #    with Pool() as pool:
    #        comparison_out = pool.starmap(lambda vs, evs: dist_sq_from_pcs_jax(vs[0], vs[1], evs[0], evs[1]), upper_triangle_vals, upper_triangle_eigenvecs)
    #        comparison_out = pool.map(do_calc, upper_triangle)

    #    for i, j in upper_triangle:
    #        comparison = comparison.at[i, j].set(dist_sq_from_pcs_jax(
    #            vals[i], vals[j], eigenvecs[i], eigenvecs[j]
    #        ))

    print("Got here!")

    pfn = jax.vmap(lambda vs, evs: dist_sq_from_pcs_jax(vs[0], vs[1], evs[0], evs[1]))
    comparison_out = pfn(upper_triangle_vals, upper_triangle_eigenvecs)
    comparison = comparison.at[jnp.triu_indices_from(comparison, k=1)].set(comparison_out)

    # Make symmetric
    return comparison + comparison.T


@jax.jit
def dist_sq_from_pcs_jax(v1, v2, xvec, yvec):
    """
    Given two matrices, calculates the distance between them.
    """
    Xt = jnp.dot(xvec, yvec.T)
    return (
        jnp.square(v1).sum()
        + jnp.square(v2).sum()
        - 2 * jnp.sum(v1 * jnp.sum(Xt * (v2 * Xt), axis=1))
    )
