import lostruct.lostruct as ls
import cyvcf2 as vcf
import sparse
import numpy as np
import itertools
import skbio.stats.ordination
import unittest

# TODO: Rewrite this to make early VCF parsing and PCA stuff one-time run

# This codeblock taken from https://docs.python.org/3/library/itertools.html
def take(n, iterable):
    "Return first n items of the iterable as a list"
    return list(itertools.islice(iterable, n))
#...

vcf_file = "chr1-filtered.vcf.gz"

class TestVcf(unittest.TestCase):

    def test_readvcf(self):
        self.assertEqual(ls.get_landmarks(vcf_file)[0], "chl_Mt")
        self.assertEqual(ls.get_samples(vcf_file)[0], "HM017-I")
        self.assertEqual(len(ls.get_samples(vcf_file)), 50)

    def test_partitionall(self):
        self.assertEqual(len(list(ls.partition_all(3, range(10)))[3]), 1)

    def test_getsnps(self):
        record = next(ls.get_snps(vcf_file, "chr1"))
        self.assertTrue(isinstance(record, vcf.Variant))

    def test_getgts(self):
        record = next(ls.get_snps(vcf_file, "chr1"))
        gts = ls.get_gts(record)
        self.assertEqual(gts.shape, (50,))
        self.assertEqual(gts[0], 2.)
        self.assertEqual(gts[-1], 2.)

    def test_parse_vcf(self):
        windows, positions = ls.parse_vcf(vcf_file, "chr1", 99, ls.Window.SNP)
        self.assertEqual(len(windows), 119)
        self.assertEqual(len(positions), len(windows))
        self.assertEqual(positions[0][0], 59864)
        self.assertTrue(isinstance(windows[0], sparse.COO))

class TestCalculations(unittest.TestCase):

    error_tolerance = 0.00000001

    def test_cov_pca(self):
        windows, _ = ls.parse_vcf(vcf_file, "chr1", 99)
        covmat, total_variance, eigenvals, eigenvecs = ls.cov_pca(windows[0].todense(), 5, 1)

        # Given our VCF test dataset VCF file, compare that the matrices are
        # the same as a successful run performed of 24 Aug 2020.
        # The values here are the result of np.abs np.sum of the outputs of cov_pca
        # Any significant deviation means something has messed up somewhere else
        self.assertAlmostEqual(np.sum(covmat), -0.4784263654778492, delta=self.error_tolerance)
        self.assertAlmostEqual(np.sum(total_variance), 0.9265612493057297, delta=self.error_tolerance)
        self.assertAlmostEqual(np.sum(eigenvals), 1.735862014813605, delta=self.error_tolerance)
        self.assertAlmostEqual(np.sum(eigenvecs), 0.13157175919284625, delta=self.error_tolerance)

    def test_eigen_windows(self):
        windows, _ = ls.parse_vcf(vcf_file, "chr1", 99)
        covmat, total_variance, eigenvals, eigenvecs = ls.eigen_windows(windows[0], 5, 1)
        self.assertAlmostEqual(np.sum(covmat), -0.4784263654778492, delta=self.error_tolerance)
        self.assertAlmostEqual(np.sum(total_variance), 0.9265612493057297, delta=self.error_tolerance)
        self.assertAlmostEqual(np.sum(eigenvals), 1.735862014813605, delta=self.error_tolerance)
        self.assertAlmostEqual(np.sum(eigenvecs), 0.13157175919284625, delta=self.error_tolerance)

    def test_l1_norm(self):
        windows, _ = ls.parse_vcf(vcf_file, "chr1", 99)
        _, _, eigenvals, _ = ls.eigen_windows(windows[0], 5, 1)
        self.assertEqual(np.sum(ls.l1_norm(eigenvals)), 5.0)

    # Many of these tests are redundant, so this one won't be...
    # It also tests dist_sq_from_pcs so we won't test that separately...
    def test_get_pcs_dists(self):
        windows, _ = ls.parse_vcf(vcf_file, "chr1", 99)
        result = list()
        for x in take(4, windows):
            result.append(ls.eigen_windows(x, 10, 1))
        result = np.vstack(result)
        pc_dists = ls.get_pc_dists(result)

        self.assertEqual(pc_dists[0][0], 0.0)
        self.assertAlmostEqual(pc_dists[0][3], 0.30474948474286145, delta=self.error_tolerance)

    def test_get_pcs_dists_fastmath(self):
        windows, _ = ls.parse_vcf(vcf_file, "chr1", 99)
        result = list()
        for x in take(4, windows):
            result.append(ls.eigen_windows(x, 10, 1))
        result = np.vstack(result)
        pc_dists = ls.get_pc_dists(result, fastmath=True)

        self.assertEqual(pc_dists[0][0], 0.0)
        self.assertAlmostEqual(pc_dists[0][3], 0.30474948474286145, delta=self.error_tolerance)

    def test_compare_to_rcode(self):
        windows, _ = ls.parse_vcf(vcf_file, "chr1", 95)
        covmat, total_variance, eigenvals, eigenvecs = ls.cov_pca(windows[0].todense(), 10, 1)

        results = np.loadtxt("lostruct-results/chr1.filtered.pca.csv", 
                                delimiter=",", 
                                skiprows=1)

        totalandvalsR = results[0][0:11]
        totalandvalsPy = np.concatenate(([total_variance], eigenvals)),
        # Comes out as 0.9999921929150888
        self.assertTrue(np.corrcoef(totalandvalsR, totalandvalsPy)[0][1] >= 0.99999)

        # Squared here, because signs are often opposite between the two analyses.
        eigenvecsR = np.square(results[0][11:61])
        eigenvecsPy = np.square(eigenvecs[0])
        # Comes out as 0.9999921929150888
        self.assertTrue(np.corrcoef(eigenvecsR, eigenvecsPy)[0][1] >= 0.99999)
        self.assertEqual(covmat.shape, (50, 50))

        mds_coords = np.loadtxt("lostruct-results/mds_coords.csv", 
                delimiter=",", skiprows=1, usecols=[2])

        result = list()
        for x in windows:
            result.append(ls.eigen_windows(x, 10, 1))
        result = np.vstack(result)
        pc_dists = ls.get_pc_dists(result)
        mds = skbio.stats.ordination.pcoa(pc_dists)
        # Comes out as 0.9971509982243156
        self.assertTrue(np.corrcoef(mds.samples['PC1'], mds_coords)[0][1] >= 0.995)

        pc_dists = ls.get_pc_dists(result, fastmath=True)
        mds = skbio.stats.ordination.pcoa(pc_dists)
        # Comes out as 0.9971509982243156
        self.assertTrue(np.corrcoef(mds.samples['PC1'], mds_coords)[0][1] >= 0.995)

