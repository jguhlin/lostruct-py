# lostruct-py

This is a reimplementation of lostruct from the original code: [lostruct](https://github.com/petrelharp/local_pca),
by [Joseph Guhlin](https://github.com/jguhlin) with assistance by [Peter Ralph](https://github.com/petrelharp).

[![Build Status](https://travis-ci.org/jguhlin/lostruct-py.svg?branch=master)](https://travis-ci.org/jguhlin/lostruct-py)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3997106.svg)](https://doi.org/10.5281/zenodo.3997106)[![codecov](https://codecov.io/gh/jguhlin/lostruct-py/branch/master/graph/badge.svg)](https://codecov.io/gh/jguhlin/lostruct-py)

# Demonstration
Please see the [Example Notebook](https://nbviewer.jupyter.org/github/jguhlin/lostruct-py/blob/master/Lostruct%20Example.ipynb)

# Installation
Lostruct-py is available on [PyPi](https://pypi.org/project/lostruct/)
```pip install lostruct``` is the easiest way to get started.

# Usage

## Input Files
Inputs should be a set of markers in BCF or VCF format. Both should be indexed as appropriate (see: bcftools index). Filtering before running this analysis is strongly suggested (Allele frequency, SNPs only, missingness, etc).

# Citing
If you use this version, plesae cite it via Zenodo [DOI: 10.5281/zenodo.3997106](https://doi.org/10.5281/zenodo.3997106)
as well as the original paper describing the method:
```
Li, Han, and Peter Ralph. "Local PCA shows how the effect of population structure differs along the genome." Genetics 211.1 (2019): 289-304.
```

## CyVCF2
This project also uses [cyvcf2](https://github.com/brentp/cyvcf2) for fast VCF processing and should be cited:

```
Brent S Pedersen, Aaron R Quinlan, cyvcf2: fast, flexible variant analysis with Python, Bioinformatics, Volume 33, Issue 12, 15 June 2017, Pages 1867–1869, https://doi.org/10.1093/bioinformatics/btx057
```

# Changes from Lostruct R package
Please note numpy and R are different when it comes to row-major vs. column-major. Essentially, many things in the python version will be transposed from R.

# Requirements
Python >= 3.6 (may work with older versions). Developed on Python 3.8.5

* numba
* numpy
* cyvcf2

CyVCF2 requires zlib-dev, libbz2-dev, libcurl-dev, liblzma-dev; numa requires libllvm.
These may be installed with `conda` or `pip`, e.g. by running `pip install -r requirements.txt`.

# Changes
See [CHANGES.MD](CHANGES.md) for the full list.

## 0.0.4
* Package name changed to lostruct
* Parallelization of get_pc_dists
* Implementation of fastmath parameter for get_pc_dists

# Tests

Tests were derived from [Medicago HapMap data](http://www.medicagohapmap.org/downloads/Mt40/Mt4.0_HapMap_README.pdf). While the software had high correlation with lostruct R the values were determined. If values begin to deviate from the method these tests will now fail.

To run tests simply do:
```
python -m nose
```

The tests furthermore require `unittest` and `scikit-bio` (and, `nose` to run them this way).

## TOX
Tox allows you run tests with multiple versions of the python interpreter in venvs. It is best to use pyenv to install multiple versions python to run before submitting pull requests to be certain tests complete successfully across all versions.

# Correlation Data
To test correlation of results between the R and Python versions we used data from the Medicago HapMap project, specifically SNPs for sister taxa chromsome 1, processed, and run with LoStruct R.

## Data
```bcftools annotate chr1-filtered-set-2014Apr15.bcf -x INFO,FORMAT | bcftools view -a -i 'F_MISSING<=0.2' | bcftools view -q 0.05 -q 0.95 -m2 -M2 -a -Oz -o chr1-filtered.vcf.gz```

## Lostruct Processing
```Rscript run_lostruct.R -t SNP -s 95 -k 10 -m 10 -i data/```

Run 21 Aug 2020, using lostruct R git hash: 444b8c64bebdf7cdd0323e7735ccadddfc1c8989

This generates the mds_coords.tsv that is used in the correlation comparison. Additionally, the existing tests cover correlation.

To test the weight generation, a random sample of weights was created and used. Output was moved to `lostruct-results/weights_mds_coords.csv` and generated with the random_weights.txt found in `test_data/random_weights.txt`.

```./run_lostruct.R -i data -t snp -s 95 -k 10 -m 10 -w random_weights.txt```

# FAQ / Notes

## Future
Currently the end-user is expected to save the outputs. But would be good to save it in a similar way to lostruct R-code. Please open an issue if you need this.

## Feature Completeness with R implementation
We are not yet feature complete with the R implementation. If something is needed please check for existing issues and comment about your need.

## PCA, MDS, PCoA
PCoA returns the same results as lostruct's MDS implementation (cmdscale). In the example Jupyter notebook you can see the correlation is R =~ 0.998. Some examples of other methods of clustering / looking at differences are included in the notebook.

## Speed and Memory
NUMBA and CyVCF2 are used for speeding up processes, and the software becomes multithreaded by default. The Sparse library is used to reduce memory requirements. parse_vcf function is multithreaded. Distance calculation is not.

### tl;dr of below
Below two options are offered, fastmath for get_pc_dists function, and method="fsvd" for pcoa. When using both you will see a performance increase and memory requirement decrease. Accuracy should decrease, but the absolute correlation we see with our test dataset remains ~0.998. Be aware when using fsvd the sign of the correlation may change.

### Fastmath
Additionally, a mode implemented Numba's "fastmath" is available. For the function get_pc_dists() set fastmath=True. This results in a ~8% speed boost with very little change in the final output (correlation to R code output remains >= 0.995). This was benchmarked on the Medicago data used in the jupyter notebook using timeit, with 100 repeats with fastmath=False and Fastmath=True.
```get_pc_dists(result, fastmath=True)```

The difference with fastmath=True and leaving it off can be seen [here](https://raw.githubusercontent.com/jguhlin/lostruct-py/master/benchmark_0.0.4-Get_PCs_Dists.svg). Note: Downloading the file will allow you to see more detailed information, as some javascript is contained in the SVG but disabled on GitHub.

If you need to limit thread usage, please see [Numba's guide](https://numba.pydata.org/numba-doc/latest/user/threading-layer.html#setting-the-number-of-threads)

### Very Large Datasets
The R implementation handles very large datasets in less memory. The problem arises with the PCoA function. A metric MDS using sklearn may work. Another alternative would be to export the data and run cmdscale in R directly.

The sklearn MDS function differs from the scikit-bio function, here we focus on the scikit-bio version.

There are two options in python for this as well:
```pcoa(method="fsvd", ...)```
Which reduces memory and increases speed, at the cost of some accuracy.

```pcoa(inplace=True, ...)```
Centers a distance matrix in-place, further reducing memory requirements.

```pcoa(number_of_dimensions=10)```
Returns only the first 10 dimensions (configurable) of the scaling. This has no real effect if method is default or manuially set to "eigh" as the eigenvalues and eigenvectors are all calculated, so all are calculated and this becomes a truncation.

You can see the difference between method="fsvd" and method="eigh" (default) [here](https://raw.githubusercontent.com/jguhlin/lostruct-py/master/benchmark_0.0.4-PCoA.svg). These are tested with a minimum of 50 rounds. Note: Downloading the file will allow you to see more detailed information, as some javascript is contained in the SVG but disabled on GitHub.

Using all three techniques, correlation is maintained although the sign may change.
```
mds = pcoa(pc_dists, method="fsvd", inplace=True, number_of_dimensions=10)
np.corrcoef(mds.samples["PC1"], mds_coords['MDS1'].to_numpy())[0][1]
-0.9978147088087447
```

For more information please see the [applicable documentation](http://scikit-bio.org/docs/0.5.4/generated/generated/skbio.stats.ordination.pcoa.html) as well as the [relevant changelog](https://github.com/biocore/scikit-bio/blob/master/CHANGELOG.md#version-054-2018-08-23). A [Zenodo entry](https://zenodo.org/record/1404403) is also available on this topic.

# References

Additional citations can be found in [CITATIONS](Citations.md) (UMAP, PHATE, Medicago HapMap).