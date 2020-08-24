# lostruct-py

This is a reimplementation of lostruct from the original code: [Lostruct](https://github.com/petrelharp/local_pca). Please cite the original paper

[![Build Status](https://travis-ci.org/jguhlin/lostruct-py.svg?branch=master)](https://travis-ci.org/jguhlin/lostruct-py) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3997106.svg)](https://doi.org/10.5281/zenodo.3997106)


# Demonstration / How to use
Please see the [Example Notebook](https://nbviewer.jupyter.org/github/jguhlin/lostruct-py/blob/master/Lostruct-py%20Example.ipynb)

# Installation
Lostruct-py is available on [PyPi](https://pypi.org/project/lostruct-py/)
```pip install lostruct-py``` is the easiest way to get started.

# Citing
Please use our DOI to cite this specific project. Also please cite the original Lostruct paper and CyVCF2. 

DOI: 10.5281/zenodo.3997106

## Original Lostruct Paper
Please cite the original lostruct paper:
```
Li, Han, and Peter Ralph. "Local PCA shows how the effect of population structure differs along the genome." Genetics 211.1 (2019): 289-304.
```

## CyVCF2
This paper also uses [cyvcf2](https://github.com/brentp/cyvcf2) for fast VCF processing and should be cited:

```
Brent S Pedersen, Aaron R Quinlan, cyvcf2: fast, flexible variant analysis with Python, Bioinformatics, Volume 33, Issue 12, 15 June 2017, Pages 1867–1869, https://doi.org/10.1093/bioinformatics/btx057
```

# Changes from Lostruct R package
Please note numpy and R are different when it comes to row-major vs. column-major. Essentially, many things in the python version will be transposed from R.

# Requirements
Python >= 3.6 (may work with older versions). Developed on Python 3.8.5

* numba
* numpy
* pandas
* scipy
* skbio
* sklearn
* cyvcf2

CyVCF2 requires zlib-dev, libbz2-dev, libcurl-dev, liblzma-dev, and probably others

Easiest to install all of these through conda

# Correlation Data
Used Medicago HapMap sister taxa chromsome 1, processed, and run with LoStruct

## Data
```bcftools annotate chr1-filtered-set-2014Apr15.bcf -x INFO,FORMAT | bcftools view -a -i 'F_MISSING<=0.2' | bcftools view -q 0.05 -q 0.95 -m2 -M2 -a -Oz -o chr1-filtered.vcf.gz```

## Lostruct Processing
```Rscript run_lostruct.R -t SNP -s 95 -k 10 -m 10 -i data/```

This generates the mds_coords.tsv that is used in the correlation comparison.

# FAQ / Notes

## Future
Currently the end-user is expected to save the outputs. But could be good to save it in a similar way to lostruct R-code. Please open an issue if you need this.

## PCA, MDS, PCoA
PCoA returns the same results as lostruct's MDS implementation (cmdscale). In the example Jupyter notebook you can see the correlation is R =~ 0.998. Some examples of other methods of clustering / looking at differences are included in the notebook.

## Casting complex values to real discards the imaginary part
This is fine.

