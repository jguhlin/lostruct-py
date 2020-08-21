# lostruct-py

This is a reimplementation of lostruct from the original code: [Lostruct](https://github.com/petrelharp/local_pca). Please cite the original paper

## Citing

### Original Lostruct Paper
Please cite the original lostruct paper:
```
 Li, Han, and Peter Ralph. "Local PCA shows how the effect of population structure differs along the genome." Genetics 211.1 (2019): 289-304.
```

### CyVCF2
This paper also uses [cyvcf2](https://github.com/brentp/cyvcf2) for fast VCF processing and should be cited:

```
Brent S Pedersen, Aaron R Quinlan, cyvcf2: fast, flexible variant analysis with Python, Bioinformatics, Volume 33, Issue 12, 15 June 2017, Pages 1867–1869, https://doi.org/10.1093/bioinformatics/btx057
```

## Requirements

* numba
* numpy
* pandas
* scipy
* skbio
* sklearn
* cyvcf2

CyVCF2 requires zlib-dev, libbz2-dev, libcurl-dev, liblzma-dev, and probably others

Easiest to install all of these through conda

# Demonstration / How to use
Please see the [Example Notebook](https://nbviewer.jupyter.org/github/jguhlin/lostruct-py/blob/master/Lostruct-py%20Example.ipynb)

## Data
```bcftools annotate chr1-filtered-set-2014Apr15.bcf -x INFO,FORMAT | bcftools view -a -i 'F_MISSING<=0.2' | bcftools view -q 0.05 -q 0.95 -m2 -M2 -a -Oz -o chr1-filtered.vcf.gz```

# Lostruct / Lostruct-py comparison
```Rscript run_lostruct.R -t SNP -s 95 -k 10 -m 10 -i data/```

# FAQ / Notes

## PCA, MDS, PCoA
PCoA returns the same results as lostruct's MDS implementation (cmdscale). In the example Jupyter notebook you can see the correlation is R =~ 0.998. Some examples of other methods of clustering / looking at differences are included in the notebook.

## Casting complex values to real discards the imaginary part
This is fine.

