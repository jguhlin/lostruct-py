# Release 0.0.4
* Split up get_pc_dists fn so numba can parallelize the majority of the work properly
* Rename python package to lostruct
* Add docstrings for most functions
* Added notes for working with large datasets
* Added more information to tutorial / jupyter notebook demonstration
* Added argument support for SNPs by window or by base pair length -- Not yet implemented
* Added support for weights - Not yet tested!
* Support for benchmarks
* Support for fastmath Numba for some functions, along with testing and correlation comparisons
* Experimentation with pcoa method="fsvd"

# Release 0.0.3
Minor typo fix for PyPi

# Release 0.0.2
* Drop pandas requirement, switch cov matrix to numpy masked array
* Implement sparse matrix support to reduce memory requirements on very large datasets
* Changes suggested by @petrelharp
* Unit testing
* Travis CI support
* Better integration into pypi
* requirements.txt

# Release 0.0.1
