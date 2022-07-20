# Benchmarks
Install pytest-benchmark pygal pygaljs to perform benchmarks.

Should be run with:
```pytest --benchmark-histogram=benchmark_version --benchmark-save=version```

Replacing version with a numeric version or a temporary hold. For example:
```pytest --benchmark-histogram=0.0.4 --benchmark-save=0.0.```

Comparisons can be done with:
```pytest-benchmark compare 0.0.4 newversion```
Where 0.0.4 is the last published version to compare against, and newversion is the same output you saved with earlier. See [the docs](https://pytest-benchmark.readthedocs.io/en/latest/comparing.html) for more information.
