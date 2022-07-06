import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="lostruct",
    version="0.0.4",
    author="Joseph Guhlin",
    author_email="joseph.guhlin@gmail.com",
    description="Re-implementation of lostruct in Python, used to compare local population structure across populations.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jguhlin/lostruct-py",
    packages=setuptools.find_packages(),
    keywords="structure, population, genomics, bioinformatics, PCA, MDS",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Python Software Foundation License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Development Status :: 4 - Beta",
    ],
    install_requires=[
        "numpy",
        "numba",
        "cyvcf2",
        "sparse",
    ],
    extras_require = {
        "jax": ["jax", "jaxlib"]

    },
    python_requires=">=3.6",
    project_urls={
        "Bug Reports": "https://github.com/jguhlin/lostruct-py/issues",
        "Source": "https://github.com/jguhlin/lostruct-py",
        "R-lang version": "https://github.com/petrelharp/local_pca",
        "Original Method Paper": "https://www.genetics.org/content/211/1/289",
    },
)
