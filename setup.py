import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="lostruct-py-jguhlin",
    version="0.0.1",
    author="Joseph Guhlin",
    author_email="joseph.guhlin@gmail.com",
    description="Re-implementation of lostruct in Python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jguhlin/lostruct-py",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)