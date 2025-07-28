from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="eidos-cosmology",
    version="1.0.0",
    author="Benjamin P. Duncan",
    author_email="founder@aerohalo.io",
    description="Emergent Information Dynamics of Observational Spacetime",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/AeroHalo/EIDOS",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20",
        "scipy>=1.7",
        "matplotlib>=3.4",
        "astropy>=5.0",
        "emcee>=3.1",
        "corner>=2.2",
        "pandas>=1.3",
        "h5py>=3.0",
        "tqdm>=4.62",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
            "black>=21.0",
            "flake8>=3.9",
            "sphinx>=4.0",
            "nbsphinx>=0.8",
            "jupyter>=1.0",
            "ipython>=7.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "eidos-analyze=scripts.reproduce_main_results:main",
            "eidos-validate=scripts.run_validation:main",
        ],
    },
)
