import numpy as np
from setuptools import setup, Extension, find_packages
from Cython.Build import cythonize

with open("README.md", "r") as fh:
    long_description = fh.read()


setup(
    name="crispr-bean",
    version="0.2.3",
    python_requires=">=3.8.0",
    author="Jayoung Ryu",
    author_email="jayoung_ryu@g.harvard.edu",
    description="Base Editor screen analysis [Bayesian Estimation of variant effect] with guide Activity Normalization",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pinellolab/crispr-bean",
    packages=find_packages(),
    ext_modules=cythonize(["bean/mapping/CRISPResso2Align.pyx"]),
    include_dirs=np.get_include(),
    setup_requires=[
        "setuptools>=18.0",
        "cython",
    ],
    scripts=[
        "bin/bean-count",
        "bin/bean-count-samples",
        "bin/bean-qc",
        "bin/bean-filter",
        "bin/bean-run",  # TODO: prevent error when extra requirements are not met.
    ],
    install_requires=[
        "numpy",
        "perturb-tools>=0.2.0",
        "anndata>=0.8.0",
        "Bio>=1.5",
        "matplotlib",
        "pandas",
        "scipy",
        "seaborn",
        "tqdm",
    ],
    extras_require={"model": ["pyBigWig", "pyro-ppl<=1.8.1", "statsmodels", "torch<2"]},
    include_package_data=True,
    package_data={
        "": [
            "bean/annotate/ldlr_exons.fa",
            "notebooks/sample_quality_report.ipynb",
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
)
