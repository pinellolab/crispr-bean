import numpy as np
from setuptools import setup, find_namespace_packages
from Cython.Build import cythonize

with open("README.md", "r") as fh:
    long_description = fh.read()


setup(
    name="crispr-bean",
    version="1.2.1",
    python_requires=">=3.8.0",
    author="Jayoung Ryu",
    author_email="jayoung_ryu@g.harvard.edu",
    description="Base Editor screen analysis [Bayesian Estimation of variant effect] with guide Activity Normalization",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pinellolab/crispr-bean",
    packages=find_namespace_packages(),
    ext_modules=cythonize(["bean/mapping/CRISPResso2Align.pyx"]),
    include_dirs=np.get_include(),
    setup_requires=[
        "setuptools>=18.0",
        "Cython",
    ],
    scripts=["bin/bean"],
    install_requires=[
        "numpy",
        "pandas",
        "scipy",
        "perturb-tools>=0.3.3",
        "matplotlib",
        "seaborn>=0.13.0",
        "tqdm",
        "bio",
        "liftover",
        "openpyxl>=3.0.10",
        "papermill>=2.4.0",
        "pyBigWig>=0.3.18",
        "pyro-ppl>=1.8.5",
        "scikit-learn",
        "statsmodels>=0.12.1",
        "ipykernel",
        "pytest-order",
        "nbconvert",
        "logomaker",
    ],
    extras_require={"model": ["pyBigWig", "pyro-ppl", "statsmodels", "torch"]},
    include_package_data=True,
    package_data={
        "bean.annotate": ["ldlr_exons.fa"],
        "bean.notebooks": [
            "sample_quality_report.ipynb",
            "profile_editing_preference.ipynb",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
)
