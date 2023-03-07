import numpy as np
import setuptools
from Cython.Build import cythonize

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="berets",
    version="0.1.3",
    python_requires=">=3.8.0",
    author="Jayoung Ryu",
    author_email="jayoung_ryu@g.harvard.edu",
    description="Base Editor with or without REporter data analysis Toolkit",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pinellolab/beret",
    packages=setuptools.find_packages(),
    ext_modules=cythonize(["beret/mapping/CRISPResso2Align.pyx"]),
    include_dirs=np.get_include(),
    scripts=["bin/beret-count", "bin/beret-count-samples", "bin/beret-qc"],
    install_requires=[
        "numpy",
        "perturb-tools>=0.0.16",
    ],
    include_package_data=True,
    package_data={
        "": [
            "beret/annotate/ldlr_exons.fa",
            "beret/notebooks/sample_quality_report.ipynb",
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
)
