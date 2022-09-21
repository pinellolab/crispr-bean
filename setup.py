from Cython.Build import cythonize
import numpy as np
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
     name='berets',  
     version='0.1.1',
     author="Jayoung Ryu",
     author_email="jayoung_ryu@g.harvard.edu",
     description="Tools for analyzing CRISPR data with REPorter edits",
     long_description=long_description,
     long_description_content_type="text/markdown",
     url="https://github.com/pinellolab/beret",
     packages=setuptools.find_packages(),
     ext_modules=cythonize(["beret/preprocessing/CRISPResso2Align.pyx"]),
     include_dirs=np.get_include(),
     scripts=["bin/beret-count",
     "bin/beret-count-samples"],
     install_requires=[
        'numpy',
        'perturb-tools>=0.0.16',
      ],
      include_package_data=True,
    package_data={'': ['beret/annotate/ldlr_exons.fa']},
     classifiers=[
         "Programming Language :: Python :: 3",
         "Operating System :: OS Independent",
     ],
 )

