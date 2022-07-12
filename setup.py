import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
     name='berets',  
     version='0.0.2',
     author="Jayoung Ryu",
     author_email="jayoung_ryu@g.harvard.edu",
     description="Tools for analyzing CRISPR data with REPorter edits",
     long_description=long_description,
     long_description_content_type="text/markdown",
     url="https://github.com/pinellolab/beret",
     packages=setuptools.find_packages(),
     scripts=["bin/beret-count",
     "bin/beret-count-samples"],
     install_requires=[
        'numpy',
        'perturb-tools',
      ],
     classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: MIT License",
         "Operating System :: OS Independent",
     ],
 )

