import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
<<<<<<< HEAD
     name='berets',  
     version='0.0.6',
=======
     name='beret',  
     version='0.5.1',
>>>>>>> 86778063f13de9baa0eff426d453afaf7fcbdd1e
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
        'perturb-tools>=0.0.16',
      ],
     classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: MIT License",
         "Operating System :: OS Independent",
     ],
 )

