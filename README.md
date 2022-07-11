# <img src="beret.svg" alt="beret" width="300"/>
**B**ase **E**diting with **Re**porter analysis **T**oolkit.  
This is an analysis toolkit for the pooled CRISPR reporter or sensor data. The reporter technique transfects cells with plasmid with not only sgRNA but with the **target sequence surrogate** which we call **reporter** or **sensor**.  
  
<img src="anbe.svg" alt="anbe" width="500"/>

## Installation 
```
git clone https://github.com/pinellolab/beret.git
cd beret/
pip install -e .
```

## Count reporter screen data  
Aligns guide with matched reporter allele counts in multiple samples.  

<img src="reporter_screen.svg" alt="reporter screen" width="700"/>  

```python
beret-count-samples         \
  --input sample_list.csv   \ # sample with lines 'R1_filepath,R2_filepath,sample_name\n'  
  -b A                      \ # base that is being edited (A/G)
  -f gRNA_library.csv       \ # sgRNA information 
  -o .                      \ # output directory    
  -a                        \ # read allele information  
  -r                        \ # read reporter edit information
  -m                        \ # read matched guide and reporter edit information  
  -t 12                     \ # number of threads  
  --name LDLvar_fullsort    \ # name of this sample run  
```

This produces `.h5ad` and `.xlsx` file with guide and per-guide allele counts.  
`.h5ad` file follows annotated matrix format compatible with `AnnData` and is based on `Screen` object in [purturb_tools](https://github.com/pinellolab/perturb-tools) and contains the per-guide allele counts.    
<img src="screendata.svg" alt="screendata" width="700"/>

## Using as python module
```
import beret as br
cdata = br.read_h5ad("beret_counts_sample.h5ad")
```

See the [tutorial](beret_test.rst) for more detail.
