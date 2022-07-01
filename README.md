# beret
Base Editing with REporter analysis Toolkit.

## Installation
Clone and do 
```
pip install -e .
```

## Count reporter screen data
```
beret-count-samples  \
  --input sample_list.csv   \ # sample with lines 'R1_filepath,R2_filepath,sample_name\n'  
  -b A  \ # base that is being edited (A/G)
  -f ../../../gRNA_info/LDLvar_gRNA_beret.csv  \ # sgRNA information 
  -n sample  \ # number of sample  
  -o .  \ # output directory    
  -a  \ # read allele information  
  -r  \ # read reporter edit information
  -m  \ # read matched guide and reporter edit information  
  -t 12  \ # number of threads  
  --name 120821_LDLvar_fullsort  \ # name of this sample run  
```

This produces `.h5ad` and `.xlsx` file.

## Analysis
```
adata.log_norm
adata.log_fold_change_aggregate("bot", "top")
``` 

## Using as python module
```
import beret as br
cdata = br.read_h5ad("beret_counts_sample.h5ad")
```

See the [tutorial](beret_test.rst) for more detail.
