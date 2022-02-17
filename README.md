# crisprep
Tool for analyzing CRISPR data with reporter edit counts.

### Read count
```
crisprep-count-samples  \
  --input sample_list.csv   \ # sample with lines 'R1_filepath,R2_filepath,sample_name\n'  
  -b A  \ # base that is being edited (A/G)
  -f ../../../gRNA_info/LDLvar_gRNA_crisprep.csv  \ # sgRNA information 
  -o .  \ # output directory    
  -a  \ # read allele information  
  -r  \ # read reporter edit information  
  -t 12  \ # number of threads  
  --name 120821_LDLvar_fullsort  \ # name of this sample run  
```

## Using as python module
```
import crisprep as cp
cdata = cp.read_h5ad("crisprepCounts.h5ad")
```

### Combine reps
```
cdata_combined = cdata_jul + cdata_oct
```

### Annotate alleles
```
from crisprep.annotate.translate_allele import *  
# TBD
```
