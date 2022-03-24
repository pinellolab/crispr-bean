# crisprep
Tool for analyzing CRISPR data with reporter edit counts.

### Installation
Clone and do 
```
pip install -e .
```

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
This is only implemented for LDLR. 
```
from crisprep.framework._supporting_fn import get_aa_alleles
allele = cp.Allele.from_str('11222248:5:-:A>G,11222252:1:-:A>G')  # abs_pos:rel_pos:strand:base_edit
get_aa_alleles(allele, include_synonymous = False)
$ ['375:Y>H']
```
