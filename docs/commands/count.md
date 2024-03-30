# `bean count[-samples]`: Count (reporter) screen data  
`bean count-samples` (or `bean count` for a single sample) maps guide into guide counts, **allowing for base transition in spacer sequence**. When the matched reporter information is provided, it can count the **target site edits** and **alleles produced by each guide**. Mapping is efficiently done based on [CRISPResso2](https://github.com/pinellolab/CRISPResso2) modified for base-edit-aware mapping.



```python
bean count-samples \
  --input sample_list.csv   `# sample with lines 'R1_filepath,R2_filepath,sample_name\n'` \
  -b A                      `# base that is being edited (A/G)` \
  -f sgRNA_info_table.csv   `# sgRNA information` \
  -o .                      `# output directory` \
  -r                        `# read edit/allele information from reporter` \
  -t 12                     `# number of threads` \
  --name my_sorting_screen  `# name of this sample run` \
```
```python
bean count --R1 R1.fq --R2 R2.fq -b A -f sgRNA_info_table.csv -r
```
By default, `bean count[-samples]` assume R1 and R2 are trimmed off of the adapter sequence. You may need to adjust the command arguments according to your read structure. 

   <img src="../imgs/sequence_struct.png" alt="Read structuren" width="600"/>  

See full detail [below](#full-parameters).

# Input file format
See :ref:`input` for input file formats.

# Output file format
`count` or `count-samples` produces `.h5ad` and `.xlsx` file with guide and per-guide allele counts.  
* `.h5ad`: This output file follows annotated matrix format compatible with `AnnData` and is based on `Screen` object in [purturb_tools](https://github.com/pinellolab/perturb-tools). See [Data Structure](#data-structure) section for more information.  
* `.xlsx`: This output file contains `.guides`, `.samples`, `.X[_bcmatch,_edits]`. (`allele_tables` are often too large to write into an Excel!)  
