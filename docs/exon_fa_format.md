# Input .fa file format for `bean-filter`
You can provide custom FASTA file with exon sequence entries. Currently only supports positive strand genes.

* Exon FASTA files can be downloaded from UCSC Genomic sequences / Table Browser: [see the instruction video](https://www.youtube.com/watch?v=T4E0Ez5Vjz8)
* You can manually format as: 
    * Header line has ` range=chrom:start-end ` and `strand=+/-` tag that is parsed.
    * fasta entry has the sequence of exons, where the first (includes 5'-UTR) and last (includes 3'-UTR) exon sequence has lower-case sequence denoting noncoding sequences.
* See the example .fa [here](../tests/data/ldlr_exons.fa).