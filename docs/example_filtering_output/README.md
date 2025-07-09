# `bean filter` output files
These are the example output of [`bean filter`](https://pinellolab.github.io/crispr-bean/filter.html).
* `...filtered_allele_stats.pdf`: Plots the number of alleles per gRNA and guides per targeted variant (element) at each allele filtering stage.
* `...filter_log.txt`: Stores number of alleles and key at each allele filtering stage.
* `...alleleFiltered.h5ad`: AnnData file with filtered alleles saved in `.uns`. 
    ```
    # This results in the following ReporterScreen object
    >>> bdata
    Genome Editing Screen comprised of n_guides x n_conditions = 7500 x 25
        guides:    'Region', 'pos', 'strand', 'sequence', 'reporter', 'barcode', '5-nt PAM', 'targetPos', 'pos_seq', 'guide_len', 'start_pos', 'target_start', 'masked_sequence', 'masked_barcode', 'target_group', 'edit_rate', 'chrom'
        samples:   'replicate', 'condition', 'gini_X', 'median_corr_X', 'median_lfc_corr.top_bot', 'median_editing_rate', 'mask', 'lower_quantile', 'upper_quantile'
        samples_m:
        samples_p:
        layers:    'X_RPM', 'X_bcmatch', 'edit_rate', 'edits', 'lognorm_counts'
        uns:       'allele_counts', 'allele_counts_0_19', 'allele_counts_0_19_noindels', 'allele_counts_0_19_noindels_A>G', 'allele_counts_0_19_noindels_A>G_translated', 'allele_counts_0_19_noindels_A>G_translated_prop0.1_0.3', 'edit_counts', 'lfc', 'lfc_corr', 'metadata', 'repguide_mask', 'target_base_change', 'target_base_changes', 'tiling'
    # Exporting table to print as example rows
    >>> pd.concat([
        bdata.uns['allele_counts_0_19_noindels_A>G_translated_prop0.1_0.3'].iloc[2000:2005], 
        bdata.uns['allele_counts_0_19_noindels_A>G_translated_prop0.1_0.3'].iloc[-5:]
        ]
        ).to_csv("example_filtered_output.csv")
    ```
    `example_filtered_output.csv` would look as follows:  
    * `aa_allele`: Formatted as `CODING_EDIT|NONCODING_EDIT`.
        * `CODING_EDIT` is formatted as `GENE_ID:AA_POS:AA_REF>AA_ALT`.
        * `NONCODING_EDIT` is formatted as `[CHROM:]GENOME_POS:EDIT_POSITION_IN_REPORTER:STRAND:REF>ALT`.

    
    |index |guide                      |aa_allele                                                                               |rep5_top|rep5_high|rep5_bulk|rep5_low|rep5_bot|rep6_top|rep6_high|rep6_bulk|rep6_low|rep6_bot|rep7_top|rep7_high|rep7_bulk|rep7_low|rep7_bot|rep8_top|rep8_high|rep8_bulk|rep8_low|rep8_bot|rep9_top|rep9_high|rep9_bulk|rep9_low|rep9_bot|
    |------|---------------------------|----------------------------------------------------------------------------------------|--------|---------|---------|--------|--------|--------|---------|---------|--------|--------|--------|---------|---------|--------|--------|--------|---------|---------|--------|--------|--------|---------|---------|--------|--------|
    |2000  |13_3374_neg                |LDLR_EXONS:630:S>S&#124;                                                                     |117     |68       |84       |26      |38      |97      |129      |109      |114     |272     |34      |18       |27       |54      |21      |20      |46       |74       |65      |33      |7       |52       |24       |56      |0       |
    |2001  |13_3375_neg                |LDLR_EXONS:630:S>S&#124;                                                                     |266     |236      |133      |83      |253     |234     |142      |201      |164     |63      |107     |128      |99       |74      |100     |121     |54       |93       |154     |76      |3       |133      |15       |82      |0       |
    |2002  |13_3375_pos                |LDLR_EXONS:630:S>G&#124;                                                                     |61      |94       |102      |2       |45      |28      |61       |1        |11      |78      |9       |4        |19       |11      |23      |12      |10       |6        |32      |17      |4       |16       |1        |2       |0       |
    |2003  |13_3376_neg                |LDLR_EXONS:630:S>S&#124;                                                                     |146     |80       |3        |28      |7       |195     |107      |42       |121     |43      |31      |19       |14       |18      |28      |19      |20       |26       |36      |70      |1       |29       |9        |24      |0       |
    |2004  |13_3376_pos                |LDLR_EXONS:630:S>G,LDLR_EXONS:632:N>D&#124;                                                  |125     |24       |47       |97      |9       |24      |1        |52       |44      |1       |9       |17       |24       |47      |17      |27      |0        |11       |15      |11      |0       |23       |2        |14      |0       |
    |11857 |Intron 1 DNaseHS 2_5704_pos|&#124;11203065:6:+:A>G,11203079:20:+:A>G                                                     |21      |48       |12       |66      |0       |16      |54       |62       |23      |16      |20      |8        |10       |11      |26      |36      |11       |7        |22      |23      |0       |43       |4        |10      |8       |
    |11858 |Intron 1 DNaseHS 2_5707_pos|&#124;11203077:15:+:A>G,11203079:17:+:A>G                                                    |0       |9        |10       |32      |39      |25      |25       |53       |9       |21      |13      |25       |16       |12      |39      |13      |28       |15       |19      |47      |0       |40       |5        |28      |16      |
    |11859 |Intron 1 DNaseHS 2_5707_pos|&#124;11203079:17:+:A>G                                                                      |199     |25       |81       |112     |71      |141     |80       |103      |100     |81      |66      |15       |14       |21      |46      |46      |58       |50       |51      |39      |8       |98       |7        |54      |0       |
    |11860 |Intron 1 DNaseHS 2_5710_pos|&#124;11203079:14:+:A>G                                                                      |68      |0        |24       |0       |7       |3       |1        |30       |40      |0       |9       |5        |0        |4       |13      |15      |4        |9        |9       |24      |5       |8        |5        |14      |9       |
    |11861 |Intron 1 DNaseHS 2_5713_pos|&#124;11203075:7:+:A>G,11203076:8:+:A>G,11203079:11:+:A>G,11203083:15:+:A>G,11203088:20:+:A>G|0       |0        |28       |4       |0       |0       |0        |0        |0       |0       |0       |1        |8        |13      |1       |16      |6        |7        |0       |5       |4       |2        |3        |0       |7       |
