# multiomics_topic_modeling workflow
# 1. Data
Processed data files can be found at : https://files.cqls.oregonstate.edu/David_Lab/M3/
Raw data files available upon request.

# 2. Process/Normalization/Filtering
## A. Files are processed as described in Tataru et. al (2022). Additionally, children still breast feeding and children under 2 years old were removed to create files:
- 16s/ps_not_norm_comp_pass_min_postDD_age_filtered.rds
- mtg/PS_kegg_fun_MTG_age_filtered.rds
- mtt/PS_kegg_fun_MTT.rds
- mbx/metabolomics_data_orgi_scale.csv

## B. Metadata was cleaned and standardized using the script ```metadata/clean_mapping_ml.R``` 


## C. Normalizations were tested and selected using the script ```processing/normalization_choices.R``` to create the files:
- 16s/ps_16s_dds_taxannotation.rds
- mtg/ps_mtg_rle_nooutliers_adjcounts.rds
- mtt/ps_mtt_rle_nooutliers_adjcounts.rds
- mbx/ps_mbx_rle_nooutliers_adjcounts_fixedmapping.rds

### D. Prevalence filtering was performed using the "Filter" chunk of script ```lda/latent_dirichlet_allocation_3```

# 3. Topic modeling
# 4. Clustering topics
# 5. Interpreting topics
# 6. Clustering samples
# 7. Differential analysis


