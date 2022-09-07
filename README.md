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

# 3. Train topic models using scripts ```lda/train_models.R``` to train the models and ```lda/train_models.sh``` to submit training jobs to a server.
### A. Evaluate the models by creating simulated data from each model, then comparing sample quantile, feature quantile, and marginal pairwise distribution between simulated and true data using script ```lda/model_evaluation.R```

# 4. Clustering topics into cross-omic topics using the "Cluster topics" chunk of script ```lda/latent_dirichlet_allocation_3```

# 5. Interpret topics based on high feature weights using the "Feature Weights MBX/MTG/MTT/16s" chunks of script ```lda/latent_dirichlet_allocation_3```.
These chunks create data for Figure 2 in the manuscript, which is the post-processed manually for a cleaner visualization.

# 6. Cluster samples using the "Full heatmap" chunck of script ```lda/latent_dirichlet_allocation_3```

# 7. Differential analysis using the "MBX/MTG/MTT/16s differences between groups" chunks of script ```lda/latent_dirichlet_allocation_3```
These chunks produce Figure 3 in the manuscript

### script ```lda/latent_dirichlet_allocation_3``` also produced all supplementary figures in the manuscript, with the exception of the model_evaluation, which is produced by the model_evalutation.R script.

* Please direct all question about reproducing code to ctataru15@gmail.com
