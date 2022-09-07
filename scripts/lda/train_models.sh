#input_file="ps_16s_dds.rds"
#omic="16s"
#filter=0.05

#input_file="ps_mtg_rle_nooutliers_adjcounts.rds"
#omic="mtg"
#filter=0.1

#input_file="ps_mtt_rle_nooutliers_adjcounts.rds"
#omic="mtt"
#filter=0.1

#input_file="ps_mbx_rle_nooutliers_adjcounts_fixedmapping.rds"
#omic="mbx"
#filter=0.1

#input_file="ps_plasma_cuda_norm.rds"
#omic="plasma"
#filter=0.0


#input_file="ps_liver_cuda_norm.rds"
#omic="liver"
#filter=0.0


#input_file="ps_brain_cuda_norm.rds"
#omic="brain"
#filter=0.0


input_file="ps_stool_cuda_norm.rds"
omic="stool"
filter=0.0

iter=50000
for i in {1..20..1}
do
	rm -rf log_lda_"$omic"_"$i"topic
	SGE_Batch -c "/nfs3/PHARM/David_Lab/christine/software/R-3.6.0/bin/Rscript latent_dirichlet_allocation.R $input_file model_"$omic"_"$i"topics_"$iter"iter_"$filter"filter.rds $filter $i $iter" -r log_lda_"$omic"_"$i"topic -q mus -P 2

done
#i=3
#SGE_Batch -c "/nfs3/PHARM/David_Lab/christine/software/R-3.6.0/bin/Rscript latent_dirichlet_allocation.R $input_file model_"$omic"_"$i"topics_"$iter"iter_"$filter"filter.rds $filter $i $iter" -r log_lda_"$omic"_"$i"topic -q mus -P 2
