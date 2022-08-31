
#!/bin/bash
## Run all dimensions with a set of criteria


gene_prev_filters=(0.1 0.5 0.8)
direct_filters=("True" "False")
metrics=("manhattan" "euclidean" "braycurtis")
layers=("True" "False")
i=1
for gene_prev_filter in ${gene_prev_filters[@]}; do
	for direct_filter in ${direct_filters[@]}; do
		for metric in ${metrics[@]}; do
			for layer2 in ${layers[@]}; do
				epochs=5000
				SGE_Batch -c "/nfs3/PHARM/David_Lab/christine/software/Python-3.6.2/bin/python3 between_omic_encoding.py  --metric $metric --layer2 $layer2 --direct_filter $direct_filter --epochs $epochs --gene_prev_filter $gene_prev_filter --load_file '' --embedding_size 10" -r mtg_metabol_10_$i -P 20

				SGE_Batch -c "/nfs3/PHARM/David_Lab/christine/software/Python-3.6.2/bin/python3 between_omic_encoding.py  --metric $metric --layer2 $layer2 --direct_filter $direct_filter --epochs $epochs --gene_prev_filter $gene_prev_filter --load_file '' --embedding_size 25" -r mtg_metabol_25_$i -P 20

				SGE_Batch -c "/nfs3/PHARM/David_Lab/christine/software/Python-3.6.2/bin/python3 between_omic_encoding.py  --metric $metric --layer2 $layer2 --direct_filter $direct_filter --epochs $epochs --gene_prev_filter $gene_prev_filter --load_file '' --embedding_size 50" -r mtg_metabol_50_$i -P 20

				SGE_Batch -c "/nfs3/PHARM/David_Lab/christine/software/Python-3.6.2/bin/python3 between_omic_encoding.py  --metric $metric --layer2 $layer2 --direct_filter $direct_filter --epochs $epochs --gene_prev_filter $gene_prev_filter --load_file '' --embedding_size 100" -r mtg_metabol_100_$i -P 20

				SGE_Batch -c "/nfs3/PHARM/David_Lab/christine/software/Python-3.6.2/bin/python3 between_omic_encoding.py  --metric $metric --layer2 $layer2 --direct_filter $direct_filter --epochs $epochs --gene_prev_filter $gene_prev_filter --load_file '' --embedding_size 200" -r mtg_metabol_200_$i -P 20

				SGE_Batch -c "/nfs3/PHARM/David_Lab/christine/software/Python-3.6.2/bin/python3 between_omic_encoding.py  --metric $metric --layer2 $layer2 --direct_filter $direct_filter --epochs $epochs --gene_prev_filter $gene_prev_filter --load_file '' --embedding_size 500" -r mtg_metabol_500_$i -P 20

				i=$((i+1))
			done
		done
	done
done
#gene_prev_filter=0.8 #include gene if it is present in at least 1-gene_prev_filter percent of samples
#direct_filter="False"
#metric="braycurtis"
#layer2="False"
#i=36

