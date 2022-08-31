import os
import itertools
import numpy as np

def getWeightDir(omic, metric, layer2, gene_prev_filter):
    #Name file    
    weight_dir = "weights/" + omic + "/"
    weight_dir = weight_dir + "prev_filt" + str(gene_prev_filter) + "/"
    if layer2 == "True":
        weight_dir = weight_dir + "layer2/"
    else:
        weight_dir = weight_dir + "layer1/"
    weight_dir = weight_dir + metric + "/"

    return(weight_dir)


def getFileName(weight_dir, dim):
    best_val_scores = []
    weight_dir_tmp = weight_dir + dim + "dim/"
    weight_files = os.listdir(weight_dir_tmp)

    #pick best val score in that folder
    val_scores = [float(i.split('val')[1][0:5]) for i in weight_files]
    best_val_score = np.min(val_scores)
    file_name = weight_files[np.argmin(val_scores)]
    return(file_name)


gene_prev_filters = [0.1, 0.5, 1.0]
direct_filters = ["True", "False"]
metrics = ['braycurtis', 'horn', 'manhattan', 'euclidean', ]
layers = ["True", "False"]
omic = "mtg_16s"
for gene_prev_filter in gene_prev_filters: 
	for metric in metrics:
		for layer in layers:
			weight_dir = getWeightDir(omic = omic, metric = metric, layer2 = layer, gene_prev_filter = gene_prev_filter)
			for dim in ['10', '25', '50', '100', '200', '500']:
				keep_file_name = getFileName(weight_dir, dim)
				weight_dir_tmp = weight_dir + dim + "dim/"
				files = os.listdir(weight_dir_tmp)
				print(keep_file_name)
				for file in files:
					if file != keep_file_name:
						os.remove(weight_dir_tmp + file)