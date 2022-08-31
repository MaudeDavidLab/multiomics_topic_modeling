import os
import itertools
import numpy as np

def getWeightDir(omic, metric):
    #Name file
    
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


omic = "s"
metrics = ['braycurtis', 'horn', 'manhattan']

for metric in metrics:
	weight_dir = "weights/autoencoder/" + omic + "/" + metric + "/"
	for dim in ['10', '25', '50', '100', '200', '500']:
		keep_file_name = getFileName(weight_dir, dim)
		weight_dir_tmp = weight_dir + dim + "dim/"
		files = os.listdir(weight_dir_tmp)
		print(keep_file_name)
		for file in files:
			if file != keep_file_name:
				os.remove(weight_dir_tmp + file)
