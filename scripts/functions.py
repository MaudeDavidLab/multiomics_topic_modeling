## dataset
#import tensorflow.keras as keras
## for Model definition/training
from tensorflow.keras.models import Model, load_model
from tensorflow.keras.layers import Input, Flatten, Dense, concatenate,  Dropout
from tensorflow.keras.optimizers import Adam
import tensorflow.keras.backend as K


from tensorflow.keras.utils import plot_model
from tensorflow.keras.callbacks import ModelCheckpoint

## required for semi-hard triplet loss:
from tensorflow.python.ops import array_ops
from tensorflow.python.ops import math_ops
from tensorflow.python.framework import dtypes
import tensorflow as tf

## for visualizing 
import matplotlib.pyplot as plt, numpy as np
from sklearn.decomposition import PCA
import pandas as pd
import seaborn as sns

import random
from sklearn.model_selection import train_test_split
from sklearn.utils import resample
from numpy.random import RandomState
from itertools import chain
from scipy.stats import spearmanr
from scipy.stats import pearsonr
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve
from sklearn.metrics import confusion_matrix
from sklearn.svm import SVC
from sklearn.linear_model import Ridge
from sklearn.linear_model import Lasso
from sklearn.model_selection import KFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import ElasticNet
from sklearn.linear_model import LogisticRegression

from collections import Counter
from mord import LogisticAT
import math
from matplotlib.offsetbox import AnchoredText
from os import listdir
import sys
from copy import copy
sys.path.append('functions.py')

        
def splitDatasets(data, mapping, state = 0, return_bool = False):
    family_ids = mapping.familyID.values
    train, test = train_test_split(list(set(family_ids)), test_size = .2, random_state = state)
    #test, val = train_test_split(test, test_size = .5, random_state = state)
    
    train_bool = [i in train for i in family_ids]
    #val_bool = [i in val for i in family_ids]
    test_bool = [i in test for i in family_ids]
    
    x_train = data.loc[train_bool, :]
    #x_val = data.loc[val_bool, :]
    x_test = data.loc[test_bool, :]

    #print((x_train.shape[0] + x_val.shape[0] + x_test.shape[0]) == data.shape[0])
    if return_bool == False:
        return(x_train, x_test)
    else:
        return(train_bool, test_bool)

def get16sData(filt):
    s_tab = pd.read_csv("../data/16s/asv_table_nooutliers.txt", sep = "\t")
    s_map = pd.read_csv("../data/16s/mapping_nooutliers.txt", sep = "\t")
    s_tab = s_tab.T
    s_map.head()

    if filt:
        keep = (s_tab == 0).astype(int).sum(axis = 0) / s_tab.shape[0] < .99
        s_tab = s_tab.loc[:, keep]
    print(s_tab.shape)

    return(s_tab, s_map)

def getMtgData(filt):
    mtg_tab = pd.read_csv("../data/mtg/asv_table_nooutliers.txt", sep = "\t")
    mtg_map = pd.read_csv("../data/mtg/mapping_nooutliers.txt", sep = "\t")
    mtg_tab = mtg_tab.T
    mtg_map.head()
    
    if filt:
        keep = ((mtg_tab == 0).astype(int).sum(axis = 0) / mtg_tab.shape[0]) < .8
        mtg_tab = mtg_tab.loc[:, keep]

    return(mtg_tab, mtg_map)

def getMttData(filt):
    mtt_tab = pd.read_csv("../data/mtt/asv_table_nooutliers.txt", sep = "\t")
    mtt_map = pd.read_csv("../data/mtt/mapping_nooutliers.txt", sep = "\t")
    mtt_tab = mtt_tab.T
    mtt_map.head()
    
    if filt:    
        keep = ((mtt_tab == 0).astype(int).sum(axis = 0) / mtt_tab.shape[0]) < .8
        mtt_tab = mtt_tab.loc[:, keep]
    print(mtt_map.shape)
    
    return(mtt_tab, mtt_map)

def brayCurtisLoss(y_true, y_pred):
    return(K.abs(K.sum(y_true - y_pred)) / K.sum(y_true + y_pred))

def getMetabolData():
    tab = pd.read_csv("../data/metabol/asv_table_nooutliers.txt", sep = "\t", index_col = 0)
    mapping = pd.read_csv("../data/metabol/mapping_nooutliers.txt", sep = "\t")
    tab = tab.iloc[:, 1:]
    tab = tab.T
    return(tab, mapping)

def getEmbedded16s():
    s_tab = pd.read_csv("../data/16s/embedded_agp.txt", sep = "\t")
    s_map = pd.read_csv("../data/16s/mapping_nooutliers.txt", sep = "\t")
    
    return(s_tab, s_map)

def getModel(file_name, metric):
    model = load_model(file_name, compile = False)
    if metric == "euclidean":
        model.compile(loss = 'mean_squared_logarithmic_error')
    if metric == "braycurtis":
        model.compile(loss = brayCurtisLoss)
    if metric == "manhattan":
        model.compile(loss = 'mean_absolute_error')
    return(model)

def encode(data, model_file):
    input_shape = data.shape[1]
    model = load_model(model_file, compile = False)
    model.compile(loss = brayCurtisLoss)
    input_vec = Input(shape= input_shape)
    x = Flatten()(input_vec)
    encoded = model.get_layer('encoder')(x)
    encoder = Model(inputs = input_vec, outputs = encoded)
    encoded = pd.DataFrame(encoder.predict(data), index = data.index)
    return(encoded)



def getMapping(data, mapping):
    mapping = pd.DataFrame()
    for i in data.index.values:
        mapping = mapping.append(mapping.loc[mapping.host_name == i, :], ignore_index = True)
    mapping = mapping.drop_duplicates(subset='host_name', keep='first')
    mapping.index = mapping.host_name
    return(mapping.loc[data.index.values, :])



def filterForProkaryotes(mbx_tab):
    metabolites_not_in_db = pd.read_csv('../data/metabol/metabolites_not_in_database.csv', index_col = 0)
    metabolites_in_db_directly = mbx_tab.columns.values[[i not in metabolites_not_in_db.index.values for i in mbx_tab.columns.values]]
    keep1 = metabolites_not_in_db.IN_MICROBIAL_DATABASE != 'FALSE'
    keep3 = metabolites_not_in_db.IN_MICROBIAL_DATABASE != 'FALSE; not in pubchem'
    keep4 = metabolites_not_in_db.IN_MICROBIAL_DATABASE != 'FALSE, not in pubchem'
    keep2 = [i == i for i in metabolites_not_in_db.IN_MICROBIAL_DATABASE] #check for NaN
    keep = [i and j and k and l for i,j,k,l in zip(keep1.values, keep2, keep3.values, keep4.values)]
    metabolites_in_db_indirectly = metabolites_not_in_db.loc[keep, :].index.values
    metabolites_keep = np.concatenate((metabolites_in_db_directly, metabolites_in_db_indirectly))
    return(mbx_tab[metabolites_keep])


#combine omics
def combineOmics(data1, data2, data3, data4):
    samples = set(data1.index.values).intersection(set(data2.index.values))
    samples = list(set(samples).intersection(set(data3.index.values)))
    samples = list(set(samples).intersection(set(data4.index.values)))
    #data = pd.concat([data1.loc[samples, :], data2.loc[samples, :], data2.loc[samples, :]], axis = 1)
    return(samples)

def formatMapping(mapping, samples_keep):
    mapping = mapping.drop_duplicates(subset= "host_name", keep = 'first')
    mapping.index = mapping.host_name
    mapping = mapping.loc[samples_keep, :]
    return(mapping)

def expandTimepoints(data, mapping, should_encode = False):
    expanded_time = pd.DataFrame()
    for host_name in np.unique(mapping.host_name):
        df_small = data.loc[mapping.host_name == host_name, :]
        if should_encode:
            df_small = encode(df_small, autoencoder_16s)
        expanded_timepoints = df_small.to_numpy().reshape((1, 3*df_small.shape[1]))[0]
        expanded_time = expanded_time.append(pd.Series(expanded_timepoints), ignore_index = True)
    expanded_time.index = np.unique(mapping.host_name)
    return(expanded_time)

def getWeightDir(omic, metric, layer2, direct_filter, gene_prev_filter):
    #Name file    
    if direct_filter:
        weight_dir =  'weights/'+ omic + "/" + "direct_filter/"
    else:
        weight_dir =  'weights/'+ omic + "/" + "no_filter/"
    weight_dir = weight_dir + "prev_filt" + str(gene_prev_filter) + "/"
    if layer2:
        weight_dir = weight_dir + "layer2/"
    else:
        weight_dir = weight_dir + "layer1/"
    weight_dir = weight_dir + metric + "/"
    print(weight_dir)
    return(weight_dir)

def getFileName(weight_dir, dim):
    best_val_scores = []
    weight_dir_tmp = weight_dir + dim + "dim/"
    weight_files = listdir(weight_dir_tmp)

    #pick best val score in that folder
    val_scores = [float(i.split('val')[1][0:5]) for i in weight_files]
    best_val_score = np.min(val_scores)
    file_name = weight_files[np.argmin(val_scores)]
    return(weight_dir_tmp + file_name)


def matchVariableSet(data, file_to_match):
    data_to_match = pd.read_csv(file_to_match)
    return(data.loc[:, data_to_match.columns.values[1:]])

def encodeBetweenOmics(data, omic, file, dim):
    str_split = file.split("_")
    gene_prev_filter = float(str_split[0])
    direct_filter = str_split[1]
    metric = str_split[2]
    layer2 = str_split[2]
 
    file_to_match = "data/" + str(omic) + "/x_test_genefilt" + str(gene_prev_filter) + ".csv"
    data_match = matchVariableSet(data, file_to_match)
    model_file = getFileName(getWeightDir(omic, metric, layer2, direct_filter, gene_prev_filter), dim)
    encoded_data = encode(data_match, model_file)
    return(encoded_data)

class Storage:
    def __init__(self, omics):
        self.s = copy(omics[0])
        self.mtg = copy(omics[1])
        self.mtt = copy(omics[2])
        self.mbx = copy(omics[3])
        self.all_omics = copy(omics[4])
    def getList(self):
        return([self.s, self.mtg, self.mtt, self.mbx, self.all_omics])
    def getOmicNames(self):
        return("16S", "MTG", "MTT", "MBX", "ALL")
    def pruneSamples(self, keep_samples):
        self.s = self.s.loc[keep_samples, :]
        self.mtg = self.mtg.loc[keep_samples, :]
        self.mtt = self.mtt.loc[keep_samples, :]
        self.mbx = self.mbx.loc[keep_samples, :]
        self.all_omics = self.all_omics.loc[keep_samples, :]
    def applyFunction(self, function):
        self.s = pd.DataFrame(function(self.s))
        self.mtg = pd.DataFrame(function(self.mtg))
        self.mtt = pd.DataFrame(function(self.mtt))
        self.mbx = pd.DataFrame(function(self.mbx))
        self.all_omic = pd.DataFrame(function(self.all_omics))
    
