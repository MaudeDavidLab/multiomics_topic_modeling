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
from sklearn.ensemble import RandomForestRegressor

from collections import Counter
from mord import LogisticAT
import math
from matplotlib.offsetbox import AnchoredText
from os import listdir
from copy import copy
import seaborn as sns

def transformOrdinal_fruit(y):
    y[[math.isnan(i) for i in y]] = np.nanmedian(y)
    target = np.array([int(i) for i in y])
    new_target = []
    for i in target:
        if i == 4:
            new_target = new_target + [2]
        elif i == 5: 
            new_target = new_target + [3]
        elif i == 6:
            new_target = new_target + [4]
        else:
            new_target = new_target + [i]
    return(np.array(new_target))

def plotRegression(y, preds):
    m, b = np.polyfit(y, preds, 1)
    r, p = pearsonr(y, preds)
    #ax.scatter(y, preds)
    #ax.plot(y, m*y + b) 
    #text_box = AnchoredText("R=" + str(round(r, 3)) + ", P=" + str(round(p, 3)), 
    #                        frameon = True, loc = 4, pad = 0.5, prop=dict(fontsize=15))
    #ax.add_artist(text_box)
    #ax.text(x = 3, y = 1, horizontalalignment = "center", verticalalignment = "bottom", s = )
    error = np.abs(1 - m)
    #error = np.sum([np.square(i - j) for i,j in zip(y, model.predict(data))])
    return(r)

def getConfusionMatrix(model, data, y, ax):
    preds =  model.predict(data)
    #sns.heatmap(confusion_matrix(y, preds), annot = True, ax = ax)
    m, b = np.polyfit(y, preds, 1)
    ax.scatter(y, preds)
    ax.plot(y, m*y + b) 
    p = pearsonr(y, preds)[1]
    

def crossValidation(x, y, mapping, model, state = 0, plot = False):
    kf = KFold(n_splits = 3, shuffle = True, random_state = 0) #
    family_ids_unique = np.unique(mapping.familyID.values)
    rs = []
    ps = []
    for train_index, test_index in kf.split(family_ids_unique):
        train_bool = [i in family_ids_unique[train_index] for i in mapping.familyID.values]
        test_bool = [i in family_ids_unique[test_index] for i in mapping.familyID.values]
        model.fit(x.loc[train_bool, :], y[train_bool])
        preds = model.predict(x.loc[test_bool, :])
        r, p = pearsonr(y[test_bool], preds)
        rs = rs + [r]
        ps = ps + [p]
    return(rs, ps)

def getModel(model_type, param):
    if model_type == "elastic_net":
        model = ElasticNet(l1_ratio = l1_ratio,  alpha=param, random_state = 0, max_iter = 10000)
    if model_type == "ordinal_log":
        model = LogisticAT(alpha = param)
    if model_type == "logistic":
        model = LogisticRegression(C = param)
    if model_type == 'rf':
        model = RandomForestRegressor(max_depth = param, random_state = 0)
    return(model)

def getParams(model_type):
    if model_type == "elastic_net":
        params = [10, 100, 1000]
    if model_type == "ordinal_log":
        params = [100000, 1000000]
    if model_type == "logistic":
        params = [0, 100, 1000, 10000, 100000, 1000000]
    if model_type == "rf":
        #params = [2, 3, 4, 5, 10, 20, 30, 40, 50]
        params = [2, 4, 6, 8, 10, 20, 30, 40, 50, 60, 70, 100, 120]
    return(params)

def getBestParam(x, y, mapping, model_type, params):
    rs_crossVal = []
    params = [2,3]
    for param in params:
        model = getModel(model_type, param)
        rs, ps = crossValidation(x, y, mapping, model)
        rs_crossVal = rs_crossVal + [np.mean(rs)]
    best_performing_index = np.argmax(rs_crossVal)
    best_param = params[best_performing_index]
    return(best_param)

def cleanMLData(x, mapping, criteria, transform = None):
    y = mapping[criteria].values
    keep = [i == i for i in y]
    y = y[keep]
    x = x.loc[keep, :]
    mapping = mapping.loc[keep, :]
    if not transform == None:
        y = transform(y)
    y = np.array([int(i) for i in y])
    return(x, y, mapping)

def getPredictions(x_train, y_train, x_test, y_test, model_type, best_param):
    #x_train, y_train, mapping_train = cleanMLData(x.iloc[train_bool, :], mapping.iloc[train_bool, :], criteria, transform)
    #x_test, y_test, mapping_test = cleanMLData(x.iloc[test_bool, :], mapping.iloc[test_bool, :], criteria, transform)
    #params = getParams(model_type)
    #best_param = getBestParam(x_train, y_train, mapping_train, model_type, params)
    model = getModel(model_type, best_param)
    model.fit(x_train, y_train)
    preds = model.predict(x_test)
    preds_train = model.predict(x_train)
    return(y_train, preds_train, y_test, preds)

def plotBars(rs_raw, rs_encode, omic_names, criteria):
    fig = plt.figure()
    df = pd.DataFrame({'raw':rs_raw, 'encode':rs_encode})
    df['omic'] = omic_names
    df.set_index('omic')
    df = df.melt('omic')
    sns.barplot(x = "omic", y = "value", hue = "variable", data = df).set_title(criteria)
    plt.ylim(-0.5, 0.7)
    
def setupFigure(row_names, col_names = ["Training",  "Testing"]):
    fig, axes = plt.subplots(nrows=5, ncols=2, figsize=(12, 9))
    for ax, col in zip(axes[0], col_names):
        ax.set_title(col, fontsize = 20)

    for ax, row in zip(axes[:,0], row_names):
        ax.set_ylabel(row, rotation=90, size='large', fontsize = 20)
    return(axes)   

def trainAndTestModel(data_storage, mapping_storage, train_bool, test_bool, criteria, model_type, transform):
    omic_names = data_storage.getOmicNames()
    #print(omic_names)
    #axes = setupFigure(row_names = omic_names)
    rs = []
    i = 0
    for x, mapping in zip(data_storage.getList(), mapping_storage.getList()):
        #print(i)
        print(omic_names[i])
        x_train, y_train, mapping_train = cleanMLData(x.iloc[train_bool, :], mapping.iloc[train_bool, :], criteria, transform)
        x_test, y_test, mapping_test = cleanMLData(x.iloc[test_bool, :], mapping.iloc[test_bool, :], criteria, transform)
        params = getParams(model_type)
        best_param = getBestParam(x_train, y_train, mapping_train, model_type, params) 
        rs_tmp = []
        if model_type == 'rf':
            for j in np.arange(7):
                y_train, preds_train, y_test, preds = getPredictions(x_train, y_train, x_test, y_test, model_type, best_param)
                r, p = pearsonr(y_test, preds)
                rs_tmp = rs_tmp + [r]
            r_test = np.mean(rs_tmp)
        else:
            #r_train = plotRegression(y_train, preds_train, axes[i, 0])
            #r_test = plotRegression(y_test, preds, axes[i, 1])
            r_train = plotRegression(y_train, preds_train)
            r_test = plotRegression(y_test, preds)
        rs = rs + [r_test]
        i = i + 1
    return(rs)

def runCompareOmicsPipeline(data_storage, encoded_storage, mapping_storage, train_bool, test_bool, criteria, model_type, transform):
    rs_raw = trainAndTestModel(data_storage, mapping_storage, train_bool, test_bool, criteria, model_type, transform)
    rs_encode = trainAndTestModel(encoded_storage, mapping_storage, train_bool, test_bool, criteria, model_type, transform)
    plotBars(rs_raw, rs_encode, data_storage.getOmicNames(), criteria)  
    return(rs_raw, rs_encode)

    
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
    