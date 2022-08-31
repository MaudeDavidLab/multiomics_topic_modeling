## dataset
#import tensorflow.keras as keras
## for Model definition/training
from tensorflow.keras.models import Model, load_model
from tensorflow.keras.layers import Input, Flatten, Dense, concatenate,  Dropout
from tensorflow.keras.optimizers import Adam
import tensorflow.keras.backend as K
from tensorflow.keras import initializers
import tensorflow.keras.losses as losses


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
from sklearn.svm import SVC
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from os import listdir
import os

import argparse

#parse parameters
parser = argparse.ArgumentParser()
parser.add_argument('--embedding_size', type = str, dest = 'embedding_size')
parser.add_argument('--metric', type = str, dest = 'metric')
parser.add_argument('--layer2', type = str, dest = 'layer2')
#parser.add_argument('--direct_filter', type = str, dest = 'direct_filter')
parser.add_argument('--epochs', type = int, dest = 'epochs')
parser.add_argument('--load_file', type = str, dest = 'load_file')
parser.add_argument('--gene_prev_filter', type = str, dest = "gene_prev_filter")
parser.add_argument('--omic', type = str, dest = "omic")

args = parser.parse_args()
print(args)

def create_base_network(input_shape, embedding_size, output_shape, layer2= False):
    """
    Base network to be shared (eq. to feature extraction).
    """
    input_vec = Input(shape= input_shape)
    x = Flatten()(input_vec)
    
    encoded = Dense(embedding_size, activation='relu', name = 'encoder', kernel_initializer=initializers.RandomNormal(stddev=0.01))(x)
    if layer2 == "True":
        second = Dense(embedding_size, activation = 'relu', name = 'second', kernel_initializer=initializers.RandomNormal(stddev=0.01))(encoded)
        decoded = Dense(output_shape, activation = 'relu', name = 'decoder', kernel_initializer=initializers.RandomNormal(stddev=0.01))(second)
    else: 
        decoded = Dense(output_shape, activation = 'relu', name = 'decoder', kernel_initializer=initializers.RandomNormal(stddev=0.01))(encoded)

    autoencoder = Model(inputs=input_vec, outputs= decoded)
    
    
    return autoencoder
def getCorrs(use, model, y, metric):
    pred = model.predict(np.array(use))
    #corrs = getCosineSim(pred)
    #corrs_base = getCosineSim(np.array(use))
    corrs = squareform(pdist(pred, metric = metric))
    corrs_base = squareform(pdist(np.array(y), metric = metric))
    return(corrs, corrs_base)

def plotCorrs(corrs, corrs_base):
    fig = plt.figure(figsize=(15, 4))
    axarr = fig.subplots(1,2)
    p = axarr[0].matshow(corrs_base)
    p2 = axarr[1].matshow(corrs)
    print(mantel(corrs_base, corrs))
    fig.colorbar(p)
    fig.colorbar(p2)
    #Print mantel test metric

def encode(data, model_file):
    input_shape = data.shape[1]
    model = load_model(model_file, compile = False)
    
    model.compile(loss = brayCurtisLoss)
    input_vec = Input(shape= input_shape)
    x = Flatten()(input_vec)
    encoded = model.get_layer('encoder')(x)
    encoder = Model(inputs = input_vec, outputs = encoded)
    encoded = pd.DataFrame(encoder.predict(data))
    return(encoded)


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



def getXData(gene_prev_filter):
    global omic
    x_train = pd.read_csv('data/' + omic + 'x_train.csv', index_col = 0)
    x_val = pd.read_csv('data/' + omic + 'x_val.csv', index_col = 0)
    x_test = pd.read_csv('data/' + omic + 'x_test.csv', index_col = 0)

    #filter gene data
    x_tmp = pd.concat([x_train, x_val, x_test], )
    x_tmp = x_tmp.loc[:, (x_tmp==0).sum() < (x_tmp.shape[0] * gene_prev_filter)]
    x_train = x_tmp.loc[x_train.index.values, :]
    x_val = x_tmp.loc[x_val.index.values, :]
    x_test = x_tmp.loc[x_test.index.values, :]

    x_train = x_train + 1
    x_val = x_val + 1
    x_test = x_test + 1
    print("Number of features MTG: " + str(x_train.shape[1]))
    return(x_train, x_val, x_test)

def getYData(direct_filter = "False"):
    #load data
    global omic
    y_train = pd.read_csv('data/' + omic + 'y_train.csv', index_col = 0)
    y_val = pd.read_csv('data/' + omic + 'y_val.csv', index_col = 0)
    y_test = pd.read_csv('data/' + omic + 'y_test.csv', index_col = 0)
    if direct_filter == "True":
        y_train = filterForProkaryotes(y_train)
        y_val = filterForProkaryotes(y_val)
        y_test = filterForProkaryotes(y_test)
    print("Number of features MBX: ", str(y_train.shape[1]))
    return(y_train, y_val, y_test)

def getWeightDir(omic, dim, metric, layer2, direct_filter, gene_prev_filter):
    #Name file    
    #if direct_filter=="True":
    #    weight_dir =  'weights/'+ omic + "/" + "direct_filter/"
    #else:
    #    weight_dir =  'weights/'+ omic + "/" + "no_filter/"
    weight_dir = "weights/" + omic + "/"
    weight_dir = weight_dir + "prev_filt" + str(gene_prev_filter) + "/"
    if layer2 == "True":
        weight_dir = weight_dir + "layer2/"
    else:
        weight_dir = weight_dir + "layer1/"
    weight_dir = weight_dir + metric + "/" + dim + "dim/"
    print(weight_dir)
    if not os.path.exists(weight_dir):
        os.makedirs(weight_dir)
    return(weight_dir)

def brayCurtisLoss(y_true, y_pred):
    return(K.abs(K.sum(y_true - y_pred)) / K.sum(y_true + y_pred))
def lamb(y):
    return(K.sum(K.square(y), axis = 1) / K.square(K.sum(y, axis = 1)))
def hornLoss(y_true, y_pred):
    print("Y shape" + str(y_true.shape))
    return(1 - 2*K.dot(y_true, K.transpose(y_pred)) / ((lamb(y_true) + lamb(y_pred)) *K.sum(y_true, axis = 1) * K.sum(y_pred, axis = 1)))

def lamb_morisita(y):
    loss = K.dot(y, (K.transpose(y)-1)) / ( K.sum(y)* (K.sum(y)- 1) )
    #print(loss)
    return(loss)
    
def morisitaLoss(y_true, y_pred):
    #y_true = tf.math.round(y_true)
    #y_pred = tf.math.round(y_pred + 1)
    print("Y shape" + str(y_true.shape))
    return(1 - K.dot(y_true, K.transpose(y_pred)) / ((lamb_morisita(y_true) + lamb_morisita(y_pred)) *K.sum(y_true, axis = 1) * K.sum(y_pred, axis = 1)))


omic = args.omic + "/"
embedding_size = args.embedding_size
metric = args.metric
layer2 = args.layer2
epochs = args.epochs
gene_prev_filter = float(args.gene_prev_filter)
x_train, x_val, x_test = getXData(gene_prev_filter = gene_prev_filter)
y_train, y_val, y_test = getYData(direct_filter = direct_filter)
weight_dir = getWeightDir(omic, embedding_size, metric, layer2, direct_filter, gene_prev_filter)


input_shape = x_train.shape[1]
output_shape = y_train.shape[1]
batch_size = x_train.shape[0]

period = 20
initial_epoch = 0
learning_rate = .0001

if args.load_file != "":
    autoencoder = load_model(weight_dir + args.load_file)
else:
    autoencoder = create_base_network(input_shape, int(embedding_size), output_shape, layer2)
print(autoencoder.summary())
if args.metric == "euclidean":
    autoencoder.compile(optimizer = Adam(learning_rate), loss = 'mean_squared_logarithmic_error')
if args.metric == "braycurtis":
    autoencoder.compile(optimizer = Adam(learning_rate), loss = brayCurtisLoss)
if args.metric == "manhattan":
    autoencoder.compile(optimizer = Adam(learning_rate), loss = 'mean_absolute_error')
if args.metric == "horn":
    losses.custom_loss = hornLoss
    autoencoder.compile(optimizer = Adam(learning_rate), loss = hornLoss)
if args.metric == "morisita":
    losses.custom_loss = morisitaLoss
    autoencoder.compile(optimizer = Adam(learning_rate), loss = morisitaLoss)


if len(os.listdir(weight_dir)) == 0:  
    filepath = weight_dir + "train{loss:.3f}_val{val_loss:.3f}.hdf5" 
    checkpoint = ModelCheckpoint(filepath, monitor='val_loss', verbose=1, save_best_only=True, period = period)
    callbacks_list = [checkpoint]

    #train model
    print(autoencoder.summary())
    print(weight_dir)

    H = autoencoder.fit(np.array(x_train), np.array(y_train),
                    validation_data = (np.array(x_val), np.array(y_val)),
                    epochs = epochs, batch_size = batch_size, initial_epoch = initial_epoch, callbacks = callbacks_list)
    plt.figure(figsize=(8,8))
    plt.plot(H.history['loss'], label='training loss')
    plt.plot(H.history['val_loss'], label='validation loss')
    plt.legend()
    plt.title('Train/validation loss')
    plt.show()
