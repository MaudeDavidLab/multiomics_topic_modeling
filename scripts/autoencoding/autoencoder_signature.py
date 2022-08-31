## dataset
#import tensorflow.keras as keras
## for Model definition/training
from tensorflow.keras.models import Model, load_model
from tensorflow.keras.layers import Input, Flatten, Dense, concatenate,  Dropout
from tensorflow.keras.optimizers import Adam
import tensorflow.keras.backend as K
import tensorflow.keras.losses as losses



from tensorflow.keras.utils import plot_model
from tensorflow.keras.callbacks import ModelCheckpoint
from keras.models import model_from_json

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
from sklearn.preprocessing import scale

import argparse
from os import listdir
import os

def create_base_network(input_shape, embedding_size):
    """
    Base network to be shared (eq. to feature extraction).
    """
    input_vec = Input(shape= input_shape)
    x = Flatten()(input_vec)
    
    encoded = Dense(embedding_size, activation='relu', name = "encoder")(x)
   

    decoded = Dense(input_shape, activation = 'relu', name = "decoder")(encoded)
    
    autoencoder = Model(inputs=input_vec, outputs= decoded)
    
    encoder = Model(inputs = input_vec, outputs = encoded)
    
    encoded_input = Input(shape=(embedding_size,))
    decoder_layer = autoencoder.layers[-1]
    decoder = Model(inputs = encoded_input, outputs = decoder_layer(encoded_input))
    
    return autoencoder, encoder, decoder


#parse parameters
parser = argparse.ArgumentParser()
parser.add_argument('--train_data_file', type = str, dest = 'train_data_file')
parser.add_argument('--val_data_file', type = str, dest = 'val_data_file')
parser.add_argument('--test_data_file', type = str, dest = 'test_data_file')
parser.add_argument('--epochs', type = int, dest = 'epochs')
parser.add_argument('--lr', type = float, dest = 'lr')
parser.add_argument('--embedding_size', type = int, dest = 'embedding_size')
parser.add_argument('--metric', type = str, dest = 'metric')
parser.add_argument('--omic', type = str, dest = 'omic')

args = parser.parse_args()
print(args)

def encode(data, model_file):
    input_shape = data.shape[1]
    model = load_model(model_file)
    input_vec = Input(shape= input_shape)
    x = Flatten()(input_vec)
    encoded = model.get_layer('dense')(x)
    encoder = Model(inputs = input_vec, outputs = encoded)
    encoded = pd.DataFrame(encoder.predict(data), index = data.index)
    return(encoded)



def getWeightDir(omic, dim, metric):
    #Name file
    weight_dir = "weights/autoencoder/" + omic + "/" + metric + "/" + str(dim) + "dim/"
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

losses.custom_loss = morisitaLoss

#load data
x_train = pd.read_csv(args.train_data_file, index_col = 0)
x_val = pd.read_csv(args.val_data_file, index_col = 0)
x_test = pd.read_csv(args.test_data_file, index_col = 0)
x_train = x_train + 1
x_val = x_val + 1
x_test = x_test + 1

#def dropLowStd(data):
#    return(data.loc[:, np.std(data, axis = 0) > 0])

#def normalize(data):
#    return((data - np.mean(data, axis = 0)) / (np.var(data, axis = 0)))
#x_train = normalize(dropLowStd(x_train))
#x_val = normalize(dropLowStd(x_val))
#x_test = normalize(dropLowStd(x_test))
#features = x_train.columns.intersection(x_val.columns)
#features = features.intersection(x_test.columns)
#x_train = x_train.loc[:, features]
#x_val = x_val.loc[:, features]
#x_test = x_test.loc[:, features]

print(x_train.shape)

input_shape = x_train.shape[1]
batch_size = 100
period = 100
initial_epoch = 0

autoencoder, encoder, decoder = create_base_network(input_shape, args.embedding_size)

if args.metric == "euclidean":
    autoencoder.compile(optimizer = Adam(args.lr), loss = 'mean_squared_logarithmic_error')
if args.metric == "braycurtis":
    losses.custom_loss = brayCurtisLoss
    autoencoder.compile(optimizer = Adam(args.lr), loss = brayCurtisLoss)
if args.metric == "manhattan":
    autoencoder.compile(optimizer = Adam(args.lr), loss = 'mean_absolute_error')
if args.metric == "horn":
    losses.custom_loss = hornLoss
    autoencoder.compile(optimizer = Adam(args.lr), loss = hornLoss)
if args.metric == "morisita":
    losses.custom_loss = morisitaLoss
    autoencoder.compile(optimizer = Adam(args.lr), loss = morisitaLoss)

weights_dir = getWeightDir(args.omic, args.embedding_size, args.metric)
if len(os.listdir(weights_dir)) == 0:
    filepath =  weights_dir + "ep{epoch:02d}_train{loss:.3f}_val{val_loss:.3f}.hdf5" 
    checkpoint = ModelCheckpoint(filepath, monitor='val_loss', verbose=1, save_best_only=False, period = period)
    callbacks_list = [checkpoint]

    #train model
    print(autoencoder.summary())

    epochs = args.epochs + initial_epoch
    H = autoencoder.fit(np.array(x_train), np.array(x_train),
                    validation_data = (np.array(x_val), np.array(x_val)),
                    epochs = epochs, batch_size = batch_size, initial_epoch = initial_epoch, callbacks = callbacks_list)
