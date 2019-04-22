import glob
import sklearn
import numpy as np
from numpy import genfromtxt
from sklearn import ensemble
#import random
import pickle

#random.seed(3.14) # sklearn use np.random to generate random_state
num_tree=100
max_depth=3

## *.dat is gene x sample; transpose them
train=genfromtxt('train.dat',delimiter='\t').T
feature_train=genfromtxt('feature_train.dat',delimiter='\t').T 
genename=genfromtxt('genename.txt',dtype='str')

for i in range(train.shape[1]):

	index=~np.isnan(train[:,i])
	X=feature_train[index,:]
	Y=train[index,i]

	if index.sum()>1: # if training samples < 2, can't rf and predict 0
		est = sklearn.ensemble.RandomForestRegressor(n_estimators=num_tree, max_depth=max_depth, random_state=0).fit(X,Y)
		filename = genename[i] + '.model'
		pickle.dump(est, open(filename, 'wb'))

	else:
		filename = genename[i] + '.fakemodel'
		pickle.dump(est, open(filename, 'wb'))


