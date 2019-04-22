import glob
import sklearn
import numpy as np
from numpy import genfromtxt
from sklearn import ensemble
import random

random.seed(3.14) # sklearn use np.random to generate random_state
num_tree=100
max_depth=5

## *.dat is gene x sample; transpose them
train=genfromtxt('train.dat',delimiter='\t').T
feature_train=genfromtxt('feature_train.dat',delimiter='\t').T 
feature_test=genfromtxt('feature_test.dat',delimiter='\t').T

num=feature_test.shape[0] # number of test samples
vec0=np.full([1,num],0,dtype='float64') # vector of zeros
pred=np.empty([0,num])

for i in range(train.shape[1]):

	if i%500==0:
		print(i)

	index=~np.isnan(train[:,i])
	X=feature_train[index,:]
	Y=train[index,i]

	if index.sum()>1: # if training samples < 2, can't rf and predict 0
		est = sklearn.ensemble.RandomForestRegressor(n_estimators=num_tree, max_depth=max_depth, random_state=0).fit(X,Y)
		vec=est.predict(feature_test)
	else:
		vec=vec0

	pred=np.vstack((pred,vec))


np.savetxt('prediction.dat', pred, delimiter='\t')

