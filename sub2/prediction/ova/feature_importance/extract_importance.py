import glob
import sklearn
import numpy as np
from numpy import genfromtxt
from sklearn import ensemble
import pickle
import os

files=os.listdir('model/')

for i in range(len(files)):

	filename='model/' + files[i]
	est= pickle.load(open(filename, 'rb'))
	filename='importance/' + files[i].split('.')[0] + '.txt'
	np.savetxt(filename, est.feature_importances_.T,delimiter='\t')


