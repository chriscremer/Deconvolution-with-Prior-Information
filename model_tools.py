
import numpy as np


from sklearn.decomposition import PCA
from sklearn.decomposition import NMF
from sklearn.decomposition import FastICA

from numpy.linalg import inv
from numpy.linalg import pinv

import weight_tools as wt

import global_variables as gv

import random

def init_model(init_type, numb_model_subpops, X):


	if init_type == 0:
		Z = np.random.rand(numb_model_subpops, len(X[0]))
	elif init_type == 1:
		# pca = PCA(n_components=numb_model_subpops)
		# pca.fit(X.T)
		# Z = pca.transform(X.T).T
		Z = X[random.sample(range(len(X)), numb_model_subpops)]
	else:
		ica = FastICA(n_components=numb_model_subpops)
		ica.fit(X.T)
		Z = ica.transform(X.T).T		

	#print 'Z shape ' + str(Z.shape)

	T = np.identity(len(Z))
	#print 'T shape ' + str(T.shape)

	TZ = np.dot(T, Z)
	#print 'TZ shape ' + str(TZ.shape)

	return TZ


def optimize_model(possible_Ws, TZ):

	for i in range(50):

		#print 'Iter ' + str(i)

		#print 'selecting W'
		#W = select_w(X, possible_Ws, TZ)
		W = wt.select_w_parallel(possible_Ws)
		X_hat = np.dot(W, TZ)
		norm = np.linalg.norm(gv.X - X_hat)
		#print '         	Norm ' + str(norm)


		#print 'optimizig TZ'
		TZ = np.dot(pinv(np.dot(W.T,W)), np.dot(W.T, gv.X))

		gv.set_current_TZ(TZ)

		new_X_hat = np.dot(W, TZ)
		new_norm = np.linalg.norm(gv.X - new_X_hat)
		#print '         	Norm ' + str(new_norm)

		if norm == new_norm:
			print 'Number of iterations to optimize= ' + str(i)
			break

	return W, TZ