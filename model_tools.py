
import numpy as np


from sklearn.decomposition import PCA
from sklearn.decomposition import NMF
from sklearn.decomposition import FastICA

from numpy.linalg import inv
from numpy.linalg import pinv

import weight_tools as wt


def init_model(init_type, numb_model_subpops, samps):


	if init_type == 0:
		Z = np.random.rand(numb_model_subpops, len(samps[0]))
	elif init_type == 1:
		pca = PCA(n_components=numb_model_subpops)
		pca.fit(samps.T)
		Z = pca.transform(samps.T).T
	else:
		ica = FastICA(n_components=numb_model_subpops)
		ica.fit(samps.T)
		Z = ica.transform(samps.T).T		

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
		norm = np.linalg.norm(X - X_hat)
		#print '         	Norm ' + str(norm)


		#print 'optimizig TZ'
		TZ = np.dot(pinv(np.dot(W.T,W)), np.dot(W.T, X))
		new_X_hat = np.dot(W, TZ)
		new_norm = np.linalg.norm(X - new_X_hat)
		#print '         	Norm ' + str(new_norm)

		if norm == new_norm:
			break

	return W, TZ