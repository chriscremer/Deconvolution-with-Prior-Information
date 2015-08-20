
import numpy as np


from sklearn.decomposition import PCA
from sklearn.decomposition import NMF
from sklearn.decomposition import FastICA

from numpy.linalg import inv
from numpy.linalg import pinv


import weight_tools as wt

import global_variables as gv

import random
import math

def init_model(init_type, numb_model_subpops, X):


	if init_type == 'Random_Values':
		Z = np.random.rand(numb_model_subpops, len(X[0]))
	elif init_type == 'Random_Samples':
		#if more components than samples, then add some made up samples
		#else select some samples 
		if len(X) < numb_model_subpops:
			random_profiles = np.random.rand(numb_model_subpops - len(X), len(X[0]))
			#concatenate the two
			Z = np.concatenate((X, random_profiles), axis=0)
		else:
			Z = X[random.sample(range(len(X)), numb_model_subpops)]
	elif init_type == 'PCA':
		pca = PCA(n_components=numb_model_subpops)
		pca.fit(X.T)
		Z = pca.transform(X.T).T		
	elif init_type == 'ICA':
		ica = FastICA(n_components=numb_model_subpops)
		ica.fit(X.T)
		Z = ica.transform(X.T).T
	else:
		print 'WHAT INITIALIZATION IS THIS?'		

	T = np.identity(len(Z))
	TZ = np.dot(T, Z)

	return TZ


def optimize_model(TZ):

	gv.set_sample_norms([-1]*len(gv.X))
	gv.set_current_W(np.zeros((len(gv.X), len(gv.TZ))))
	max_iters = 200
	printed = 0
	best_W = None
	best_TZ = None
	best_norm = -1
	for i in range(max_iters):

		

		#print 'Iter ' + str(i)

		#print 'selecting W'
		# W = select_w(X, possible_Ws, TZ)

		#this is the restricted way
		# W = wt.select_w_parallel()

		#this is the non-constrained way
		# W = np.dot(pinv(np.dot(TZ,TZ.T)), np.dot(TZ, gv.X.T))
		# W = W.T

		#this is the new way
		W = wt.select_w_parallel_NEW()
		gv.set_current_W(W)
		X_hat = np.dot(W, TZ)
		norm = np.linalg.norm(gv.X - X_hat)
		gv.norms = np.linalg.norm(gv.X - X_hat, axis=1)

		#W = np.dot(pinv(np.identity(len(TZ))*0.01 + np.dot(TZ,TZ.T)), np.dot(TZ, gv.X.T))
		#print W.shape

		#print 'optimizig TZ'
		TZ = np.dot(pinv(np.dot(W.T,W)), np.dot(W.T, gv.X))
		gv.set_current_TZ(TZ)
		new_X_hat = np.dot(W, TZ)
		new_norm = np.linalg.norm(gv.X - new_X_hat)
		gv.norms = np.linalg.norm(gv.X - new_X_hat, axis=1)

		# dif = gv.X - new_X_hat
		# # qwer = []
		# qwert = []
		# for qwe in range(len(gv.X)):
		# 	# print np.linalg.norm(dif[qwe])
		# 	# qwer.append(np.linalg.norm(dif[qwe]))
		# 	for asd in range(len(gv.X[0])):
		# 		qwert.append(dif[qwe][asd]*1.0*dif[qwe][asd])

		# print len(qwert)
		# print sum(gv.norms)

		# print 'sums'
		# # print sum(qwert)
		# print math.sqrt(sum(qwert))
		# # print sum(qwer)
		
		# print new_norm
		# print np.linalg.norm(dif)


		# print 'other test'

		# a = np.arange(9) - 4
		# print np.linalg.norm(a)
		# sum1 = []
		# for tt in a:
		# 	sum1.append(tt * 1.0 * tt)

		# print sum(sum1)
		# print math.sqrt(sum(sum1))

		# adsf






		# for row in range(len(TZ)):
		# 	for col in range(len(TZ[row])):
		# 		if TZ[row][col] < 0:
		# 			print 'HEGATIVE'


		#print '         	Norm ' + str(new_norm)

		print 'Opt ' + str(i) + ' norm after W ' + str(norm) + ' norm after Z ' + str(new_norm)



		if best_norm > new_norm or best_norm == -1:
			best_norm = new_norm
			best_W = W
			best_TZ = TZ

		if norm == new_norm:
			print '# iters until optimized= ' + str(i)
			printed = 1
			break

		
	if printed == 0:
		print '# iters went to limit: ' + str(max_iters)
	return best_W, best_TZ