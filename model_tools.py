
import numpy as np


from sklearn.decomposition import PCA
from sklearn.decomposition import NMF
from sklearn.decomposition import FastICA

from sklearn.linear_model import Lasso

from sklearn.covariance import EmpiricalCovariance as ec

from numpy.linalg import inv
from numpy.linalg import pinv
from numpy.linalg import det

from scipy.stats import multivariate_normal as mn
from scipy.optimize import nnls

import weight_tools as wt

import global_variables as gv

import random
import math as m

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

	# T = np.identity(len(Z))
	# TZ = np.dot(T, Z)

	#init W
	# W = []
	# n_profiles = len(Z)
	# for i in range(len(X)):
	# 	# weight_vector = np.zeros(n_profiles)
	# 	w_i = [0]*n_profiles
	# 	w_i[i%n_profiles] = 1.0
	# 	W.append(w_i)
	# W = np.array(W)

	return Z


def optimize_model():

	gv.set_sample_norms([-1]*len(gv.X))
	gv.set_current_W(np.zeros((len(gv.X), len(gv.TZ))))
	max_iters = 200
	printed = 0
	best_W = None
	best_TZ = None
	best_norm = -1
	for itera in range(max_iters):

		#print 'Iter ' + str(i)

		#print 'selecting W'
		# W = select_w(X, possible_Ws, TZ)

		#this is the restricted way
		# W = wt.select_w_parallel()

		#this is the non-constrained way
		# W = np.dot(pinv(np.dot(gv.TZ,gv.TZ.T)), np.dot(gv.TZ, gv.X.T))
		# W = W.T

		#this is the new way
		# W = wt.select_w_parallel_NEW()
		# print np.array(W).shape

		#LASSO
		# problem since ... i tried working it out. aaahh
		# W = []
		# for i in range(len(gv.X)):
		# 	model = Lasso()
		# 	model.fit(gv.TZ.T, gv.X[i].T)
		# 	W_i = model.coef_
		# 	W.append(W_i)
		# W = np.array(W)
		# print np.array(W).shape
	
		#My new least squares but Wi done individually. 
		# W = []
		# for i in range(len(gv.X)):
		# 	# Wi = (XiZt)*(ZZt)-1
		# 	W_i = np.dot(np.dot(gv.X[i], gv.TZ.T), pinv(np.dot(gv.TZ, gv.TZ.T))) 
		# 	W.append(W_i)
		# W = np.array(W)	

		#Same but with L2 reg
		# W = []
		# for i in range(len(gv.X)):
		# 	# Wi = (XiZt)*(ZZt+I)-1
		# 	W_i = np.dot(np.dot(gv.X[i], gv.TZ.T), pinv(np.dot(gv.TZ, gv.TZ.T) + np.identity(len(gv.TZ)))) 
		# 	W.append(W_i)
		# W = np.array(W)

		# nnls
		W = []
		for i in range(len(gv.W)):
			W_i = nnls(gv.TZ.T, gv.X[i])[0]
			W.append(W_i)
		W = np.array(W)
		
		#if no assignments to a profile, give it the one with the worst norm
		for i in range(len(W.T)):
			if sum(W.T[i]) == 0:
				X_hat = np.dot(W, gv.TZ)
				dif = gv.X - X_hat
				norms = np.linalg.norm(dif, axis=1)
				worst_sample = list(norms).index(max(norms))
				for j in range(len(W[worst_sample])):
					if j == i:
						W[worst_sample][j] = 1.0
					else:
						W[worst_sample][j] = 0.0

		for i in range(len(W.T)):
			if sum(W.T[i]) == 0:
				print 'ZERO'

		gv.set_current_W(W)
		X_hat = np.dot(gv.W, gv.TZ)
		dif = gv.X - X_hat
		norm = np.linalg.norm(dif)
		gv.norms = np.linalg.norm(dif, axis=1)


		#print 'optimizig TZ'
		#TZ = np.dot(pinv(np.dot(gv.W.T,gv.W)), np.dot(gv.W.T, gv.X))
		TZ = []
		for d in range(len(gv.X.T)):
			# TZ_d = np.reshape(nnls(W, gv.X.T[d])[0], (len(gv.TZ),1))
			TZ_d = nnls(W, gv.X.T[d])[0]
			TZ.append(TZ_d)
		TZ = np.array(TZ).T
		print TZ.shape


		# TZ_x = np.reshape(nnls(W, gv.X.T[0])[0], (3,1))
		# TZ_y = np.reshape(nnls(W, gv.X.T[1])[0], (3,1))
		# TZ = np.concatenate((TZ_x, TZ_y), axis=1)

		gv.set_current_TZ(TZ)
		new_X_hat = np.dot(gv.W, gv.TZ)
		dif = gv.X - new_X_hat
		new_norm = np.linalg.norm(dif)
		gv.norms = np.linalg.norm(dif, axis=1)





		print 'Opt ' + str(itera) + ' norm after W ' + str(norm)  + ' norm after Z ' + str(new_norm) 

		if best_norm > new_norm or best_norm == -1:
			best_norm = new_norm
			best_W = gv.W
			best_TZ = gv.TZ

		if (norm - new_norm) < .0000000001:
			print '# iters until optimized= ' + str(itera)
			printed = 1
			break

		
	if printed == 0:
		print '# iters went to limit: ' + str(max_iters)
	return best_W, best_TZ






def optimize_model_with_loglike(TZ):

	gv.set_sample_norms([-1]*len(gv.X))
	gv.set_current_W(np.zeros((len(gv.X), len(gv.TZ))))
	max_iters = 200
	printed = 0
	best_W = None
	best_TZ = None
	best_norm = -1
	for itera in range(max_iters):

		

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
		dif = gv.X - X_hat
		norm = np.linalg.norm(dif)
		gv.norms = np.linalg.norm(dif, axis=1)


		# cov_est = ec(store_precision=False)
		# cov_est.fit(dif)
		# print cov_est.covariance_
		# print cov_est.score(gv.X)
		# asfd
		# np.diag(np.diag(cov))

		#ln p(Xn|X_hat, cov)
		# cov = (1.0/len(gv.X)) * np.dot((dif).T, (dif)) 
		cov = np.identity(len(dif[0]))
		# cov = np.diag(np.diag(cov))
		inv_cov = pinv(cov)
		# print inv_cov[0][:10]
		cov_det = det(cov)
		# print cov_det
		n_samps = len(gv.X)
		n_dims = len(gv.X[0])
		first_term = (-1./2.)*n_samps*n_dims*m.log(2*m.pi) 
		# second_term = (-1./2.)*n_samps*m.log(cov_det)
		second_term = (-1./2.)*n_samps

		sum1 = 0
		for i in range(n_samps):
			
			val = np.dot(np.dot(dif[i].T, inv_cov), dif[i]) 
			sum1 += val 
			# print val

		# log_likelihood = first_term + second_term + (-1./2.)*sum1
		log_likelihood = (-1./2.)*sum1
		print 'Log like ' + str(log_likelihood) + ' after W'
		





		# 	like = mn.logpdf(gv.X[i], mean=X_hat[i], cov=cov, allow_singular=True)
		# 	# print like
		# 	sum1 += like
		# print 'Log Like ' + str(sum1)


		#W = np.dot(pinv(np.identity(len(TZ))*0.01 + np.dot(TZ,TZ.T)), np.dot(TZ, gv.X.T))
		#print W.shape

		#print 'optimizig TZ'
		TZ = np.dot(pinv(np.dot(W.T,W)), np.dot(W.T, gv.X))
		gv.set_current_TZ(TZ)
		new_X_hat = np.dot(W, TZ)
		dif = gv.X - new_X_hat
		new_norm = np.linalg.norm(dif)
		gv.norms = np.linalg.norm(dif, axis=1)






		#ln p(Xn|X_hat, cov)
		# cov = (1.0/len(gv.X)) * np.dot((dif).T, (dif)) 
		cov = np.identity(len(cov))
		# cov = np.diag(np.diag(cov))
		inv_cov = pinv(cov)
		cov_det = det(cov)
		# print cov_det
		first_term = (-1./2.)*n_samps*n_dims*m.log(2*m.pi) 
		# second_term = (-1./2.)*n_samps*m.log(cov_det)
		second_term = (-1./2.)*n_samps

		sum1 = 0
		for i in range(n_samps):
			
			sum1 += np.dot(np.dot(dif[i].T, inv_cov), dif[i]) 

		# log_likelihood = first_term + second_term + (-1./2.)*sum1
		log_likelihood = (-1./2.)*sum1
		print 'Log like ' + str(log_likelihood) + ' after Z'



		# cov = (1.0/len(gv.X)) * np.dot((dif).T, (dif)) 
		# np.diag(np.diag(cov))
		# sum1 = 0
		# for i in range(len(gv.X)):
		# 	like = mn.logpdf(gv.X[i], mean=new_X_hat[i], cov=cov, allow_singular=True)
		# 	# print like
		# 	sum1 += like
		# print 'Log Like ' + str(sum1)



		print 'Opt ' + str(itera) + ' norm after W ' + str(norm) + ' norm after Z ' + str(new_norm)



		if best_norm > new_norm or best_norm == -1:
			best_norm = new_norm
			best_W = W
			best_TZ = TZ

		if norm == new_norm:
			print '# iters until optimized= ' + str(itera)
			printed = 1
			break

		
	if printed == 0:
		print '# iters went to limit: ' + str(max_iters)
	return best_W, best_TZ



def optimize_model_with_ls(TZ):

	gv.set_sample_norms([-1]*len(gv.X))
	gv.set_current_W(np.zeros((len(gv.X), len(gv.TZ))))
	max_iters = 200
	printed = 0
	best_W = None
	best_TZ = None
	best_norm = -1
	for itera in range(max_iters):

		#print 'Iter ' + str(i)

		#print 'selecting W'
		# W = select_w(X, possible_Ws, TZ)

		#this is the restricted way
		# W = wt.select_w_parallel()

		#this is the non-constrained way
		# W = np.dot(pinv(np.dot(TZ,TZ.T)), np.dot(TZ, gv.X.T))
		# W = W.T

		print 'freqs ' + str(gv.freqs[0])

		#this is the new way
		W = wt.select_w_parallel_NEW()

		print 'old W ' + str(W[0])

		gv.set_current_W(W)
		X_hat = np.dot(W, TZ)
		dif = gv.X - X_hat
		norm = np.linalg.norm(dif)
		gv.norms = np.linalg.norm(dif, axis=1)


		ls = 0
		if ls == 1:
			W_ls = []
			lambda1 = 1.0*len(gv.X[0])/len(gv.TZ)
			print 'lambda ' + str(lambda1)
			Z_dot_Zt = np.dot(gv.TZ, gv.TZ.T)
			ZZ_plus_lambda = Z_dot_Zt + (lambda1*np.identity(len(Z_dot_Zt)))
			ZZ_plus_lambda_inv = pinv(ZZ_plus_lambda)

			for i in range(len(W)):
				W_i = np.dot((np.dot(gv.X[i], gv.TZ.T) + (lambda1*W[i])), (ZZ_plus_lambda_inv))
				W_ls.append(W_i)
			W = np.array(W_ls)
			gv.set_current_W(W)
			X_hat = np.dot(W, TZ)
			dif = gv.X - X_hat
			norm = np.linalg.norm(dif)
			gv.norms = np.linalg.norm(dif, axis=1)
			print 'new W ' + str(W[0])


		if ls == 2:
			W_ls = []
			for i in range(len(W)):
				#get the indexes of the non-zero weights
				non_zero_indexes = []
				for j in range(len(W[i])):
					if W[i][j] > 0.0:
						non_zero_indexes.append(j)
				temp_TZ = gv.TZ[non_zero_indexes]
				# print 'temp TZ shape ' + str(temp_TZ.shape)
				new_W = np.dot(pinv(np.dot(temp_TZ,temp_TZ.T)), np.dot(temp_TZ, gv.X[i].T))
				new_W = new_W.T
				# print 'new W shape ' + str(new_W.shape)
				#check to see if there are any negatives, if there are, dont use it, for now
				neg = 0
				for j in range(len(new_W)):
					if new_W[j] < 0.0:
						neg = 1
						continue
				if neg == 0:
					W_i = np.zeros(len(W[i]))
					for j in range(len(new_W)):
						W_i[non_zero_indexes[j]] = new_W[j]
					W_ls.append(W_i)
				else:
					W_ls.append(W[i])
			W = np.array(W_ls)
			gv.set_current_W(W)
			X_hat = np.dot(W, TZ)
			dif = gv.X - X_hat
			norm = np.linalg.norm(dif)
			gv.norms = np.linalg.norm(dif, axis=1)
			print 'new W ' + str(W[0])


		#print 'optimizig TZ'
		TZ = np.dot(pinv(np.dot(W.T,W)), np.dot(W.T, gv.X))
		gv.set_current_TZ(TZ)
		new_X_hat = np.dot(W, TZ)
		dif = gv.X - new_X_hat
		new_norm = np.linalg.norm(dif)
		gv.norms = np.linalg.norm(dif, axis=1)

		print 'Opt ' + str(itera) + ' norm after W ' + str(norm) + ' norm after Z ' + str(new_norm)

		if best_norm > new_norm or best_norm == -1:
			best_norm = new_norm
			best_W = W
			best_TZ = TZ

		if abs(norm - new_norm) < 0.001:
			print '# iters until optimized= ' + str(itera)
			printed = 1
			break

		
	if printed == 0:
		print '# iters went to limit: ' + str(max_iters)
	return best_W, best_TZ



def match_profiles(TZ, real_profiles):

	#Get the norms between all pairs
	norm_matrix = []
	for learned_profile in range(len(TZ)):
		this_profile = []
		for hidden_profile in range(len(real_profiles)):
			this_profile.append(np.linalg.norm(TZ[learned_profile] - real_profiles[hidden_profile]))
		norm_matrix.append(this_profile)

	for list1 in norm_matrix:
		print str(['%.2f' % elem for elem in list1])

	#Hungarian Algorithm
	from munkres import Munkres, print_matrix
	m = Munkres()
	indexes = m.compute(norm_matrix)
	indexes2 = [x[1] for x in indexes]
	print indexes2
	# print_matrix(norm_matrix, msg='Lowest cost through this matrix:')
	total = 0
	for row, column in indexes:
	    value = norm_matrix[row][column]
	    value2 = np.linalg.norm(TZ[row] - real_profiles[column])
	    total += value
	    print '(%d, %d) -> %f' % (row, column, value2)
	print 'total cost: %f' % total


	return indexes2

