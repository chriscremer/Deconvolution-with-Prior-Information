


import numpy as np

from scipy.optimize import nnls

from nnls import nnls_with_prior

import itertools

from cvxopt import matrix as mtrx
from cvxopt import solvers


def init_Z(X, n_components):

	#Initialize Z
	# Z = X[random.sample(range(len(X)), self.n_components)]
	Z = []
	used = []
	for i in range(n_components):
		if len(Z) == 0:
			#so that they start at the same place
			rand = 0
			# rand = random.randint(0,len(X)-1)
			Z.append(X[rand])
			used.append(rand)
		else:
			#select sample that is furthest from mean of already selected samples
			Z_array = np.array(Z)
			mean1 = [np.mean(x) for x in Z_array.T]
			furthest_samp = -1
			furthest_dist = -1
			for samp in range(len(X)):
				if samp in used:
					continue
				else:
					distance = sum(abs(X[samp] - mean1))
				if distance > furthest_dist or furthest_dist == -1:
					furthest_dist = distance
					furthest_samp = samp
			Z.append(X[furthest_samp])
			used.append(furthest_samp)
	Z = np.array(Z)
	return Z
	

def use_all_components(X, W, Z):

	#If no assignments to a profile, give it the one with the worst norm
	count = 0 #Numb of zero profiles found
	calculated = 0 #To know if I've aleady calculated the difference matrix
	for i in range(len(W.T)):
		if sum(W.T[i]) == 0:
			if calculated == 0:
				calculated = 1
				dif = X - np.dot(W, Z)
				norms = np.linalg.norm(dif, axis=1)
				sort_indexes = np.argsort(norms)[::-1] #descending order
			for j in range(len(W[sort_indexes[count]])):
				if j == i:
					W[sort_indexes[count]][j] = 1.0
				else:
					W[sort_indexes[count]][j] = 0.0
			count += 1

	return W


def optimize_z_using_nnls(X, W):

	#Optimize Z
	Z = []
	for d in range(len(X.T)):
		Z_d = nnls(W, X.T[d])[0]
		Z.append(Z_d)
	Z = np.array(Z).T

	return Z


def optimize_z_and_w_using_nnls_given_X_and_Z(X, Z, max_iters, tol):

	converged = 0
	for iter1 in range(max_iters):

		#Optimize W
		W = []
		for i in range(len(X)):
			W_i = nnls(Z.T, X[i])[0]
			W.append(W_i)
		W = np.array(W)

		W = use_all_components(X, W, Z)

		#scale W for each sample so that sum = 1
		for i in range(len(W)):
			W[i] = W[i]/sum(W[i])

		norm = np.linalg.norm(X - np.dot(W, Z))
		# print 'norm ' + str(norm)

		#Optimize Z
		Z = []
		for d in range(len(X.T)):
			Z_d = nnls(W, X.T[d])[0]
			Z.append(Z_d)
		Z = np.array(Z).T

		new_norm = np.linalg.norm(X - np.dot(W, Z))
		# print 'new_norm ' + str(new_norm)

		if (norm - new_norm) < tol:
			# print '# iters until optimized= ' + str(iter1)
			converged = 1
			break

	if converged == 0:
		print 'Did not converge1'

	return W, Z


def assign_W_to_top_5_components(X, freqs, Z):

	new_W = []
	for samp in range(len(X)):

		#get 5 closest profiles to the sample
		dist_to_profiles = []
		for p in range(len(Z)):
			# dist_to_profiles.append(sum(abs(X[samp] - Z[p])))
			dist_to_profiles.append(np.linalg.norm(X[samp] - Z[p]))
		index_of_closest_profiles = np.argsort(dist_to_profiles)[:5]
		Z_reduced = Z[index_of_closest_profiles]

		#get non zero freqs
		non_zero_freqs = []
		for freq in range(len(freqs[samp])):
			if freqs[samp][freq] > 0:
				non_zero_freqs.append(freqs[samp][freq])

		#if more than 5 freqs, remove lower freqs and re-normalize
		if len(non_zero_freqs) > 5:
			sorted_freqs = sorted(non_zero_freqs)[::-1]
			top_5 = np.array(sorted_freqs[:5])
			top_5 = top_5 / sum(top_5)
			non_zero_freqs = list(top_5)

		#if less than 5, add some zeros
		if len(non_zero_freqs) < 5:
			while len(non_zero_freqs) < 5:
				non_zero_freqs.append(0.0)

		#find permuations and remove the duplicates
		perms = list(set(itertools.permutations(non_zero_freqs)))
		best_perm = -1
		best_perm_error = -1
		for perm in range(len(perms)):

			perm1 = np.array(perms[perm])
			X_i_hat = np.dot(perm1, Z_reduced)
			# perm_error = sum(abs(X_i_hat - X[samp]))
			perm_error = np.linalg.norm(X_i_hat - X[samp])

			if perm_error < best_perm_error or best_perm_error == -1:
				best_perm_error = perm_error
				best_perm = perm1

		#convert indexes back to the W vector indexes
		W_i = np.zeros(len(Z))
		for index in range(len(index_of_closest_profiles)):
			W_i[index_of_closest_profiles[index]] = best_perm[index]
		new_W.append(W_i)

	W = np.array(new_W)

	W = use_all_components(X, W, Z)

	return W


def fit_prior_freqs(W, freqs):


	#TODO
	#match the Ws

	new_W = []
	for samp in range(len(W)):
		norm_matrix = []
		for freq in range(len(freqs[samp])):
			this_samp = []
			for w in range(len(W[samp])):
				this_samp.append(abs(freqs[samp][freq] - W[samp][w]))
			norm_matrix.append(this_samp)
		# for list1 in norm_matrix:
		# 	print str(['%.2f' % elem for elem in list1])
		from munkres import Munkres, print_matrix
		m = Munkres()
		indexes = m.compute(norm_matrix)
		indexes2 = [x[1] for x in indexes]

		new_W_i = [0]*len(W[samp])
		for i in range(len(indexes2)):
			new_W_i[indexes2[i]] = freqs[samp][i]
		new_W.append(new_W_i)

	W = np.array(new_W)

	return W


def nnls_with_prior_for_each_samp(X,W,Z,lambda1):

	W_new = []
	for i in range(len(X)):
		# print i
		# print W[i]
		W_i = nnls_with_prior(Z.T,X[i],W[i],lambda1)
		W_new.append(W_i)
		# print W_i

	W = np.array(W_new)

	W = use_all_components(X, W, Z)

	#scale W for each sample so that sum = 1
	for i in range(len(W)):
		W[i] = W[i]/sum(W[i])

	return W


def optimize_z_and_w_using_nnls_and_LP_given_X_and_Z(X, Z, max_iters, tol):


	N = len(X)
	K = len(Z)
	D = len(X[0])


	c = [0.]*K
	c.extend([1.]*D)
	c.extend([1.]*K)
	c = np.array(c)
	print 'c shape ' + str(c.shape)
	# print c
	c = mtrx(c)
	# c = matrix(c, (c.shape[0], 1))
	# print c.size[0]
	# print c.size[1]
	# print type(c)
	



	# print Z.shape
	G_row = np.concatenate((Z.T, np.identity(D)), axis=1) 
	G_row = np.concatenate((G_row, np.zeros((D,K))), axis=1)
	# print G_row.shape
	G_row2 = np.concatenate((-Z.T, np.identity(D)), axis=1) 
	G_row2 = np.concatenate((G_row2, np.zeros((D,K))), axis=1)
	# print G_row2.shape
	G = np.concatenate((G_row,G_row2),axis=0)
	# print G.shape

	G_row3 = np.concatenate((np.identity(K), np.zeros((K,D))), axis=1) 
	G_row3 = np.concatenate((G_row3, -np.identity(K)), axis=1)
	# print G_row3.shape
	G = np.concatenate((G,G_row3),axis=0)
	# print G.shape

	G_row4 = np.concatenate((-np.identity(K), np.zeros((K,D))), axis=1) 
	G_row4 = np.concatenate((G_row4, -np.identity(K)), axis=1)
	# print G_row4.shape
	G = np.concatenate((G,G_row4),axis=0)
	# print G.shape

	G_row5 = np.concatenate((-np.identity(K), np.zeros((K,D))), axis=1) 
	G_row5 = np.concatenate((G_row5, np.zeros((K,K))), axis=1)
	# print G_row5.shape
	G = np.concatenate((G,G_row5),axis=0)
	print 'G shpae ' + str(G.shape)
	G = mtrx(G)

	converged = 0
	for iter1 in range(max_iters):

		#Optimize W
		W = []
		for i in range(len(X)):
			# W_i = nnls(Z.T, X[i])[0]

			h = list(X[i])
			h.extend(-X[i])
			h.extend([0]*K)
			h.extend([0]*K)
			h.extend([0]*K)
			h = np.array(h)

			h = mtrx(h)

			# print type(c)
			# print (mtrx)
			# print c.typecode
			# print c.size[1]

			solvers.options['show_progress'] = False
			sol = solvers.lp(c, G, h)
			# print sol
			W_i = np.array([x for x in sol['x'][:5]])
			# print W_i
			# print 'done ' + str(i)

			W.append(W_i)
		W = np.array(W)

		print W.shape

		W = use_all_components(X, W, Z)

		#scale W for each sample so that sum = 1
		for i in range(len(W)):
			W[i] = W[i]/sum(W[i])

		norm = np.linalg.norm(X - np.dot(W, Z))
		# print 'norm ' + str(norm)

		#Optimize Z
		Z = []
		for d in range(len(X.T)):
			Z_d = nnls(W, X.T[d])[0]
			Z.append(Z_d)
		Z = np.array(Z).T

		new_norm = np.linalg.norm(X - np.dot(W, Z))
		# print 'new_norm ' + str(new_norm)

		if (norm - new_norm) < tol:
			# print '# iters until optimized= ' + str(iter1)
			converged = 1
			break

	if converged == 0:
		print 'Did not converge1'

	return W, Z




class Abs_Error_Abs_Reg():
	"""
	Deconvolution Incorporating Frequency Information (Difi)
	A simple linear model with latent variables.
	"""


	def __init__(self, n_components=3, tol=1e-3, max_iter=200, rand_inits=1, lambda1=None):
		self.n_components = n_components
		self.tol = tol
		self.max_iter = max_iter
		self.rand_inits = rand_inits
		self.lambda1 = lambda1
		self.components_ = None
		self.norm = None
		self.W = None
		self.Z = None


	def fit(self, X):
		"""
		Fit the DIFI model to X

		Parameters
		----------
		X : array-like, shape (n_samples, n_features)
			Training data.

		Returns
		-------
		self
		"""

		best_norm = -1
		best_components = -1
		for inits in range(self.rand_inits):
			# print 'Rand Init ' + str(inits)

			#Initialize Z
			Z = init_Z(X, self.n_components)

			#Optimize model with alternating NNLS
			# W, Z = optimize_z_and_w_using_nnls_given_X_and_Z(X, Z, self.max_iter, self.tol)

			W, Z = optimize_z_and_w_using_nnls_and_LP_given_X_and_Z(X, Z, self.max_iter, self.tol)


			#Select W using top 5
			# W = assign_W_to_top_5_components(X, freqs, Z)

			#Inject freq data by fitting it to W
			# W = fit_prior_freqs(W, freqs)

			#Optimize model with alternating NNLS
			# W, Z = optimize_z_and_w_using_nnls_given_X_and_W(X, W, self.max_iter, self.tol)

			#Optimize model with alternating NNLS
			# W, Z = optimize_z_and_w_using_nnls_given_X_and_W(X, W, self.max_iter, self.tol)

			#Allow to deviate
			# W = regression_with_prior(X, W, Z)


			# print np.linalg.norm(X - np.dot(W, Z))

			# print 'W1 ' + str(W[0])

			# to_beat = -1
			# for i in range(self.max_iter):
			# 	# print 'iter ' + str(i)

			# 	#Inject freq data by fitting it to W
			# 	# W = fit_prior_freqs(W, freqs)
			# 	# W = assign_W_to_top_5_components(X, freqs, Z)
			# 	# print 'W2 ' + str(W[0])
				
			# 	# print np.linalg.norm(X - np.dot(W, Z))

			# 	#Let it deviate
			# 	if self.lambda1 == None:
			# 		lambda1 = len(X[0])
			# 	else:
			# 		lambda1 = self.lambda1

			# 	W = nnls_with_prior_for_each_samp(X,W,Z,lambda1)

			# 	W = lp_for_each_samp(X,Z,lambda1)


			# 	# print 'W3 ' + str(W[0])

			# 	norm = np.linalg.norm(X - np.dot(W, Z))
			# 	if to_beat - norm < self.tol and to_beat != -1:
			# 		# print 'done'
			# 		break
			# 	to_beat = norm

			# 	# print np.linalg.norm(X - np.dot(W, Z))
			# 	best_W = W

			# 	#Optimize Z
			# 	Z = optimize_z_using_nnls(X, W)

			# 	# print np.linalg.norm(X - np.dot(W, Z))

			# 	# print
			# 	if i == (self.max_iter - 1):
			# 		print 'max reached!!'

			# print str(i) + ' iterations'


		# print Z.T[:5].T

		self.W = W
		self.Z = Z

		return self


	def transform(self, X, freqs):
		"""Solve for W given Z and X

		Parameters
		----------
		X : array-like, shape (n_samples, n_features)

		Returns
		-------
		X_new : array-like, shape (n_samples, n_components)
			The latent variables of X.
		"""

		print 'This needs to be implemented'

		return None
