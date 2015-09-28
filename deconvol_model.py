
"""
Deconvolution
A model for deconvolution of mixed samples.
"""

# Author: Chris Cremer <ccremer@cs.toronto.edu>

import numpy as np
import random
from scipy.optimize import nnls

def init_Z(X, n_components):

	#Initialize Z
	# Z = X[random.sample(range(len(X)), self.n_components)]
	Z = []
	used = []
	for i in range(n_components):
		if len(Z) == 0:
			rand = random.randint(0,len(X)-1)
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
		


class Deconvol():
	"""
	Deconvol
	A simple linear model with latent variables.
	"""


	def __init__(self, n_components=3, tol=1e-3, max_iter=200, rand_inits=1):
		self.n_components = n_components
		self.tol = tol
		self.max_iter = max_iter
		self.rand_inits = rand_inits
		self.components_ = None
		self.norm = None


	def fit(self, X):
		"""
		Fit the Deconvol model to X

		Parameters
		----------
		X : array-like, shape (n_samples, n_features)
			Training data.

		Returns
		-------
		self
		"""

		best_norm = -1
		lowest_W = -1
		best_components = -1
		for inits in range(self.rand_inits):
			print 'Rand Init ' + str(inits)

			#Initialize Z
			Z = init_Z(X, self.n_components)

			#Optimize model
			converged = 0
			for iter1 in range(self.max_iter):

				#Optimize W
				W = []
				for i in range(len(X)):
					W_i = nnls(Z.T, X[i])[0]
					W.append(W_i)
				W = np.array(W)

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


				dif = X - np.dot(W, Z)
				norm = np.linalg.norm(dif)

				#Optimize Z
				Z = []
				for d in range(len(X.T)):
					Z_d = nnls(W, X.T[d])[0]
					Z.append(Z_d)
				Z = np.array(Z).T

				dif = X - np.dot(W, Z)
				new_norm = np.linalg.norm(dif)


				if (norm - new_norm) < self.tol:
					print '# iters until optimized= ' + str(iter1)
					converged = 1
					break

			if converged == 0:
				print 'Did not converge. Went to ' + str(self.max_iter) + ' iterations'

			# if new_norm < best_norm or best_norm == -1:
			W_sum = sum(sum(abs(W)))
			if W_sum < lowest_W or lowest_W == -1:
				best_norm = new_norm
				best_components = Z
				lowest_W = W_sum

		# self.W = W
		self.components_ = best_components
		self.norm = best_norm

		return self


	def transform(self, X):
		"""Solve for W given Z and X

		Parameters
		----------
		X : array-like, shape (n_samples, n_features)

		Returns
		-------
		X_new : array-like, shape (n_samples, n_components)
			The latent variables of X.
		"""

		W = []
		for i in range(len(X)):
			W_i = nnls(self.components_.T, X[i])[0]
			W.append(W_i)
		W = np.array(W)

		return W





class DIFI():
	"""
	Deconvolution Including Frequency Information (DIFI)
	A simple linear model with latent variables.
	"""


	def __init__(self, n_components=3, tol=1e-3, max_iter=200, rand_inits=1):
		self.n_components = n_components
		self.tol = tol
		self.max_iter = max_iter
		self.rand_inits = rand_inits
		self.components_ = None
		self.norm = None


	def fit(self, X, freqs):
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
		lowest_W = -1
		best_components = -1
		for inits in range(self.rand_inits):
			print 'Rand Init ' + str(inits)

			#Initialize Z
			Z = init_Z(X, self.n_components)

			#Optimize model
			converged = 0
			for iter1 in range(self.max_iter):

				#Optimize W

				#Check all permutations: #profiles P #freqs
				#Instead
				#Check permutations of closest x profiles: (#freqs + x) P #freqs













				W = []
				for i in range(len(X)):
					W_i = nnls(Z.T, X[i])[0]
					W.append(W_i)
				W = np.array(W)

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


				dif = X - np.dot(W, Z)
				norm = np.linalg.norm(dif)

				#Optimize Z
				Z = []
				for d in range(len(X.T)):
					Z_d = nnls(W, X.T[d])[0]
					Z.append(Z_d)
				Z = np.array(Z).T

				dif = X - np.dot(W, Z)
				new_norm = np.linalg.norm(dif)


				if (norm - new_norm) < self.tol:
					print '# iters until optimized= ' + str(iter1)
					converged = 1
					break

			if converged == 0:
				print 'Did not converge. Went to ' + str(self.max_iter) + ' iterations'

			# if new_norm < best_norm or best_norm == -1:
			W_sum = sum(sum(abs(W)))
			if W_sum < lowest_W or lowest_W == -1:
				best_norm = new_norm
				best_components = Z
				lowest_W = W_sum

		# self.W = W
		self.components_ = best_components
		self.norm = best_norm

		return self


	def transform(self, X):
		"""Solve for W given Z and X

		Parameters
		----------
		X : array-like, shape (n_samples, n_features)

		Returns
		-------
		X_new : array-like, shape (n_samples, n_components)
			The latent variables of X.
		"""

		W = []
		for i in range(len(X)):
			W_i = nnls(self.components_.T, X[i])[0]
			W.append(W_i)
		W = np.array(W)

		return W