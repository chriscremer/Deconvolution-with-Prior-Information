


import numpy as np
from numpy.linalg import inv
from numpy import dot
from scipy.stats import linregress
from scipy.optimize import nnls
from sklearn.linear_model import LinearRegression


def nnls_with_prior(X, y, prior, lambda1):

	print1 = 0

	m = len(X)
	n = len(X[0])
	passive_set = []
	active_set = range(n)
	sol1 = [0]*n
	# lagr = dot(X.T,(y - dot(X, sol1)))
	grad = dot(X.T,(y - dot(X, sol1))) - lambda1*(sol1 - prior)

	if print1 == 1:
		print '1- sol1 ' + str(sol1)
		print '2- grad ' + str(grad)

	tol = .00001
	while (len(active_set) != 0 and max(grad[active_set]) > tol):

		if print1 == 1:
			print '3- sol1 ' + str(sol1)
			print '4- grad ' + str(grad)
			print 'max active set grad ' + str(max(grad[active_set]))
			print '5- active set ' + str(active_set)
			print '6- passive set ' + str(passive_set)

		# var_to_add = list(grad).index(max(grad[active_set]))
		var_to_add = list(grad).index(max(grad[active_set]))

		# max_grad = -1
		# max_grad_index = -1
		# for i in range(len(grad)):
		# 	if i in active_set:
		# 		if grad[i] > max_grad or max_grad == -1:
		# 			max_grad_index = i
		# 			max_grad = grad[i]


		if print1 == 1:
			print '6.1- var_to_add ' + str(var_to_add)
			# print 'new var to add!!! ' + str(max_grad_index)

		passive_set.append(var_to_add)
		active_set.remove(var_to_add)

		if print1 == 1:
			print '7- active set ' + str(active_set)
			print '8- passive set ' + str(passive_set)


		# print 'active set ' + str(active_set)
		# print 'passive set ' + str(passive_set)

		X_passive_only = X.T[passive_set].T
		prior_passive_only = prior.T[passive_set].T

		numer = dot(X_passive_only.T, y) + lambda1*prior_passive_only
		denom = dot(X_passive_only.T, X_passive_only) + lambda1*(np.identity(len(X_passive_only.T)))
		sol2 = dot(inv(denom), numer)
		# pseudo_inv = dot(inv(dot(X_passive_only.T, X_passive_only)), X_passive_only.T)
		# sol2 = dot(pseudo_inv, y)

		if print1 == 1:
			print '9- sol2 ' + str(sol2)

		#inner loop here
		while (min(sol2) <= 0):

			if print1 == 1:
				print 'inner'

			# print 'C ' + str(sol2)

			#add the zeros in, this will help with indexing
			sol2_full = [0]*n
			for i in range(len(passive_set)):
				sol2_full[passive_set[i]] = sol2[i]
			sol2_full = np.array(sol2_full)

			if print1 == 1:
				print 'sol2_full ' + str(sol2_full)

			#find min alpha
			a = -1
			for i in range(len(sol2_full)):
				# if i in passive_set:
				if sol2_full[i] < 0:
					a_i = sol1[i] / (sol1[i] - sol2_full[i])
					if a_i < a or a == -1:
						a = a_i
						if print1 == 1:
							print i
							print 'alpha ' + str(a)
							# print sol1[i]
							# print sol2_full[i]
							# print (sol1[i] - sol2_full[i])
			#not sure why only fnnls has this negative
			# a = -a
			sol1 = sol1 + a*(sol2_full - sol1)

			# print sol1

			#C4
			# passive_set = []
			# active_set = []
			# for i in range(len(sol1)):
			# 	if sol1[i] < tol:
			# 		active_set.append(i)
			# 	else:
			# 		passive_set.append(i)

			for i in passive_set:
				if sol1[i] < tol:
					if print1 == 1:
						print 'to active set ' + str(i)
						print 'because its sol1 value is ' + str(sol1[i])
					active_set.append(i)
					passive_set.remove(i)

			# print 'active set ' + str(active_set)
			# print 'passive set ' + str(passive_set)

			
			#C5
			X_passive_only = X.T[passive_set].T
			prior_passive_only = prior.T[passive_set].T
			numer = dot(X_passive_only.T, y) + lambda1*prior_passive_only
			denom = dot(X_passive_only.T, X_passive_only) + lambda1*(np.identity(len(X_passive_only.T)))
			sol2 = dot(inv(denom), numer)
			if print1 == 1:
				print 'inner- sol2 ' + str(sol2)
			#C6???
			# print 'D ' + str(sol2)


			# asd

			# pseudo_inv = dot(inv(dot(X_passive_only.T, X_passive_only)), X_passive_only.T)
			# sol2 = dot(pseudo_inv, y)

		# print '3'

		# print 'active set ' + str(active_set)
		# print 'passive set ' + str(passive_set)


		sol1 = [0]*n
		for i in range(len(passive_set)):
			sol1[passive_set[i]] = sol2[i]
		sol1 = np.array(sol1)
		# grad = dot(X.T,(y - dot(X, sol1)))
		grad = dot(X.T,(y - dot(X, sol1))) - lambda1*(sol1 - prior)

	if print1 == 1:
		print '98- sol1 ' + str(sol1)
		print '99- grad ' + str(grad)
		print

	return sol1




def nnls_no_prior(X, y):

	m = len(X)
	n = len(X[0])
	passive_set = []
	active_set = range(n)
	sol1 = [0]*n
	lagr = dot(X.T,(y - dot(X, sol1)))

	tol = .001
	while (len(active_set) != 0 and max(lagr[active_set]) > tol):

		var_to_add = list(lagr).index(max(lagr[active_set]))
		passive_set.append(var_to_add)
		active_set.remove(var_to_add)
		X_passive_only = X.T[passive_set].T

		pseudo_inv = dot(inv(dot(X_passive_only.T, X_passive_only)), X_passive_only.T)
		sol2 = dot(pseudo_inv, y)
		# print 'AA ' + str(sol1)
		# print 'A ' + str(sol2)

		#inner loop here
		while (min(sol2) <= 0):
			#add the zeros in, this will help with indexing
			sol2_full = [0]*n
			for i in range(len(passive_set)):
				sol2_full[passive_set[i]] = sol2[i]
			sol2_full = np.array(sol2_full)
			# print 'D ' + str(sol2_full)
			a = -1
			for i in range(len(sol2_full)):
				if i in passive_set:
					a_i = sol1[i] / (sol1[i] - sol2_full[i])
					if a_i < a or a == -1:
						a = a_i
			#not sure why only fnnls has this negative
			a = -a
			sol1 = sol1 + a*(sol2_full - sol1)
			# print 'C ' + str(sol1)
			passive_set = []
			active_set = []
			for i in range(len(sol1)):
				if sol1[i] < tol:
					active_set.append(i)
				else:
					passive_set.append(i)
			X_passive_only = X.T[passive_set].T
			pseudo_inv = dot(inv(dot(X_passive_only.T, X_passive_only)), X_passive_only.T)
			sol2 = dot(pseudo_inv, y)
			# print 'B ' + str(sol2)

		sol1 = [0]*n
		for i in range(len(passive_set)):
			sol1[passive_set[i]] = sol2[i]
		sol1 = np.array(sol1)
		lagr = dot(X.T,(y - dot(X, sol1)))


	return sol1





if __name__ == "__main__":

	x = [[1,2,3],
		[2,3,1],
		[5,8,1],
		[7,7,8],
		[9,6,2],
		[1,5,8],
		[2,4,4],
		[6,3,0],
		[2,2,7],
		[2,1,2]
	]

	y = [1,2,3,4,4,5,6,7,7,7]

	#LS
	model = LinearRegression()
	model.fit(x,y)
	print model.coef_


	#NNLS
	print nnls(x, y)[0]

	X = np.array(x)
	y = np.array(y)

	
	#My implementation of NNLS
	print nnls_no_prior(X,y)

	#My NNLS with prior

	print nnls_with_prior(X, y, np.array([.50,-.1,.30]), 10000.)



	'''
	m = len(X)
	n = len(X[0])
	passive_set = []
	active_set = range(n)
	sol1 = [0]*n
	w = dot(X.T,(y - dot(X, sol1)))

	tol = .001
	while (len(active_set) != 0 and max(w[active_set]) > tol):

		var_to_add = list(w).index(max(w[active_set]))
		passive_set.append(var_to_add)
		active_set.remove(var_to_add)
		X_passive_only = X.T[passive_set].T

		pseudo_inv = dot(inv(dot(X_passive_only.T, X_passive_only)), X_passive_only.T)
		sol2 = dot(pseudo_inv, y)
		# print 'AA ' + str(sol1)
		# print 'A ' + str(sol2)

		#inner loop here
		while (min(sol2) <= 0):
			#add the zeros in, this will help with indexing
			sol2_full = [0]*n
			for i in range(len(passive_set)):
				sol2_full[passive_set[i]] = sol2[i]
			sol2_full = np.array(sol2_full)
			# print 'D ' + str(sol2_full)
			a = -1
			for i in range(len(sol2_full)):
				if i in passive_set:
					a_i = sol1[i] / (sol1[i] - sol2_full[i])
					if a_i < a or a == -1:
						a = a_i
			#not sure why only fnnls has this negative
			a = -a
			sol1 = sol1 + a*(sol2_full - sol1)
			# print 'C ' + str(sol1)
			passive_set = []
			active_set = []
			for i in range(len(sol1)):
				if sol1[i] < tol:
					active_set.append(i)
				else:
					passive_set.append(i)
			X_passive_only = X.T[passive_set].T
			pseudo_inv = dot(inv(dot(X_passive_only.T, X_passive_only)), X_passive_only.T)
			sol2 = dot(pseudo_inv, y)
			# print 'B ' + str(sol2)

		sol1 = [0]*n
		for i in range(len(passive_set)):
			sol1[passive_set[i]] = sol2[i]
		sol1 = np.array(sol1)
		w = dot(X.T,(y - dot(X, sol1)))

	print sol1
	'''