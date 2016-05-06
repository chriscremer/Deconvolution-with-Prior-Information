


import numpy as np
import math
from scipy.optimize import nnls

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

def scale_to_one(W):
	#scale W for each sample so that sum = 1
	for i in range(len(W)):
		W[i] = W[i]/sum(W[i])

	return W

def rearrange(subpops, Z, W):

	#Hungarian Algorithm to match predicted to real Zi
	norm_matrix = []
	for learned_profile in range(len(Z)):
		this_samp = []
		for real_profile in range(len(subpops)):
			# this_samp.append(sum(abs(Z[learned_profile] - subpops[real_profile])))
			this_samp.append(np.linalg.norm(Z[learned_profile] - subpops[real_profile]))
		norm_matrix.append(this_samp)
	from munkres import Munkres, print_matrix
	m = Munkres()
	indexes = m.compute(norm_matrix)
	indexes2 = [x[1] for x in indexes]

	rearranged_fractions = W.T[indexes2].T
	# rearranged_profiles = subpops[indexes2]

	return rearranged_fractions

def optimize_z_using_nnls(X, W):

	#Optimize Z
	Z = []
	for d in range(len(X.T)):
		Z_d = nnls(W, X.T[d])[0]
		Z.append(Z_d)
	Z = np.array(Z).T

	return Z

def optimize_w_using_nnls(X, Z):

	#Optimize W
	W = []
	for i in range(len(X)):
		W_i = nnls(Z.T, X[i])[0]
		W.append(W_i)
	W = np.array(W)

	W = use_all_components(X, W, Z)
	W = scale_to_one(W)

	return W

def optimize_w_using_L1_L1(X, Z, lmbd):


	N = len(X)
	K = len(Z)
	D = len(X[0])

	c = [0.]*K
	c.extend([1.]*D)
	c.extend([lmbd*1.]*K)
	# c = np.array(c)
	c = mtrx(c)

	G_row = np.concatenate((Z.T, -np.identity(D)), axis=1) 
	G_row = np.concatenate((G_row, np.zeros((D,K))), axis=1)
	
	G_row2 = np.concatenate((-Z.T, -np.identity(D)), axis=1) 
	G_row2 = np.concatenate((G_row2, np.zeros((D,K))), axis=1)
	G = np.concatenate((G_row,G_row2),axis=0)

	G_row3 = np.concatenate((np.identity(K), np.zeros((K,D))), axis=1) 
	G_row3 = np.concatenate((G_row3, -np.identity(K)), axis=1)
	G = np.concatenate((G,G_row3),axis=0)

	G_row4 = np.concatenate((-np.identity(K), np.zeros((K,D))), axis=1) 
	G_row4 = np.concatenate((G_row4, -np.identity(K)), axis=1)
	G = np.concatenate((G,G_row4),axis=0)

	G_row5 = np.concatenate((-np.identity(K), np.zeros((K,D))), axis=1) 
	G_row5 = np.concatenate((G_row5, np.zeros((K,K))), axis=1)
	G = np.concatenate((G,G_row5),axis=0)
	G = mtrx(G)

	#Sum to one
	Aeq = [1.]*K
	Aeq.extend([0.]*(D+K))
	Aeq = np.array(Aeq)
	Aeq = mtrx(Aeq).T
	# print Aeq.size
	beq = [1.]
	beq = mtrx(beq)

	#Optimize W
	W = []
	for i in range(len(X)):

		h = list(X[i])
		h.extend(-X[i])
		h.extend([0]*K)
		h.extend([0]*K)
		h.extend([0]*K)
		# h = np.array(h)
		h = mtrx(h)

		solvers.options['show_progress'] = False
		# print i
		sol = solvers.lp(c, G, h, A=Aeq, b=beq)
		# print str(i) + ' done'
		W_i = np.array([x for x in sol['x'][:K]])

		W.append(W_i)
	W = np.array(W)

	W = use_all_components(X, W, Z)
	# W = scale_to_one(W)

	return W

def optimize_z_using_L1(X, W):


	N = len(X)
	K = len(W[0])
	D = len(X[0])

	c = [0.]*K
	c.extend([1.]*N)
	c = mtrx(c)

	G_row = np.concatenate((W, -np.identity(N)), axis=1) 	
	G_row2 = np.concatenate((-W, -np.identity(N)), axis=1) 
	G = np.concatenate((G_row,G_row2),axis=0)
	G_row3 = np.concatenate((-np.identity(K), np.zeros((K,N))), axis=1) 
	G = np.concatenate((G,G_row3),axis=0)
	G = mtrx(G)

	#Optimize W
	Z_T = []
	for i in range(len(X.T)):

		h = list(X.T[i])
		h.extend(-X.T[i])
		h.extend([0]*K)
		h = mtrx(h)

		solvers.options['show_progress'] = False
		# print i
		sol = solvers.lp(c, G, h)
		# print str(i) + ' done'

		# print len(sol['x'])


		Z_i = np.array([x for x in sol['x'][:K]])

		# print len(Z_i)

		Z_T.append(Z_i)
	Z = np.array(Z_T).T

	# W = use_all_components(X, W, Z)
	# W = scale_to_one(W)

	return Z

def optimize_w_using_L2_L1(X, Z, lmbd):

	N = len(X)
	K = len(Z)
	D = len(X[0])

	P_row0 = np.concatenate((2*np.dot(Z,Z.T), np.zeros((K,K))), axis=1) 
	P_row1 = np.concatenate((np.zeros((K,K)), np.zeros((K,K))), axis=1) 
	P = np.concatenate((P_row0,P_row1),axis=0)
	P = mtrx(P)

	G_row1 = np.concatenate((np.identity(K), -np.identity(K)), axis=1) 
	G_row2 = np.concatenate((np.identity(K), -np.identity(K)), axis=1) 
	G = np.concatenate((G_row1,G_row2),axis=0)

	G_row3 = np.concatenate((-np.identity(K), np.zeros((K,K))), axis=1) 
	G = np.concatenate((G,G_row3),axis=0)
	G = mtrx(G)

	h = [0.]*(3*K)
	h = mtrx(h)

	Aeq = [1.]*K
	Aeq.extend([0.]*K)
	Aeq = np.array(Aeq)
	Aeq = mtrx(Aeq).T
	# print Aeq.size
	beq = [1.]
	beq = mtrx(beq)

	#Optimize W
	W = []
	for i in range(len(X)):

		q = np.concatenate((-2*np.dot(Z,X[i]), lmbd*np.ones(K)), axis=0)
		q = mtrx(q)

		solvers.options['show_progress'] = False
		# print i
		sol = solvers.qp(P, q, G, h, A=Aeq, b=beq)
		# print str(i) + ' done'
		W_i = np.array([x for x in sol['x'][:K]])

		W.append(W_i)
	W = np.array(W)

	# W = use_all_components(X, W, Z)
	# W = scale_to_one(W)

	return W

def optimize_w_using_L1_L2(X, Z, lmbd):

	N = len(X)
	K = len(Z)
	D = len(X[0])

	P = 2*np.identity(K+D)
	P = mtrx(P)

	q = [0.]*K
	q.extend([lmbd*1.]*D)
	q = mtrx(q)

	G_row1 = np.concatenate((Z.T, -np.identity(D)), axis=1) 
	G_row2 = np.concatenate((-Z.T, -np.identity(D)), axis=1) 
	G = np.concatenate((G_row1,G_row2),axis=0)

	G_row3 = np.concatenate((-np.identity(K), np.zeros((K,D))), axis=1) 
	G = np.concatenate((G,G_row3),axis=0)
	G = mtrx(G)

	Aeq = [1.]*K
	Aeq.extend([0.]*D)
	Aeq = np.array(Aeq)
	Aeq = mtrx(Aeq).T
	# print Aeq.size
	beq = [1.]
	beq = mtrx(beq)

	#Optimize W
	W = []
	for i in range(len(X)):

		h = list(X[i])
		h.extend(-X[i])
		h.extend([0]*K)
		h = mtrx(h)

		solvers.options['show_progress'] = False
		# print i
		sol = solvers.qp(P, q, G, h, A=Aeq, b=beq)
		# print str(i) + ' done'
		W_i = np.array([x for x in sol['x'][:K]])

		W.append(W_i)
	W = np.array(W)

	W = use_all_components(X, W, Z)
	# W = scale_to_one(W)

	return W

def optimize_using_L2_L1(A, B, lmbd):

	#Solves Xi = argmin (A*Xi - Bi)^2 + |x| for i...N
	#since B = AX

	N = len(B)
	D = len(B[0])
	K = len(A.T)

	P_row0 = np.concatenate((2*np.dot(A.T,A), np.zeros((K,K))), axis=1) 
	P_row1 = np.concatenate((np.zeros((K,K)), np.zeros((K,K))), axis=1) 
	P = np.concatenate((P_row0,P_row1),axis=0)
	P = mtrx(P)

	G_row1 = np.concatenate((np.identity(K), -np.identity(K)), axis=1) 
	G_row2 = np.concatenate((-np.identity(K), -np.identity(K)), axis=1) 
	G = np.concatenate((G_row1,G_row2),axis=0)

	G_row3 = np.concatenate((-np.identity(K), np.zeros((K,K))), axis=1) 
	G = np.concatenate((G,G_row3),axis=0)
	G = mtrx(G)

	h = [0.]*(3*K)
	h = mtrx(h)

	# Aeq = [1.]*K
	# Aeq.extend([0.]*K)
	# Aeq = np.array(Aeq)
	# Aeq = mtrx(Aeq).T

	# # print Aeq.size


	# beq = [1.]
	# beq = mtrx(beq)

	# print beq.size

	#Optimize X
	X = []
	for i in range(len(B)):

		q = np.concatenate((-2*np.dot(A.T,B[i]), lmbd*np.ones(K)), axis=0)
		q = mtrx(q)

		solvers.options['show_progress'] = False
		# print i
		sol = solvers.qp(P, q, G, h)
		# print str(i) + ' done'
		X_i = np.array([x for x in sol['x'][:K]])

		X.append(X_i)
	X = np.array(X)


	return X

def L1_norm(M):

	sum1 = 0
	for i in range(len(M)):
		for j in range(len(M.T)):
			sum1 += abs(M[i][j])

	return sum1

def L1_norm_vector(v):

	sum1 = sum(abs(v))

	return sum1

def model_1(X, n_comps, tol, real_weights, subpops):

	#Init Z
	Z = init_Z(X, n_comps)

	last_error = -1
	current_error = -1

	#While error decreases
	while (last_error - current_error) > tol or last_error == -1:

		#Solve for W
		W = optimize_w_using_nnls(X, Z)

		#debug
		# rearranged_fractions = rearrange(subpops, Z, W)
		# print rearranged_fractions[0]
		# print real_weights[0]

		#Solve for Z
		Z = optimize_z_using_nnls(X, W)

		last_error = current_error
		current_error = np.linalg.norm(X - np.dot(W,Z))
		print current_error

	return W, Z

def model_L1_L1(X, n_comps, tol, real_weights, subpops, lmbd):

	#Init Z
	Z = init_Z(X, n_comps)

	last_error = -1
	current_error = -1

	#While error decreases
	while (last_error - current_error) > tol or last_error == -1:

		#Solve for W
		W = optimize_w_using_L1_L1(X, Z, lmbd)
		# print W[0]
		# print 'vs'
		# W = optimize_w_using_L1_L1(X, Z, 100000)
		# print W[0]
		# sdfa



		#debug
		# rearranged_fractions = rearrange(subpops, Z, W)
		# print rearranged_fractions[0]
		# print real_weights[0]

		# current_error2 = np.linalg.norm(X - np.dot(W,Z), 1) + lmbd*np.linalg.norm(W, 1)
		current_error2 = L1_norm(X - np.dot(W,Z)) + lmbd*L1_norm(W)
		# print 'current_error2 ' + str(current_error2)
		print 'current_error2 ' + str(current_error2)

		#Solve for Z

		# print 'X.T'
		# for i in range(len(X.T)):
		# 	print X.T[i]

		# print 'prev Z error'
		# for i in range(len(X.T)):
		# 	print Z.T[i]
		# 	print str(i) + ' ' + str(L1_norm_vector(X.T[i] - np.dot(W, Z.T[i])))

		Z = optimize_z_using_L1(X, W)

		# print 'new Z error'
		# for i in range(len(X.T)):
		# 	print Z.T[i]
		# 	print str(i) + ' ' + str(L1_norm_vector(X.T[i] - np.dot(W, Z.T[i])))

		last_error = current_error
		# current_error = sum(sum(abs(X - np.dot(W,Z))))
		# current_error = np.linalg.norm(X - np.dot(W,Z), 1) + lmbd*np.linalg.norm(W, 1)
		current_error = L1_norm(X - np.dot(W,Z)) + lmbd*L1_norm(W)
		# print 'current_error ' + str(current_error)
		print 'current_error ' + str(current_error)

	return W, Z

def model_L2_L1(X, n_comps, tol, real_weights, subpops, lmbd):

	#Init Z
	Z = init_Z(X, n_comps)

	last_error = -1
	current_error = -1

	#While error decreases
	while (last_error - current_error) > tol or last_error == -1:

		#Solve for W
		W = optimize_w_using_L2_L1(X, Z, lmbd)
		# W = optimize_using_L2_L1(Z.T, X, 100)
		W = use_all_components(X, W, Z)
		

		#debug
		# rearranged_fractions = rearrange(subpops, Z, W)
		# print rearranged_fractions[0]
		# print real_weights[0]

		#Solve for Z
		# Z = optimize_z_using_nnls(X, W)
		Z = optimize_using_L2_L1(W, X.T, 0).T


		last_error = current_error
		# current_error = sum(sum(abs(X - np.dot(W,Z))))
		current_error = np.linalg.norm(X - np.dot(W,Z)) + lmbd*L1_norm(W)
		print current_error

		# M = X - np.dot(W,Z)
		# sum1 = 0
		# for i in range(len(M)):
		# 	for j in range(len(M.T)):
		# 		sum1 += (M[i][j])**2
		# print sum1
		# print math.sqrt(sum1)



	# W = scale_to_one(W)

	return W, Z

def model_L1_L2(X, n_comps, tol, real_weights, subpops, lmbd):

	#Init Z
	Z = init_Z(X, n_comps)

	last_error = -1
	current_error = -1

	#While error decreases
	while (last_error - current_error) > tol or last_error == -1:

		#Solve for W
		W = optimize_w_using_L1_L2(X, Z, lmbd)

		#debug
		# rearranged_fractions = rearrange(subpops, Z, W)
		# print rearranged_fractions[0]
		# print real_weights[0]

		#Solve for Z
		Z = optimize_z_using_L1(X, W)

		last_error = current_error
		# current_error = sum(sum(abs(X - np.dot(W,Z))))
		current_error = L1_norm(X - np.dot(W,Z)) + lmbd*np.linalg.norm(W)

		print current_error

	return W, Z


