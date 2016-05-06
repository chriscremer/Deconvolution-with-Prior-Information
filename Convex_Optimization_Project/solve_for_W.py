
import numpy as np
import math
from scipy.optimize import nnls

from cvxopt import matrix as mtrx
from cvxopt import solvers


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

	# #Sum to one
	# Aeq = [1.]*K
	# Aeq.extend([0.]*(D+K))
	# Aeq = np.array(Aeq)
	# Aeq = mtrx(Aeq).T
	# # print Aeq.size
	# beq = [1.]
	# beq = mtrx(beq)

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
		# sol = solvers.lp(c, G, h, A=Aeq, b=beq)
		sol = solvers.lp(c, G, h)
		# print str(i) + ' done'
		W_i = np.array([x for x in sol['x'][:K]])

		W.append(W_i)
	W = np.array(W)

	# W = use_all_components(X, W, Z)
	# W = scale_to_one(W)

	return W



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

	# Aeq = [1.]*K
	# Aeq.extend([0.]*K)
	# Aeq = np.array(Aeq)
	# Aeq = mtrx(Aeq).T
	# # print Aeq.size
	# beq = [1.]
	# beq = mtrx(beq)

	#Optimize W
	W = []
	for i in range(len(X)):

		q = np.concatenate((-2*np.dot(Z,X[i]), lmbd*np.ones(K)), axis=0)
		q = mtrx(q)

		solvers.options['show_progress'] = False
		# print i
		sol = solvers.qp(P, q, G, h)
		# sol = solvers.qp(P, q, G, h, A=Aeq, b=beq)
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

	# P = 2*np.identity(K+D)
	# print P.shape
	# P = mtrx(P)


	P_row0 = np.concatenate((2*lmbd*np.identity(K), np.zeros((K,D))), axis=1) 
	P_row1 = np.concatenate((np.zeros((D,K)), np.zeros((D,D))), axis=1) 
	P = np.concatenate((P_row0,P_row1),axis=0)
	P = mtrx(P)

	q = [0.]*K
	q.extend([1.]*D)
	q = mtrx(q)

	G_row1 = np.concatenate((Z.T, -np.identity(D)), axis=1) 
	G_row2 = np.concatenate((-Z.T, -np.identity(D)), axis=1) 
	G = np.concatenate((G_row1,G_row2),axis=0)

	G_row3 = np.concatenate((-np.identity(K), np.zeros((K,D))), axis=1) 
	G = np.concatenate((G,G_row3),axis=0)
	G = mtrx(G)

	# Aeq = [1.]*K
	# Aeq.extend([0.]*D)
	# Aeq = np.array(Aeq)
	# Aeq = mtrx(Aeq).T
	# # print Aeq.size
	# beq = [1.]
	# beq = mtrx(beq)

	#Optimize W
	W = []
	for i in range(len(X)):

		h = list(X[i])
		h.extend(-X[i])
		h.extend([0]*K)
		h = mtrx(h)

		solvers.options['show_progress'] = False
		# print i
		# sol = solvers.qp(P, q, G, h, A=Aeq, b=beq)
		sol = solvers.qp(P, q, G, h)

		# print str(i) + ' done'
		W_i = np.array([x for x in sol['x'][:K]])

		W.append(W_i)
	W = np.array(W)

	# W = use_all_components(X, W, Z)
	# W = scale_to_one(W)

	return W




def optimize_w_using_L2_L2(X, Z, lmbd):

	N = len(X)
	K = len(Z)
	D = len(X[0])

	P = 2*(np.dot(Z,Z.T) + lmbd*np.identity(K))
	P = mtrx(P)

	G = -np.identity(K)
	G = mtrx(G)

	h = [0.]*(K)
	h = mtrx(h)

	#Optimize W
	W = []
	for i in range(len(X)):

		q = -2*np.dot(Z,X[i])
		q = mtrx(q)

		solvers.options['show_progress'] = False
		# print i
		sol = solvers.qp(P, q, G, h)
		# sol = solvers.qp(P, q, G, h, A=Aeq, b=beq)
		# print str(i) + ' done'
		W_i = np.array([x for x in sol['x'][:K]])

		W.append(W_i)
	W = np.array(W)

	# W = use_all_components(X, W, Z)
	# W = scale_to_one(W)

	return W


def optimize_w_using_L2_sum_to_1(X, Z):

	N = len(X)
	K = len(Z)
	D = len(X[0])

	P = 2*(np.dot(Z,Z.T))
	P = mtrx(P)

	G = -np.identity(K)
	G = mtrx(G)

	h = [0.]*(K)
	h = mtrx(h)

	Aeq = [1.]*K
	Aeq = mtrx(Aeq).T
	beq = [1.]
	beq = mtrx(beq)

	#Optimize W
	W = []
	for i in range(len(X)):

		q = -2*np.dot(Z,X[i])
		q = mtrx(q)

		solvers.options['show_progress'] = False
		# print i
		# sol = solvers.qp(P, q, G, h)
		sol = solvers.qp(P, q, G, h, A=Aeq, b=beq)
		# print str(i) + ' done'
		W_i = np.array([x for x in sol['x'][:K]])

		W.append(W_i)
	W = np.array(W)

	# W = use_all_components(X, W, Z)
	# W = scale_to_one(W)

	return W

