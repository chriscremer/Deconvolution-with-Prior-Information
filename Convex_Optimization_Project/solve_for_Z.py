

import numpy as np
import math
from scipy.optimize import nnls

from cvxopt import matrix as mtrx
from cvxopt import solvers

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




def optimize_z_using_nnls(X, W):

	#Optimize Z
	Z = []
	for d in range(len(X.T)):
		Z_d = nnls(W, X.T[d])[0]
		Z.append(Z_d)
	Z = np.array(Z).T

	return Z