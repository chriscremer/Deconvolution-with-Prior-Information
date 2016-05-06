

import numpy as np

def model(X, Z, max_iters, tol):


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