



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



