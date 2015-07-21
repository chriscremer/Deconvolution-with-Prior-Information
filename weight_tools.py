
import numpy as np
import itertools


def get_possible_Ws(freqs):
	'''
	Given the frequencies of the subpopulations within a sample, this
	function wil return a list of lists of the possible assignments 
	of the frequencies/weights to each subpopulation. 

	TODO: I only look at combining 2 frequencies. Should allow more.
	'''

	possible_Ws = []

	for i in range(len(freqs)):

		#find permuations and remove the duplicates
		perms = list(set(itertools.permutations(freqs[i])))

		#now find all combinations of sums of freqs (combine 2)
		combs = list(itertools.combinations(range(len(freqs[i])), 2))
		for j in range(len(combs)):
			#see if this combo has any zeros
			#if it does skip it
			skip = 0
			for index in combs[j]:
				if freqs[i][index] == 0.0:
					skip =1 
					break
			if skip == 1:
				continue
			#combine the freqencies and replace one with a zero
			new_freq = list(freqs[i])
			new_freq[combs[j][0]] = new_freq[combs[j][0]] + new_freq[combs[j][1]]
			new_freq[combs[j][1]] = 0.0
			new_perms = list(set(itertools.permutations(new_freq)))
			perms.extend(new_perms)
			
		#now find all combinations of sums of freqs (combine 3)
		combs = list(itertools.combinations(range(len(freqs[i])), 3))
		for j in range(len(combs)):
			#see if this combo has any zeros
			#if it does skip it
			skip = 0
			for index in combs[j]:
				if freqs[i][index] == 0.0:
					skip =1 
					break
			if skip == 1:
				continue
			#combine the freqencies and replace one with a zero
			new_freq = list(freqs[i])
			new_freq[combs[j][0]] = new_freq[combs[j][0]] + new_freq[combs[j][1]] + new_freq[combs[j][2]]
			new_freq[combs[j][1]] = 0.0
			new_freq[combs[j][2]] = 0.0
			new_perms = list(set(itertools.permutations(new_freq)))
			perms.extend(new_perms)

		#print len(combs)
		#print combs

		perms = list(set(perms))

		possible_Ws.append(perms)

	return possible_Ws




def select_w(X, possible_Ws, TZ):
	'''
	Given the profile matrix of the subpopulations, find the best assignment
	of frequencies to subpopulations. This is done by minimizing the norm.
	'''

	W = []

	for i in range(len(possible_Ws)):

		best_perm = []
		best_norm = -1

		for perm in possible_Ws[i]:

			perm1 = np.array(perm)
			X_hat = np.dot(TZ.T, perm1)
			norm = np.linalg.norm(X[i] - X_hat)

			if norm < best_norm or best_norm == -1:
				best_norm = norm
				best_perm = np.array(perm)

		W.append(best_perm)

	return np.reshape(np.array(W), (len(W),len(W[0]))) 




def select_w_parallel():
	'''
	Same as select_w but this function used the multiple cores
	'''

	try:

		#samp_indexes = [i for i in range(len(possible_Ws))]
		samp_indexes = range(len(possible_Ws))
		numb_cpus = mp.cpu_count()
		#print 'numb of cpus' + str(numb_cpus)
		pool = mp.Pool()

		W = pool.map(doWork, samp_indexes)
		#print 'W ' + str(len(W)) + ' ' + str(len(W[0]))

	except KeyboardInterrupt:
		print 'keyboard interruption'

	#print np.array(W).shape
	#return np.reshape(np.array(W), (len(W)))

	return np.array(W)



def doWork(samp_index):
	'''
	This is the function called by select_w_parallel.
	It takes as imput the sample index that it will work on.
	'''

	#print 'Doing sample ' + str(samp_index)

	best_perm = []
	best_norm = -1
	for perm in possible_Ws[samp_index]:

		perm1 = np.array(perm)
		X_hat = np.dot(TZ.T, perm1)
		norm = np.linalg.norm(X[samp_index] - X_hat)

		if norm < best_norm or best_norm == -1:
			best_norm = norm
			best_perm = np.array(perm)

	return best_perm