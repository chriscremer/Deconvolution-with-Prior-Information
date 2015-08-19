
import numpy as np
import itertools

import multiprocessing as mp


def get_possible_Ws(freqs, allow_combos=True):
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


		if allow_combos:
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

	import global_variables as gv
	#samp_indexes = [i for i in range(len(possible_Ws))]
	samp_indexes = range(len(gv.X))
	numb_cpus = mp.cpu_count()
	#numb_cpus = 10
	#print 'numb of cpus' + str(numb_cpus)
	pool = mp.Pool()

	#TODO
	#need to make sure that the order it returns is the correct order
	#the get(99) allows me to exit while running
	W = pool.map_async(doWork, samp_indexes).get(99999)

	#I think this fixed the resource problem, not sure how/why
	pool.close()
	pool.join()
	#print 'W ' + str(len(W)) + ' ' + str(len(W[0]))

	#print np.array(W).shape
	#return np.reshape(np.array(W), (len(W)))

	return np.array(W)



def doWork(samp_index):
	'''
	This is the function called by select_w_parallel.
	It takes as imput the sample index that it will work on.
	'''

	#print 'Doing sample ' + str(samp_index)

	import global_variables as gv

	best_perm = []
	best_norm = -1
	for perm in gv.Ws[samp_index]:

		perm1 = np.array(perm)
		X_hat = np.dot(gv.TZ.T, perm1)
		norm = np.linalg.norm(gv.X[samp_index] - X_hat)

		if norm < best_norm or best_norm == -1:
			best_norm = norm
			best_perm = np.array(perm)

	return best_perm





def same_numb_of_entries(numb_model_subpops, freqs):
	'''
	Make all frequencies have same number of entries
	'''

	new_freqs = []
	#for each sample
	for j in range(len(freqs)):
		freq = list(freqs[j])
		#if less than number of model, add zeros
		while len(freq) < numb_model_subpops:
			freq.append(0.0)
		#if  more than number of model,then take largest
		if len(freq) > numb_model_subpops:
			freq = sorted(freq, reverse=True)
			while len(freq) > numb_model_subpops:
				freq.pop(len(freq) - 1)
		#scale so sums to 1
		freq_sum = sum(freq)
		if freq_sum != 1.0:
			for i in range(len(freq)):
				freq[i] = freq[i] / float(freq_sum)
		np.random.shuffle(freq)
		new_freqs.append(freq)

	new_freqs = np.array(new_freqs)
	#print 'new_freqs shape ' + str(new_freqs.shape)

	return new_freqs



def same_numb_of_entries_no_shuffle(numb_model_subpops, freqs):
	'''
	Make all frequencies have same number of entries, WITHOUT SHUFFLING, used for comparing at the end
	'''

	start_freqs = []
	#for each sample
	for j in range(len(freqs)):
		freq = list(freqs[j])
		#if less than number of model, add zeros
		while len(freq) < numb_model_subpops:
			freq.append(0.0)
		#if  more than number of model,then take largest
		if len(freq) > numb_model_subpops:
			freq = sorted(freq, reverse=True)
			while len(freq) > numb_model_subpops:
				freq.pop(len(freq) - 1)
		#scale so sums to 1
		freq_sum = sum(freq)
		if freq_sum != 1.0:
			for i in range(len(freq)):
				freq[i] = freq[i] / float(freq_sum)
		#np.random.shuffle(freq)
		start_freqs.append(freq)

	start_freqs = np.array(start_freqs)
	#print 'start_freqs shape ' + str(start_freqs.shape)

	return start_freqs


def match_components_to_profiles(TZ, subpops):
	'''
	Match the components to their profiles so comparing makes sense
	so for each actual profile, find the row of TZ that is most similar to it

	Need to account for when #subpops != #model components
	'''
	# possible_component_order = list(set(itertools.permutations(range(len(TZ)))))
	# best_norm_sum = -1
	# best_order = -1
	# #for each possible ordering of the components
	# for i in range(len(possible_component_order)):
	# 	#keep track of the sum of the norms
	# 	norm_sum = 0
	# 	#for each profile 
	# 	for profile_index in range(len(subpops)):
	# 		#add to the norm sum
	# 		#the norm of the difference betweem the profile and the corresponding component given this order

	# 		print 'profile index ' + str(profile_index)
	# 		print 'len subpops ' + str(len(subpops))
	# 		print 'i ' + str(i)
	# 		print 'len possible_component_order ' + str(len(possible_component_order))
	# 		print 'len possible_component_order at i' + str(len(possible_component_order[i]))

	# 		#error comes from there being more actual subpops than modeled .
	# 		#i 

	# 		norm_sum += np.linalg.norm(subpops[profile_index] - TZ[possible_component_order[i][profile_index]])

	# 	if norm_sum < best_norm_sum or best_norm_sum == -1:
	# 		best_norm_sum = norm_sum
	# 		best_order = i
	# component_order = possible_component_order[best_order]
	# return component_order


	#Will start over
	#For each model component (row of TZ), find the profile (row of subpops) that minimizes norm
	#then return a list of indexes refering to the index of the profiles 
	#so its TZ[component_order]
	#component_order is list of indexes refering to profiles

	#NO
	#its the other way around
	#for each profile, find the component that models it best

	component_order = []
	norm_sum = 0
	for i in range(len(subpops)):
		best_norm = -1
		best_index = -1

		for j in range(len(TZ)):
			norm = np.linalg.norm(subpops[i] - TZ[j])

			if norm < best_norm or best_norm == -1:
				best_norm = norm
				best_index = i

		component_order.append(best_index)
		norm_sum += best_norm

	return component_order, norm_sum





def select_w_parallel_NEW():

	import global_variables as gv

	samp_indexes = range(len(gv.X))
	numb_cpus = mp.cpu_count()

	pool = mp.Pool()
	W = pool.map_async(doWork_NEW, samp_indexes).get(99999)
	pool.close()
	pool.join()

	return np.array(W)


def doWork_NEW(samp_index):

	import global_variables as gv
	#get its freqs
	fs = gv.freqs[samp_index]
	#descending
	fs = sorted(fs, reverse=True)
	n_profiles = len(gv.TZ)
	weight_vector = np.zeros(n_profiles)
	current_v = weight_vector
	current_norm = -1
	for f in fs:

		if f == 0.0:
			break

		best_norm = -1

		for i in range(n_profiles):
			v = current_v.copy()
			v[i] = f
			v_scaled = v / sum(v)
			X_hat = np.dot(gv.TZ.T, v_scaled)
			norm = np.linalg.norm(gv.X[samp_index] - X_hat)
			if norm < best_norm or best_norm == -1:
				best_norm = norm
				best_v = v
				best_scaled_v = v_scaled

		current_v = best_v
		current_norm = best_norm
		current_scaled_v = best_scaled_v

	if best_norm == -1:
		print 'IT SHOULDNT BE -1'

	if best_norm > gv.norms[samp_index] and gv.norms[samp_index] != -1:
		return gv.W[samp_index]
	else:
		gv.norms[samp_index] = best_norm
		return current_scaled_v