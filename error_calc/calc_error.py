

import numpy as np

def match_W_error(W, Z, fractions, subpops):

	#match W instead for each samp
	from munkres import Munkres, print_matrix
	list_to_store_W_errors = []
	list_to_store_Z_errors = []
	for i in range(len(W)):
		# print 'W[i] ' + str(W[i])
		# print 'fractions[i] ' + str(fractions[i])
		norm_matrix = []
		for j in range(len(fractions[i])):
			dist_to_this_weight = []
			for k in range(len(W[i])):
				dist_to_this_weight.append(np.linalg.norm(W[i][k]-fractions[i][j]))

				# print str(j) + ' to ' + str(k) + ' ' + str(np.linalg.norm(W[i][j]-fractions[i][k]))
			norm_matrix.append(dist_to_this_weight)


		m = Munkres()
		indexes = m.compute(norm_matrix)
		indexes2 = [x[1] for x in indexes]

		# print indexes2

		rearranged_Wi = W[i].T[indexes2].T
		rearranged_Z = Z[indexes2]

		# print rearranged_Wi

		Wi_error = np.linalg.norm(rearranged_Wi - fractions[i])

		#for Z error, only compare Zi that are non zero in real fractions
		individual_Z_errors = []
		for j in range(len(fractions[i])):

			if fractions[i][j] > 0:

				#compare Z[correct index] and subpop[j]

				# print (np.linalg.norm(Z[j]-subpops[j]))
				# print np.linalg.norm(subpops[j])

				individual_Z_errors.append((np.linalg.norm(rearranged_Z[j]-subpops[j])) / np.linalg.norm(subpops[j]))

				# print individual_Z_errors

		Zi_error = sum(individual_Z_errors) / len(individual_Z_errors)


		list_to_store_W_errors.append(Wi_error)
		list_to_store_Z_errors.append(Zi_error)

	W_error = np.mean(list_to_store_W_errors)
	Z_error = np.mean(list_to_store_Z_errors)

	return W_error, Z_error


def match_Z_error(W, Z, fractions, subpops):

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

	rearranged_fractions = fractions.T[indexes2].T
	rearranged_profiles = subpops[indexes2]

	W_error = np.linalg.norm(rearranged_fractions - W)
	Z_error = np.linalg.norm(rearranged_profiles - Z) / np.linalg.norm(rearranged_profiles)

	return W_error, Z_error

