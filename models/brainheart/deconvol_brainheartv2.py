


#add  directory to path to get my packages there
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
#parentdir = os.path.dirname(currentdir)
#print sys.path
#sys.path.insert(0,parentdir) 
sys.path.insert(0,currentdir+'/simulated_data')
import make_convoluted_data

import numpy as np
import scipy
import csv

import itertools

import plot_bar_chart as pbc
import read_TCGA_gene_data as rtgd
import weight_tools as wt
import model_tools as mt

import global_variables as gv



from scipy.stats import multivariate_normal as mn

def main():

	data_file = '../affymetrix_data2/data/data.txt'

	X = np.genfromtxt(data_file, delimiter='\t', skip_header=1, usecols=range(1,34))
	X = X.T
	# gv.set_X_global(X)

	print X.shape

	with open(data_file, 'rU') as f:
		reader = csv.reader(f, delimiter='\t')
		for row in reader:
			header = row
			break

	header.pop(0)
	#print header


	pure_heart = X[0]
	pure_brain = X[-1]

	#make varying purity samples
	X = []
	purities = [[0.0, 1.0],
				[0.1, 0.9],
				[0.2, 0.8],
				[0.3, 0.7],
				[0.4, 0.6],
				[0.5, 0.5],
				[0.6, 0.4],
				[0.7, 0.3],
				[0.8, 0.2],
				[0.9, 0.1],
				[1.0, 0.0]]
	for purity in purities:
		X.append(pure_heart*purity[0] + pure_brain*purity[1])

	X = np.array(X)
	gv.set_X_global(X)
	print X.shape

	freqs = np.array(purities)


	# cellularities = []
	# hearts = np.array([x for x in X[0:3]])
	# brains = np.array([x for x in X[-4:-1]])

	# print 'hearts ' + str(hearts.shape)
	# print 'brains ' + str(brains.shape)

	# for samp in header:
	# 	if 'mix1' in samp:
	# 		cellularities.append([0.0, 1.0])
	# 	elif 'mix2' in samp:
	# 		cellularities.append([0.05, 0.95])
	# 	elif 'mix3' in samp:
	# 		cellularities.append([0.10, 0.90])
	# 	elif 'mix4' in samp:
	# 		cellularities.append([0.25, 0.75])
	# 	elif 'mix5' in samp:
	# 		cellularities.append([0.5, 0.5])
	# 	elif 'mix6' in samp:
	# 		cellularities.append([0.75, 0.25])
	# 	elif 'mix7' in samp:
	# 		cellularities.append([0.9, 0.1])
	# 	elif 'mix8' in samp:
	# 		cellularities.append([0.95, 0.05])
	# 	elif 'mix9' in samp:
	# 		cellularities.append([1.0, 0.0])


	# cellularities = np.array(cellularities)
	# freqs = cellularities
	# print cellularities

	#PARAMETERS
	# plot_file_name = '../plots/four_subpops.pdf'
	# numb_samps = 50
	# numb_feats = 10000
	min_components = 2
	max_components = 2
	# numb_of_contributing_profiles = 2
	numb_of_iterations = 1
	numb_of_iters_to_remove_local_minima = 2
	init_types = ['Random_Samples']

	#list of lists, len(means) = #types, len(means[0]) = #components
	# means = [[] for x in init_types]
	# stds = [[] for x in init_types]

	for numb_components in range(min_components, max_components+1):

		print '\nNumb of Components =' + str(numb_components)
		#this is the number of profiles that exist in the simulated data
		#if -1 then set them equal
		# if numb_of_contributing_profiles == -1:
		# 	numb_subpops = numb_components
		# else:
		# 	numb_subpops = numb_of_contributing_profiles

		#store the scores for each iter for each type
		# scores = [[] for x in init_types]

		#iterate x times and get average
		for iteration in range(numb_of_iterations):

			print 'Iter ' + str(iteration)

			#Make data
			# samps, freqs, subpops = make_convoluted_data.run_and_return(numb_subpops, numb_feats, numb_samps)

			#make all frequencies have same number of entries
			new_freqs = wt.same_numb_of_entries(numb_components, freqs)
			#make all frequencies have same number of entries, WITHOUT SHUFFLING, used for comparing at the end
			# start_freqs = wt.same_numb_of_entries_no_shuffle(numb_components, freqs)

			# X = samps
			# gv.set_X_global(X)

			possible_Ws = wt.get_possible_Ws(new_freqs)
			gv.set_Ws_global(possible_Ws)

			#profile_norm_store = []
			#freqs_norm_store = []


			#for each initialization type
			for init_type in range(len(init_types)):

				print '-------Testing ' + init_types[init_type] + '--------'

				#try the same initilization type multiple times to prevent really bad local minimums
				best_norm = -1
				for try1 in range(numb_of_iters_to_remove_local_minima):

					#Initializing model
					TZ = mt.init_model(init_types[init_type], numb_components, X)
					gv.set_current_TZ(TZ)
					#print 'Optimizing model..'
					W, TZ = mt.optimize_model(possible_Ws, TZ)

					X_hat = np.dot(W, TZ)
					norm = np.linalg.norm(X - X_hat)
					if norm < best_norm or best_norm == -1:
						best_norm = norm
						best_W = W
						best_TZ = TZ

				W = best_W
				print W
				TZ = best_TZ
				X_hat = np.dot(W, TZ)


				#compare profiles to real pure samples
				# norm_order1 = 0
				# for i in range(len(brains)):
				# 	norm_order1 += np.linalg.norm(TZ[0] - brains[i])
				# 	norm_order1 += np.linalg.norm(TZ[1] - hearts[i])
				# norm_order2 = 0
				# for i in range(len(brains)):
				# 	norm_order2 += np.linalg.norm(TZ[1] - brains[i])
				# 	norm_order2 += np.linalg.norm(TZ[0] - hearts[i])

				# print 'norm_order1 ' + str(norm_order1)
				# print 'norm_order2 ' + str(norm_order2)


				#for each profile, print norm with each sample
				for profile in range(len(TZ)):
					best_cor = 0
					worst_cor = 999
					best_samp = ''
					worst_samp = ''
					for samp in range(len(X)):
						#print 'Samp ' + str(samp) + ' norm ' + str(np.linalg.norm(profile - X[samp]))
						#print 'Samp ' + str(samp) + ' cor ' + str(scipy.stats.spearmanr(TZ[profile], X[samp])[0])
						cor = scipy.stats.spearmanr(TZ[profile], X[samp])[0]
						if cor > best_cor:
							best_cor = cor
							best_samp = samp
						if cor < worst_cor:
							worst_cor = cor
							worst_samp = samp
					print 'Profile ' + str(profile) + ' closest samp is ' + str(best_samp) + ' with cor of ' + str(best_cor)
					print 'Profile ' + str(profile) + ' furthest samp is ' + str(worst_samp) + ' with cor of ' + str(worst_cor)

				for profile in range(len(TZ)):
					best_norm = 999
					worst_norm = 0
					best_samp = ''
					worst_samp = ''
					for samp in range(len(X)):
						#print 'Samp ' + str(samp) + ' norm ' + str(np.linalg.norm(profile - X[samp]))
						#print 'Samp ' + str(samp) + ' cor ' + str(scipy.stats.spearmanr(TZ[profile], X[samp])[0])
						norm = np.linalg.norm(TZ[profile] - X[samp])
						if norm < best_norm:
							best_norm = norm
							best_samp = samp
						if norm > worst_norm:
							worst_norm = norm
							worst_samp = samp
					print 'Profile ' + str(profile) + ' best norm samp is ' + str(best_samp) + ' with norm of ' + str(best_norm)
					print 'Profile ' + str(profile) + ' worst norm samp is ' + str(worst_samp) + ' with norm of ' + str(worst_norm)


				# for profile in range(len(TZ)):
				# 	heart_sum = 0
				# 	for samp in range(len(hearts)):
				# 		heart_sum += scipy.stats.spearmanr(TZ[profile], hearts[samp])[0]
				# 	brain_sum = 0
				# 	for samp in range(len(brains)):
				# 		brain_sum += scipy.stats.spearmanr(TZ[profile], brains[samp])[0]
				# 	print 'Profile ' + str(profile) + ' average brain cor is ' + str(brain_sum/3.0)
				# 	print 'Profile ' + str(profile) + ' average heart cor is ' + str(heart_sum/3.0)



				# for samp in range(len(header)):
				# 	if 'mix1' in header[samp]:
				# 		print 'Samp ' + str(samp) + ' = mix1 (0, 1)'
				# 		# cellularities.append([0.0, 1.0])
				# 	elif 'mix2' in header[samp]:
				# 		print 'Samp ' + str(samp) + ' = mix2 (0.05, 0.95)'
				# 		# cellularities.append([0.05, 0.95])
				# 	elif 'mix3' in header[samp]:
				# 		print 'Samp ' + str(samp) + ' = mix3 (0.10, 0.90)'
				# 		# cellularities.append([0.10, 0.90])
				# 	elif 'mix4' in header[samp]:
				# 		print 'Samp ' + str(samp) + ' = mix4 (0.25, 0.75)'
				# 		# cellularities.append([0.25, 0.75])
				# 	elif 'mix5' in header[samp]:
				# 		print 'Samp ' + str(samp) + ' = mix5 (0.5, 0.5)'
				# 		# cellularities.append([0.5, 0.5])
				# 	elif 'mix6' in header[samp]:
				# 		print 'Samp ' + str(samp) + ' = mix6 (0.75, 0.25)'
				# 		# cellularities.append([0.75, 0.25])
				# 	elif 'mix7' in header[samp]:
				# 		print 'Samp ' + str(samp) + ' = mix7 (0.9, 0.1)'
				# 		# cellularities.append([0.9, 0.1])
				# 	elif 'mix8' in header[samp]:
				# 		print 'Samp ' + str(samp) + ' = mix8 (0.95, 0.05)'
				# 		# cellularities.append([0.95, 0.05])
				# 	elif 'mix9' in header[samp]:
				# 		print 'Samp ' + str(samp) + ' = mix9 (1.0, 0.0)'
						# cellularities.append([1.0, 0.0])

				#Likelihood = P(D|W,Z) = \prod P(d|W,Z) assuming IID
				# P(d|W,Z) = N(d|WZ, cov)

				#cov = numpy.zeros((len(TZ),len(TZ)))
				# cov = np.identity(len(X[0]))

				# L = 1
				# for samp in range(len(X)):
				# 	this_samp_L = mn.pdf(X[samp], mean=X_hat[samp], cov=cov)

				# 	print this_samp_L

				# 	aadfs

				# 	L = L * this_samp_L



				# gaussian_pdf = mn.pdf(x, mean=2.5, cov=cov)



				######################################################### 
				# match the components to their profiles so comparing makes sense
				# so for each actual profile, find the row of TZ that is most similar to it
				# component_order, norm_sum = wt.match_components_to_profiles(TZ, subpops)
				# print 'Sum profile norm ' + str(norm_sum)

				# scores[init_type].append(norm_sum)

				#profile_norm_store.append(norm_sum)

				#########################################################
				#print how close (norm) each component is to each matching profile
				# sum1=0
				# for profile_index in range(len(subpops)):
				# 	norm1 = np.linalg.norm(subpops[profile_index] - TZ[component_order[profile_index]])
				# 	#print 'Profile ' + str(profile_index) + ' norm ' + str(norm1)
				# 	sum1 += norm1
				# print 'Sum profile norm ' + str(sum1)
				# profile_norm_store.append(sum1)
				#########################################################

				#print 

				#########################################################
				#this prints the actual frequencies of samples and predicted assignment of frequencies
				# print 'Actual - Predicted'
				# for i in range(len(freqs)):
				# 	#only print first 10
				# 	if i > 9:
				# 		break
				# 	print str(['%.2f' % elem for elem in start_freqs[i]]) + '  ' + str(['%.2f' % elem for elem in [W[i][x] for x in component_order]])
				#########################################################


				#########################################################
				#print average norm of assigned frequencies vs actual
				# sum1 = 0
				# for i in range(len(freqs)):
				# 	dif_array = start_freqs[i] - [W[i][x] for x in component_order]
				# 	sum1 += np.linalg.norm(dif_array)
				# print 'Sum freq assignemtn norm ' + str(sum1)
				# freqs_norm_store.append(sum1)
				#########################################################


		#Given type and numb of components, this the avg performance
		# mean_error = np.mean(profile_norm_store)
		# mean_std = np.std(profile_norm_store)
		# #mean_freq_dif = np.mean(freqs_norm_store)

		# means[init_type].append(mean_error)
		# stds[init_type].append(mean_std)

		#take averages for each iter and type and put into lists
		# for t in range(len(init_types)):
		# 	mean = np.mean(scores[t])
		# 	std = np.std(scores[t])
		# 	means[t].append(mean)
		# 	stds[t].append(std)


	# pbc.plot_bar_chart(init_types, means, stds, [str(x) for x in range(min_components, max_components+1)], plot_file_name)

	print '\nDONE'






if __name__ == "__main__":

	main()



