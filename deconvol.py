


#add  directory to path to get my packages there
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
#parentdir = os.path.dirname(currentdir)
#print sys.path
#sys.path.insert(0,parentdir) 
sys.path.insert(0,currentdir+'/simulated_data')
import make_convoluted_data

import numpy as np

import itertools

import plot_bar_chart as pbc
import read_TCGA_gene_data as rtgd
import weight_tools as wt
import model_tools as mt

import global_variables as gv

def main():

	#PARAMETERS
	plot_file_name = '../plots/dif_inits_best_init.pdf'
	numb_samps = 50
	numb_feats = 10000
	min_components = 2
	max_components = 5
	numb_of_contributing_profiles = -1
	numb_of_iterations = 10
	numb_of_iters_to_remove_local_minima = 5
	init_types = ['Random_Values', 'Random_Samples']

	#list of lists, len(means) = #types, len(means[0]) = #components
	means = [[] for x in init_types]
	stds = [[] for x in init_types]

	for numb_components in range(min_components, max_components+1):

		print '\nNumb of Components =' + str(numb_components)
		#this is the number of profiles that exist in the simulated data
		#if -1 then set them equal
		if numb_of_contributing_profiles == -1:
			numb_subpops = numb_components
		else:
			numb_subpops = numb_of_contributing_profiles

		#store the scores for each iter for each type
		scores = [[] for x in init_types]

		#iterate x times and get average
		for iteration in range(numb_of_iterations):

			print 'Iter ' + str(iteration)

			#Make data
			samps, freqs, subpops = make_convoluted_data.run_and_return(numb_subpops, numb_feats, numb_samps)

			#make all frequencies have same number of entries
			new_freqs = wt.same_numb_of_entries(numb_components, freqs)
			#make all frequencies have same number of entries, WITHOUT SHUFFLING, used for comparing at the end
			start_freqs = wt.same_numb_of_entries_no_shuffle(numb_components, freqs)

			X = samps
			gv.set_X_global(X)

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
				TZ = best_TZ

				######################################################### 
				# match the components to their profiles so comparing makes sense
				# so for each actual profile, find the row of TZ that is most similar to it
				component_order, norm_sum = wt.match_components_to_profiles(TZ, subpops)
				print 'Sum profile norm ' + str(norm_sum)

				scores[init_type].append(norm_sum)

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
		for t in range(len(init_types)):
			mean = np.mean(scores[t])
			std = np.std(scores[t])
			means[t].append(mean)
			stds[t].append(std)


	pbc.plot_bar_chart(init_types, means, stds, [str(x) for x in range(min_components, max_components+1)], plot_file_name)

	print '\nDONE'






if __name__ == "__main__":

	main()



