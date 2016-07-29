


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

import time

from scipy.stats import multivariate_normal as mn

def main():


	#make my own data
	# data = []
	# purities = []
	# dimensions = 100
	# component_1 = np.random.rand(dimensions) 
	# component_2 = np.random.rand(dimensions)
	# component_3 = np.random.rand(dimensions)
	# component_4 = np.random.rand(dimensions)
	# hidden_profiles = [component_1, component_2, component_3, component_4]
	# for i in range(1,10,3):
	#     for j in range(1,10,3):
	#     	for k in range(1,10,3):
	#     		for l in range(1,10,3):
	# 		        norm_i = i*1.0/(i+j+k+l)
	# 		        norm_j = j*1.0/(i+j+k+l)
	# 		        norm_k = k*1.0/(i+j+k+l)
	# 		        norm_l = l*1.0/(i+j+k+l)
	# 		        # qwer = component_1*norm_i + component_2*norm_j
	# 		        # qwer = np.array(qwer)
	# 		        # print 'qwer shape ' + str(qwer.shape)
	# 		        # afds
	# 		        purities.append([norm_i, norm_j, norm_k, norm_l])
	# 		        data.append(component_1*norm_i + component_2*norm_j + component_3*norm_k + component_4*norm_l)

	# data = np.array(data)
	# print 'Data shape ' + str(data.shape)
	# X = data

	# X = np.array(X)
	# gv.set_X_global(X)
	# print X.shape

	# freqs = np.array(purities)
	# gv.set_freqs_global(freqs)

	#PARAMETERS
	# plot_file_name = '../plots/four_subpops.pdf'
	# numb_samps = 50
	# numb_feats = 10000

	k = 20
	min_components = k
	max_components = k
	# numb_of_contributing_profiles = 2
	numb_of_iterations = 1
	numb_of_iters_to_remove_local_minima = 2
	init_types = ['Random_Samples']

	#list of lists, len(means) = #types, len(means[0]) = #components
	# means = [[] for x in init_types]
	# stds = [[] for x in init_types]


	start = time.time()

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
		scores_to_avg = []
		scores_to_avg2 = []
		scores_to_avg3 = []
		for iteration in range(numb_of_iterations):

			print 'Iter ' + str(iteration)

			scores = []
			scores2 = []
			scores3 = []
			# p_h_list = [.1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0]
			# n_samples = [10, 20, 50, 100, 200, 500, 1000]
			# n_dimensions = [10, 100, 500, 1000, 5000, 10000, 20000, 30000]
			n_dimensions = [10, 100, 500, 1000, 5000, 10000]

			for n_dim in n_dimensions:

				print 'Testing ' + str(n_dim)	+ ' dimensions'	
				#Make data
				percent_hidden = 0.5
				numb_subpops = k
				numb_feats = n_dim
				numb_samps = 100
				X, freqs, real_profiles = make_convoluted_data.run_and_return(numb_subpops, numb_feats, numb_samps, percent_hidden)

				gv.set_X_global(X)
				gv.set_freqs_global(freqs)

				print X.shape
				print freqs.shape
				print 'Example freq ' + str(freqs[0])
				print real_profiles.shape


				#make all frequencies have same number of entries
				# new_freqs = wt.same_numb_of_entries(numb_components, freqs)
				#make all frequencies have same number of entries, WITHOUT SHUFFLING, used for comparing at the end
				# start_freqs = wt.same_numb_of_entries_no_shuffle(numb_components, freqs)

				# X = samps
				# gv.set_X_global(X)

				# possible_Ws = wt.get_possible_Ws(new_freqs)
				# gv.set_Ws_global(possible_Ws)

				#profile_norm_store = []
				#freqs_norm_store = []


				#for each initialization type
				for init_type in range(len(init_types)):

					# print '-------Testing ' + init_types[init_type] + '--------'

					#try the same initilization type multiple times to prevent really bad local minimums
					best_norm = -1
					for try1 in range(numb_of_iters_to_remove_local_minima):

						#Initializing model
						TZ = mt.init_model(init_types[init_type], numb_components, X)
						gv.set_current_TZ(TZ)
						print 'Optimizing model..'
						W, TZ = mt.optimize_model(TZ)

						X_hat = np.dot(W, TZ)
						norm = np.linalg.norm(X - X_hat)
						if norm < best_norm or best_norm == -1:
							best_norm = norm
							best_W = W
							best_TZ = TZ

					W = best_W
					#print W
					TZ = best_TZ
					X_hat = np.dot(W, TZ)


					sum_of_minimums = []
					sum_of_avgs = []
					sum_of_dif = []
					for learned_profile in range(len(TZ)):
						hidden_norms = []
						for hidden_profile in range(len(real_profiles)):
							hidden_norms.append(np.linalg.norm(TZ[learned_profile] - real_profiles[hidden_profile]) / len(TZ[0]))
						val, idx = min((val, idx) for (idx, val) in enumerate(hidden_norms))
						sum_of_minimums.append(val)
						sum_of_avgs.append(np.mean(hidden_norms))
						sum_of_dif.append(np.mean(hidden_norms) - val)
						# print 'Predicted k ' + str(learned_profile) + ' avg norm with hidden profiles ' + str(np.mean(hidden_norms)) + ' min norm ' + str(val) + ' with real profile' + str(idx)
					avg_minimums = np.mean(sum_of_minimums)
					avg_avgs = np.mean(sum_of_avgs)
					avg_difs = np.mean(sum_of_dif)

					scores.append(avg_minimums)
					scores2.append(avg_avgs)
					scores3.append(avg_difs)

			scores_to_avg.append(scores)
			scores_to_avg2.append(scores2)
			scores_to_avg3.append(scores3)

		#scores to avg is a list of lists. with the scores for each iteration
		#need to avg the scores 
		scores_to_plot = []
		for ii in range(len(scores_to_avg[0])):
			sum1 = 0.0
			for jj in range(len(scores_to_avg)):
				sum1 += scores_to_avg[jj][ii]
			scores_to_plot.append(sum1/len(scores_to_avg))

		scores_to_plot2 = []
		for ii in range(len(scores_to_avg2[0])):
			sum1 = 0.0
			for jj in range(len(scores_to_avg2)):
				sum1 += scores_to_avg2[jj][ii]
			scores_to_plot2.append(sum1/len(scores_to_avg2))


		scores_to_plot3 = []
		for ii in range(len(scores_to_avg3[0])):
			sum1 = 0.0
			for jj in range(len(scores_to_avg3)):
				sum1 += scores_to_avg3[jj][ii]
			scores_to_plot3.append(sum1/len(scores_to_avg3))


			

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
				# for profile in range(len(TZ)):
				# 	best_cor = 0
				# 	worst_cor = 999
				# 	best_samp = ''
				# 	worst_samp = ''
				# 	for samp in range(len(X)):
				# 		#print 'Samp ' + str(samp) + ' norm ' + str(np.linalg.norm(profile - X[samp]))
				# 		#print 'Samp ' + str(samp) + ' cor ' + str(scipy.stats.spearmanr(TZ[profile], X[samp])[0])
				# 		cor = scipy.stats.spearmanr(TZ[profile], X[samp])[0]
				# 		if cor > best_cor:
				# 			best_cor = cor
				# 			best_samp = samp
				# 		if cor < worst_cor:
				# 			worst_cor = cor
				# 			worst_samp = samp
				# 	print 'Profile ' + str(profile) + ' closest samp is ' + str(best_samp) + ' with cor of ' + str(best_cor)
				# 	print 'Profile ' + str(profile) + ' furthest samp is ' + str(worst_samp) + ' with cor of ' + str(worst_cor)

				# for profile in range(len(TZ)):
				# 	best_norm = 999
				# 	worst_norm = 0
				# 	best_samp = ''
				# 	worst_samp = ''
				# 	for samp in range(len(X)):
				# 		#print 'Samp ' + str(samp) + ' norm ' + str(np.linalg.norm(profile - X[samp]))
				# 		#print 'Samp ' + str(samp) + ' cor ' + str(scipy.stats.spearmanr(TZ[profile], X[samp])[0])
				# 		norm = np.linalg.norm(TZ[profile] - X[samp])
				# 		if norm < best_norm:
				# 			best_norm = norm
				# 			best_samp = samp
				# 		if norm > worst_norm:
				# 			worst_norm = norm
				# 			worst_samp = samp
				# 	print 'Profile ' + str(profile) + ' best norm samp is ' + str(best_samp) + ' with norm of ' + str(best_norm)
				# 	print 'Profile ' + str(profile) + ' worst norm samp is ' + str(worst_samp) + ' with norm of ' + str(worst_norm)


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
	text = '#Samples = ' + str(numb_samps) + '\n#Mixing Ratio = ' + str(percent_hidden) + '\n#Components = ' + str(k) + '\n#Iterations = ' + str(numb_of_iterations) + '\n#Random Restarts = ' + str(numb_of_iters_to_remove_local_minima) 
	# pbc.plot_line_with_text(scores_to_plot3, p_h_list, 'L2-Norm', 'Fraction Non-Zero', '../plots/norm_vs_frac_non_zero9.png', text)
	# pbc.plot_two_lines_with_text(scores_to_plot, scores_to_plot2, p_h_list, 'L2-Norm', 'Fraction Non-Zero', '../plots/norm_vs_frac_non_zero9.png', text)
	# pbc.plot_three_lines_with_text(scores_to_plot, scores_to_plot2, scores_to_plot3, n_dimensions, 'L2-Norm', 'Number of Dimensions', '../plots/norm_vs_dimensions.png', text)
	pbc.plot_three_lines_with_text_with_xlog_scale(scores_to_plot, scores_to_plot2, scores_to_plot3, n_dimensions, 'L2-Norm', 'Number of Dimensions', '../plots/norm_vs_dimensions5.png', text)


	print '\nDONE'
	end = time.time()
	print end - start






if __name__ == "__main__":

	main()



