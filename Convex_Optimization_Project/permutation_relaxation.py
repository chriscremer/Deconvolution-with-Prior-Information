
import numpy as np

import os,sys,inspect
home = os.path.expanduser('~')
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
# sys.path.insert(0,'..')
# sys.path.insert(0,'../simulated_data')
sys.path.insert(0, home+'/plotting')

import nice_plot

import make_real_simulated_data as mrsd
import project_tools as pt
import alternating_decomposition2 as ad

from sklearn.decomposition import NMF

import solve_for_W as sfw
import solve_for_Z as sfz


from scipy.optimize import nnls


##########################################
# Experiment Set Up
##########################################

D=1000

#Models to test
# models = ['alt_NNLS', 'model_2', 'model_3', 'model_4']
# models = ['alt_NNLS', 'model_L2_L1', 'model_L1_L1']
# models = ['alt_NNLS', 'model_L1_L1']
# models = ['alt_NNLS', 'model_L1_L1', 'model_L1_L2', 'model_L2_L1']
# models = ['model_L1_L1', 'model_L1_L2', 'model_L2_L1', 'model_L2_L2']
models = ['Projection', 'Projection_Fit_Prior', 'Exhaustive_Search']

# models = ['model_L1', 'model_L2']
# models = ['model_L2_L1']
#Number of iterations to average results
n_iters = 10
#Tolarence of convergence
tol = .001
#Regularization strength
lmbd = 0
#Independent variable, see how models perform when this variable changes
# var_name = 'Noise'
# var_list = [0.001,5.,10.,15.]
# var_list = [0.01,.2,.4, .6,.8,1.,]
# var_list = [0.001,.2,.4,.6, .8, 1.]
# var_list = [0.001,.1,.2,.3,.4]
# var_list = [0.001,.5,1.]
# var_name = 'Regularization Strength'
# var_list = [1., 5.,10., 50., 100., 500., 1000.]
# var_list = [1., 5.,10., 50., 100.]
# var_list = [10., 100., 1000.]
# var_name = '# Samples'
# var_list = [1,5,10,20]
var_name = '# Components'
var_list = [5,6,7,8,9,10]
#Data Parameters
n_samps = 50
# n_comps = 10
# p_of_zero=.7
n_comps = 10
p_of_zero=.3
noise1 = 5.
n_outliers = 0
#For storing the results
results_W = []
results_Z = []
#Name of plot file
save_as = 'plots/perm_relax.png'

##########################################
# Begin Experiment
##########################################

for iter_ in range(n_iters):
	print 
	print '\nIter: ' + str(iter_)
	print
	#Store the result for each model for each variable value
	W_error_for_each_model = [[] for x in models]
	Z_error_for_each_model = [[] for x in models]

	for var in var_list:

		print '\nIter: ' + str(iter_) + ' Var: ' + str(var)

		##########################################
		# Make Data
		# subpops: gene expression of the real hidden profiles (KxD matrix)
		# fractions: the fraction of each hidden profile in each sample (NxK matrix)
		# X: gene expressions of the mixed samples (NxD matrix)
		##########################################

		subpops, fractions, X = mrsd.make_and_return(n_subpops=var, n_samps=n_samps, probability_of_zero=p_of_zero, noise=noise1, n_outliers=n_outliers)
		print 'Shape of X: ' + str(X.shape)
		print 'Shape of subpops: ' + str(subpops.shape)
		print 'Shape of fractions: ' + str(fractions.shape)

		X, subpops = pt.preprocess_data(X, subpops, dimensions_to_keep=D, normalize=True)
		print 'Shape of X after preprocessing: ' + str(X.shape) +'\n'


		##########################################
		# Run the models and store results
		##########################################

		Z = subpops
		print fractions[0]

		for m in range(len(models)):

			#solve for W
			W = []

			if models[m] == 'Projection':
				print models[m] 
				for samp in range(len(fractions)):
					# print X[samp].shape
					# print subpops.T.shape
					W.append(nnls(subpops.T, X[samp])[0])
				print W[0]

					# if samp % 3 == 0:
					# 	print fractions[samp]
					# 	print nnls(subpops.T, X[samp])[0]
					# 	print


			if models[m] == 'Projection_Fit_Prior':
				print models[m] 
				for samp in range(len(fractions)):
					non_zero_fractions = []
					for i in range(len(fractions[samp])):
						if fractions[samp][i] > 0:
							non_zero_fractions.append(fractions[samp][i])	
					non_zero_fractions.sort()
					non_zero_fractions = non_zero_fractions[::-1]


					W_nnls = nnls(subpops.T, X[samp])[0]

					non_zero_weights = []
					for i in range(len(W_nnls)):
						if W_nnls[i] > 0:
							non_zero_weights.append(W_nnls[i])	

					fitted_weight_vector = [0.]*len(W_nnls)
					indexes_of_values = np.argsort(W_nnls)[::-1]
					for i in range(len(non_zero_weights)):
						if i >= len(non_zero_fractions):
							break
						fitted_weight_vector[indexes_of_values[i]] = non_zero_fractions[i]
					W.append(fitted_weight_vector)

				print W[0]

					# if samp % 1 == 0:
					# 	print fractions[samp]
					# 	print fitted_weight_vector
					# 	print

					# if len(non_zero_weights) < len(non_zero_fractions):
					# 	print fractions[samp]
					# 	print fitted_weight_vector
					# 	print 'yooyoyo'
					# 	fasfd

			if models[m] == 'Exhaustive_Search':
				print models[m] 
				# for samp in range(len(fractions)):	

				W = pt.assign_W_to_top_5_components(X, fractions, Z)
				print W[0]
				# W.append(W_i)



			##########################################
			# Evaluate error and store it
			##########################################


			L2_error_list= []
			for samp in range(len(fractions)):
				L2_error_list.append(np.linalg.norm(fractions[samp] - W[samp]))
				# print 'error ' + str(samp) + ' '+ str(np.linalg.norm(fractions[samp] - W[samp]))

			Avg_L2_error = np.mean(L2_error_list)


			W_error_for_each_model[m].append(Avg_L2_error)
			# Z_error_for_each_model[m].append(Z_error)

	#This iteration is complete. Store this iteration's results
	results_W.append(W_error_for_each_model)
	# results_Z.append(Z_error_for_each_model)


##########################################
# Analyze the results
##########################################

#Get the average and standard deviation over the iterations for each model
# W_error_avg, Z_error_avg, W_error_std, Z_error_std = pt.get_avg_std(results_W, results_Z, models, var_list, n_iters)

W_error_avg = [[] for x in models]
W_error_std = [[] for x in models]
for m in range(len(models)):
	for var in range(len(var_list)):
		var_for_each_iter_W = []

		for ite in range(n_iters):
			var_for_each_iter_W.append(results_W[ite][m][var])

		W_error_avg[m].append(np.mean(var_for_each_iter_W))
		W_error_std[m].append(np.std(var_for_each_iter_W))


##########################################
# Plot the analysis
##########################################

# nice_plot.plot_lines_2graphs_errorbars2(var_list, W_error_avg, W_error_std, models, var_name, 'W Spearman Cor', 'na', 1, 0)
# nice_plot.plot_lines_2graphs_errorbars2(var_list, Z_error_avg, Z_error_std, models, var_name, 'Z Spearman Cor', save_as, 0, 1)

# nice_plot.plot_lines_2graphs_errorbars2(var_list, W_error_avg, W_error_std, models, var_name, 'W L2 Norm', 'na', 1, 0)
nice_plot.plot_lines_2graphs_errorbars2(var_list, W_error_avg, W_error_std, models, var_name, 'L2 Norm', save_as, 0, 1)

# nice_plot.plot_lines_2graphs_errorbars3_logscale(var_list, W_error_avg, W_error_std, models, 'Regularization Strength', 'W Average Spearman Cor', 'na', 1, 0, xlim=[1,1001])
# nice_plot.plot_lines_2graphs_errorbars3_logscale(var_list, Z_error_avg, Z_error_std, models, 'Regularization Strength', 'Z Average Spearman Cor', save_as, 0, 1, xlim=[1,1001])

# nice_plot.plot_lines_2graphs_errorbars3_logscale(var_list, W_error_avg, W_error_std, models, 'Regularization Strength', 'W L2 Norm', 'na', 1, 0, xlim=[1,1001])
# nice_plot.plot_lines_2graphs_errorbars3_logscale(var_list, Z_error_avg, Z_error_std, models, 'Regularization Strength', 'W L2 Norm', save_as, 0, 1, xlim=[1,1001])






