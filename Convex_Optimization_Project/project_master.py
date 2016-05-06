
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

##########################################
# Experiment Set Up
##########################################

D=25

#Models to test
# models = ['alt_NNLS', 'model_2', 'model_3', 'model_4']
# models = ['alt_NNLS', 'model_L2_L1', 'model_L1_L1']
# models = ['alt_NNLS', 'model_L1_L1']
# models = ['alt_NNLS', 'model_L1_L1', 'model_L1_L2', 'model_L2_L1']
models = ['model_L1_L1', 'model_L2_L1']
# models = ['model_L2_L1']
#Number of iterations to average results
n_iters = 1
#Tolarence of convergence
tol = .001
#Regularization strength
lmbd = D
#Independent variable, see how models perform when this variable changes
var_name = 'Noise'
var_list = [0.01,2.,4.]
# var_list = [0.,100., 1000.]
#Data Parameters
n_samps = 40
n_comps = 5
p_of_zero=.3
noise1 = .5
#For storing the results
results_W = []
results_Z = []
#Name of plot file
save_as = 'plots/nov27_1.png'

##########################################
# Begin Experiment
##########################################

for iter_ in range(n_iters):
	print 
	print 'Iter: ' + str(iter_)
	print
	#Store the result for each model for each variable value
	W_error_for_each_model = [[] for x in models]
	Z_error_for_each_model = [[] for x in models]

	for var in var_list:

		print 'Iter: ' + str(iter_) + ' Var: ' + str(var)

		##########################################
		# Make Data
		# subpops: gene expression of the real hidden profiles (KxD matrix)
		# fractions: the fraction of each hidden profile in each sample (NxK matrix)
		# X: gene expressions of the mixed samples (NxD matrix)
		##########################################

		subpops, fractions, X = mrsd.make_and_return(n_subpops=n_comps, n_samps=n_samps, probability_of_zero=p_of_zero, noise=var)
		print 'Shape of X: ' + str(X.shape)
		print 'Shape of subpops: ' + str(subpops.shape)
		print 'Shape of fractions: ' + str(fractions.shape)

		##########################################
		# Preprocess Data
		# To speed up computation, we'll only use the first 500 genes
		# Also, normalize the data by diving each gene by its mean
		##########################################

		X, subpops = pt.preprocess_data(X, subpops, dimensions_to_keep=D, normalize=True)
		print 'Shape of X after preprocessing: ' + str(X.shape) +'\n'

		##########################################
		# Run the models and store results
		##########################################

		# lmbd = var

		for m in range(len(models)):

			if models[m] == 'NMF':
				print models[m]
				decomposer = NMF(n_components=n_comps)
				decomposer.fit(X)
				W = decomposer.transform(X)
				Z = decomposer.components_ 

			if models[m] == 'alt_NNLS':
				print models[m]
				W, Z = ad.model_1(X, n_comps, tol, fractions, subpops)

			if models[m] == 'model_L1_L1':
				print models[m]
				W, Z = ad.model_L1_L1(X, n_comps, tol, fractions, subpops, lmbd)

			if models[m] == 'model_L2_L1':
				print models[m]
				W, Z = ad.model_L2_L1(X, n_comps, tol, fractions, subpops, lmbd)

			if models[m] == 'model_L1_L2':
				print models[m]
				W, Z = ad.model_L1_L2(X, n_comps, tol, fractions, subpops, lmbd)

			##########################################
			# Evaluate error and store it
			##########################################

			W_error, Z_error = pt.match_Z_error(W, Z, fractions, subpops)
			print 'W Spearman Cor = ' + str(W_error)
			print 'Z Spearman Cor = ' + str(Z_error) + '\n'

			W_error_for_each_model[m].append(W_error)
			Z_error_for_each_model[m].append(Z_error)

	#This iteration is complete. Store this iteration's results
	results_W.append(W_error_for_each_model)
	results_Z.append(Z_error_for_each_model)


##########################################
# Analyze the results
##########################################

#Get the average and standard deviation over the iterations for each model
W_error_avg, Z_error_avg, W_error_std, Z_error_std = pt.get_avg_std(results_W, results_Z, models, var_list, n_iters)

##########################################
# Plot the analysis
##########################################

nice_plot.plot_lines_2graphs_errorbars2(var_list, W_error_avg, W_error_std, models, var_name, 'W Spearman Cor', 'na', 1, 0)
nice_plot.plot_lines_2graphs_errorbars2(var_list, Z_error_avg, Z_error_std, models, var_name, 'Z Spearman Cor', save_as, 0, 1)

# nice_plot.plot_lines_2graphs_errorbars3_logscale(var_list, Z_error_avg, Z_error_std, models, 'Regularization Strength', 'W Average Spearman Cor', 'na', 1, 0, xlim=[1,1001])
# nice_plot.plot_lines_2graphs_errorbars3_logscale(var_list, W_error_avg, W_error_std, models, 'Regularization Strength', 'Z Average Spearman Cor', save_as, 0, 1, xlim=[1,1001])