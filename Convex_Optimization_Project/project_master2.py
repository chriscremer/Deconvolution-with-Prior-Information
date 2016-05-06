
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


D=1000

#Models to test
# models = ['alt_NNLS', 'model_2', 'model_3', 'model_4']
# models = ['alt_NNLS', 'model_L2_L1', 'model_L1_L1']
# models = ['alt_NNLS', 'model_L1_L1']
# models = ['alt_NNLS', 'model_L1_L1', 'model_L1_L2', 'model_L2_L1']
# models = ['model_L1_L1', 'model_L2_L1']
# models = ['model_L2_L1']
models = ['Projection', 'Projection_Fit_Prior', 'Exhaustive_Search']
# models = ['Projection', 'Projection_Fit_Prior']

#Number of iterations to average results
n_iters = 5
#Tolarence of convergence
tol = .001
#Regularization strength
# lmbd = D
#Independent variable, see how models perform when this variable changes
var_name = 'Noise'
var_list = [.01,.2,.4,.6,.8,1.]
# var_name = 'Number of Components'
# var_list = [5,10, 15]
# var_list = [0.,100., 1000.]
#Data Parameters
n_samps = 50
n_comps = 5
p_of_zero=.5
# noise1 = .2
n_outliers=0
#For storing the results
results_W = []
results_Z = []
#Name of plot file
save_as = 'plots/dec13_deconv.png'

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

		subpops, fractions, X = mrsd.make_and_return(n_subpops=n_comps, n_samps=n_samps, probability_of_zero=(1.-(4./n_comps)), noise=var, n_outliers=n_outliers)
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

			if models[m] == 'Projection':
				print models[m]
				W, Z = ad.Projection_model(X, n_comps, tol)

			if models[m] == 'Projection_Fit_Prior':
				print models[m]
				W, Z = ad.Projection_fit_model(X, n_comps, tol, fractions)

			if models[m] == 'Exhaustive_Search':
				print models[m]
				W, Z = ad.Exhaustive_Search_model(X, n_comps, tol, fractions)


			##########################################
			# Evaluate error and store it
			##########################################

			W_error, Z_error = pt.match_Z_error(W, Z, fractions, subpops)
			print 'W L2 Norm = ' + str(W_error)
			print 'Z L2 Norm = ' + str(Z_error) + '\n'

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

nice_plot.plot_lines_2graphs_errorbars2(var_list, W_error_avg, W_error_std, models, var_name, 'W L2 Norm', 'na', 1, 0)
nice_plot.plot_lines_2graphs_errorbars2(var_list, Z_error_avg, Z_error_std, models, var_name, 'Z L2 Norm', save_as, 0, 1)

# nice_plot.plot_lines_2graphs_errorbars3_logscale(var_list, Z_error_avg, Z_error_std, models, 'Regularization Strength', 'W Average Spearman Cor', 'na', 1, 0, xlim=[1,1001])
# nice_plot.plot_lines_2graphs_errorbars3_logscale(var_list, W_error_avg, W_error_std, models, 'Regularization Strength', 'Z Average Spearman Cor', save_as, 0, 1, xlim=[1,1001])