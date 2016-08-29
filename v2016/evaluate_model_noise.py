
from os.path import expanduser
home = expanduser("~")


import numpy as np
from munkres import Munkres, print_matrix

import make_data as md

import DIFI_Aug2016 as m0
import nsubpopsDIFI_Aug2016 as m1
import NMF_Aug2016 as m2
# import skNMF_Aug2016 as m3
# import regNMF_Aug2016 as m3
import sparseNMF_Aug2016 as m3

from sklearn.decomposition import PCA

# import plotly.plotly as py
# import plotly.graph_objs as go

import plotly
import plotly.graph_objs as go


#Experiment settings
average_over_x_iters = 2
tol = .001
max_iters = 400


#Simulated data settings
n_subpops = 10
n_components = n_subpops

n_samps = 100
probability_of_zero = .80
noise = .1

#Variable
noise_amount=[0.00001, .05, .1, .15, .2 ,.25]




#Run models on it
models = [
			m0.DIFI(n_components=n_components, tol=tol, max_iters=max_iters),
			# m1.nsubpopsDIFI(n_components=n_components, tol=tol, max_iters=max_iters),
			# m2.NMF(n_components=n_components, tol=tol, max_iters=max_iters),
			m3.sparseNMF(n_components=n_components, tol=tol, max_iters=max_iters),

		]


model_names = ['DIFI', 'sparseNMF']		

W_L1_error = [[[] for x in models] for y in range(average_over_x_iters)]

for iter_to_avg in range(average_over_x_iters):

	for noise1 in noise_amount:


		print '\nIter', iter_to_avg, 'var', noise1

		#Make data
		X, subpops, fractions = md.make_and_return(n_subpops, n_samps, probability_of_zero, noise1)
		# print X.shape
		# print subpops.shape
		# print fractions.shape


		#Scale by mean
		mean_list = [np.mean(x) for x in subpops.T]
		for gene in range(len(subpops.T)):
			if mean_list[gene] == 0:
				continue
			subpops.T[gene] = subpops.T[gene]/mean_list[gene]
		pca = PCA(n_components=5)
		pca.fit(subpops)
		pca_subpops = pca.transform(subpops)
		# print pca_subpops.shape


		#remove 0s from fractions
		cleaned_fractions = []
		for samp in range(len(fractions)):
			frac_list = []
			for ji in range(len(fractions[samp])):
				if fractions[samp][ji] > 0:
					frac_list.append(fractions[samp][ji])
			cleaned_fractions.append(frac_list)
		#Split fractions 
		fractions_split = md.split_fractions(cleaned_fractions)




		for m_ in range(len(models)):

			print models[m_]


			W,Z = models[m_].deconvolve(X, fractions_split)

			# print W[:5]
			# print Z[0][:10]

			#Scale by mean
			mean_list = [np.mean(x) for x in Z.T]
			for gene in range(len(Z.T)):
				if mean_list[gene] == 0:
					continue
				Z.T[gene] = Z.T[gene]/mean_list[gene]


			#Calculate error

			#First need to match up the Ws with the real fractions, use Hungarian Algorithm
			print 'Matching..'
			norm_matrix = []
			for learned_weight in range(len(W.T)):
				this_samp = []
				for real_fraction in range(len(fractions.T)):
					# this_samp.append(np.linalg.norm(Z[learned_profile] - subpops[real_profile]))
					this_samp.append(sum(abs(W.T[learned_weight] - fractions.T[real_fraction])))
				norm_matrix.append(this_samp)
			# for list1 in norm_matrix:
			# 	print str(['%.2f' % elem for elem in list1])
			m = Munkres()
			indexes = m.compute(norm_matrix)
			indexes2 = [x[1] for x in indexes]
			rearranged_fractions = fractions.T[indexes2].T
			# rearranged_profiles = subpops[indexes2]
			rearranged_profiles = pca_subpops[indexes2]


			W_error = sum(sum(abs(rearranged_fractions - W)))
			# Z_error = sum(sum(abs(rearranged_profiles - Z))) / sum(sum(abs(rearranged_profiles)))
			Z_error = sum(sum(abs(rearranged_profiles - pca.transform(Z)))) #/ sum(sum(abs(rearranged_profiles)))


			#Save data and plot it
			print W_error
			print Z_error

			W_L1_error[iter_to_avg][m_].append(W_error)




#PLOT result

avgs = [[] for x in models]
stds = [[] for x in models]
#Get average of iterations and std
for m in range(len(models)):

	for var in range(len(noise_amount)):

		values = []

		for iter1 in range(average_over_x_iters):

			values.append(W_L1_error[iter1][m][var])

		avgs[m].append(np.mean(values))
		stds[m].append(np.std(values))


print avgs
print stds

data = []
for m in range(len(models)):

	# Create traces
	data.append(go.Scatter(
	    x = noise_amount,
	    y = avgs[m],
	    error_y=dict(
	        type='data',
	        array=stds[m],
	        visible=True
	        ),
	    mode = 'lines',
	    name = model_names[m]
	))


# # Plot and embed in ipython notebook!
# py.iplot(data, filename='line-mode')




# layout = go.Layout(
# 	xaxis=dict(
# 		title='PC1'
# 		# range=[0,l],
# 		# type='linear',
# 		# dtick=20,
# 	),
# 	yaxis=dict(
# 		title='PC2'
# 		# range=[0,l],
# 		# type='linear',
# 		# dtick=20,
# 		# autorange=False,
# 		# showgrid=True,
# 		# zeroline=False,
# 		# showline=False,
# 		# autotick=True,
# 		# ticks='',
# 		# showticklabels=False
# 	),

# 	# width=950,
# 	# height=800
# )

# print 'plotting'

# fig = go.Figure(data=data, layout=layout)
plotly.offline.plot(data, filename=home + '/storage/noise.html')






