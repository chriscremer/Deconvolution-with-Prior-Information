


#to plot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np


def plot_bar_chart(types, means, stds, xaxis_labels, save_as_name):

	# the x locations for the groups
	ind = np.arange(len(xaxis_labels))
	# the width of the bars
	width = 0.1

	fig, ax = plt.subplots()

	#random = ax.bar(ind-width, means[0], width, color='purple', yerr=stds[0])
	#pca = ax.bar(ind, means[1], width, color='red', yerr=stds[1])
	#ica = ax.bar(ind+width, means[2], width, color='blue', yerr=stds[2])

	colors = ['red', 'blue', 'green', 'purple']

	for i in range(len(types)):
		ax.bar(ind+(float(width*i)), means[i], width, yerr=stds[i], label=types[i], color=colors[i])

	# add some text for labels, title and axes ticks
	ax.set_ylabel('Average Norm Error')
	ax.set_xlabel('Number of Components')
	#ax.set_title('Scores by group and gender')
	ax.set_xticks(ind+width)
	ax.set_xticklabels(xaxis_labels)
	#ax.legend( (random[0], pca[0], ica[0]), ('Random', 'PCA', 'ICA') )
	#ax.legend( (random[0], pca[0]), types, fontsize=5 )
	ax.legend(fontsize=7)


	plt.savefig(save_as_name)
	print 'Saved plot'


def plot_bar_chart_no_stds(types, means, xaxis_labels, save_as_name):

	# the x locations for the groups
	ind = np.arange(len(xaxis_labels))
	# the width of the bars
	width = 0.1

	fig, ax = plt.subplots()

	#random = ax.bar(ind-width, means[0], width, color='purple', yerr=stds[0])
	#pca = ax.bar(ind, means[1], width, color='red', yerr=stds[1])
	#ica = ax.bar(ind+width, means[2], width, color='blue', yerr=stds[2])

	colors = ['red', 'blue', 'green', 'purple']

	for i in range(len(types)):
		ax.bar(ind+(float(width*i)), means[i], width, label=types[i], color=colors[i])

	# add some text for labels, title and axes ticks
	ax.set_ylabel('Average Norm Error')
	ax.set_xlabel('Number of Components')
	#ax.set_title('Scores by group and gender')
	ax.set_xticks(ind+width)
	ax.set_xticklabels(xaxis_labels)
	#ax.legend( (random[0], pca[0], ica[0]), ('Random', 'PCA', 'ICA') )
	#ax.legend( (random[0], pca[0]), types, fontsize=5 )
	ax.legend(fontsize=7)


	plt.savefig(save_as_name)
	print 'Saved plot'


def plot_line(y_values, x_values, y_label, x_label, save_as_name):

	plt.plot(x_values, y_values, 'bo')
	plt.xlabel(x_label)
	plt.ylabel(y_label)


	plt.savefig(save_as_name)
	print 'Saved plot'



def plot_line_with_text(y_values, x_values, y_label, x_label, save_as_name, text):

	fig = plt.figure()
	ax1 = fig.add_axes((.1,.4,.8,.5))
	ax1.plot(x_values, y_values, 'b.-')
	fig.text(.1,.15,text, size='x-small')
	plt.xlabel(x_label)
	plt.ylabel(y_label)


	plt.savefig(save_as_name)
	print 'Saved plot'

def plot_two_lines_with_text(y_values, y2_values, x_values, y_label, x_label, save_as_name, text):

	fig = plt.figure()
	ax1 = fig.add_axes((.1,.4,.8,.5))
	ax1.plot(x_values, y_values, 'b.-')
	ax1.plot(x_values, y2_values, 'r.-')
	fig.text(.1,.15,text, size='x-small')

	plt.xlabel(x_label)
	plt.ylabel(y_label)


	plt.savefig(save_as_name)
	print 'Saved plot'




def plot_three_lines_with_text(y_values, y2_values, y3_values, x_values, y_label, x_label, save_as_name, text):

	fig = plt.figure()
	ax1 = fig.add_axes((.1,.4,.8,.5))
	ax1.plot(x_values, y_values, 'b.-')
	ax1.plot(x_values, y2_values, 'r.-')
	ax1.plot(x_values, y3_values, 'g.-')
	fig.text(.1,.15,text, size='x-small')

	plt.xlabel(x_label)
	plt.ylabel(y_label)


	plt.savefig(save_as_name)
	print 'Saved plot'


	

def plot_three_lines_with_text_with_xlog_scale(y_values, y2_values, y3_values, x_values, y_label, x_label, save_as_name, text):

	fig = plt.figure()
	ax1 = fig.add_axes((.1,.4,.8,.5))
	ax1.plot(x_values, y_values, 'b.-')
	ax1.plot(x_values, y2_values, 'r.-')
	ax1.plot(x_values, y3_values, 'g.-')
	fig.text(.1,.15,text, size='x-small')

	plt.xlabel(x_label)
	plt.ylabel(y_label)

	ax1.set_xscale('log')

	plt.savefig(save_as_name)
	print 'Saved plot'




def plot_visualize_learning(save_as_name, X, hidden_profiles, W, Z):


	plt.scatter(X.T[0], X.T[1])


	a = len(hidden_profiles)

	colors = np.random.rand(a)
	# area = np.pi * (15 * np.random.rand(a))**2 # 0 to 15 point radiuses
	plt.scatter(hidden_profiles.T[0], hidden_profiles.T[1], s=[120]*a, c=colors, alpha=0.5)

	plt.scatter(Z.T[0], Z.T[1], s=[100]*a, c=colors, alpha=0.5, marker='*')


	plt.savefig(save_as_name)
	print 'Saved plot'


def plot_visualize_learning_iters(save_as_name, X, hidden_profiles, W0, W1, Z0, Z1):


	a = len(hidden_profiles)
	colors = np.random.rand(a)
	fig = plt.figure(1)

	plt.subplot(221)
	plt.scatter(X.T[0], X.T[1])	
	plt.scatter(hidden_profiles.T[0], hidden_profiles.T[1], s=[120]*a, c=colors, alpha=0.5)
	plt.scatter(Z0.T[0], Z0.T[1], s=[100]*a, c=colors, alpha=0.5, marker='*')

	plt.subplot(223)
	plt.scatter(X.T[0], X.T[1])	
	plt.scatter(hidden_profiles.T[0], hidden_profiles.T[1], s=[120]*a, c=colors, alpha=0.5)
	plt.scatter(Z1.T[0], Z1.T[1], s=[100]*a, c=colors, alpha=0.5, marker='*')

	plt.subplot(122)
	plt.tick_params(
	axis='both',          # changes apply to the x-axis
	which='both',      # both major and minor ticks are affected
	bottom='off',      # ticks along the bottom edge are off
	top='off',         # ticks along the top edge are off
	labelbottom='off',
	right='off',
	left='off',
	labelleft='off') # labels along the bottom edge are off
	plt.text(0, 0, 'Learned freqs', size='x-small')
	plt.text(0, .1, str(W0), size='x-small')
	plt.text(0, 0.5, 'Real freqs', size='x-small')
	plt.text(0, .6, str(W1), size='x-small')

	plt.savefig(save_as_name)
	print 'Saved plot'