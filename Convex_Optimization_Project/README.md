#Convex Optimization - Fall 2015 Project

###Permutation Relaxation for Tumour Deconvolution with Prior Proportion Information

Sequenced tumour sample data is used to make predictions about patient prognosis. The accuracy of the predictions are often stifled by the heterogeneous nature of the tumour samples. In an attempt to address this problem, tumour samples can be deconvoluted computationally. With the advent of tumour evolution modelling, we now have access to new information that can be used to improve deconvolution. This information provides the number of different cell populations and their proportions within each sample. Since the assignment of these proportions to specific cell populations is not known, we are faced with the combinatorial optimization problem of finding the optimal permutation that minimizes the residual error. Here, I show that this problem can be approximated by a relaxation to a convex problem and fitting to the nearest permutation.
