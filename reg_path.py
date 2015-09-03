




import numpy as np
from sklearn import linear_model
from sklearn.linear_model import LassoLars 
from sklearn.linear_model import LinearRegression 

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from numpy.linalg import pinv




X = [	[1.,2.],
		[1.,1.],
		[3.8, 4.]]
y = [1.,1.]
X = np.array(X)
y = np.reshape(np.array(y), (2,1))
print 'X shape ' + str(X.T.shape)
print 'y shape ' + str(y.shape)

# alphas, _, coefs = linear_model.lars_path(X.T, y, method='lasso', verbose=True)
clf = LassoLars(alpha=0.001, fit_intercept=False)
clf.fit(X.T, y)
print 'coef' + str(clf.coef_)
print 'alphas' + str(clf.alphas_)
print 'active' + str(clf.active_)
print 'path' + str(clf.coef_path_)


coefs = clf.coef_path_

xx = np.sum(np.abs(coefs.T), axis=1)
xx /= xx[-1]

plt.plot(xx, coefs.T)
ymin, ymax = plt.ylim()
plt.vlines(xx, ymin, ymax, linestyle='dashed')
plt.xlabel('|coef| / max|coef|')
plt.ylabel('Coefficients')
plt.title('LASSO Path')
plt.axis('tight')

plt.savefig('../plots/lars_path')
print 'Saved plot'


afdsf





















Xqwerwq = [[1.,1.],
		[1.,2.],
		[1.,0.],
		[1.3, 1.]]
Xqwerwq = np.array(Xqwerwq)

Z = [
		[1.,2.],
		[1.,0.],
		[1.3, 1.]]		

X = [1,1]


# test_x = [[1], [2], [3]]
# test_y = [[2], [4], [6]]
# test_x = [[1.,1.,1.3], [2.,0.,1.]]
# test_y = [[1.], [1.]]
# # test_x = [[1.], [2.]]
# # test_y = [[1.], [1.]]
# test_x = np.array(test_x)
# test_y = np.array(test_y)
# clf = LinearRegression(fit_intercept=False)
# clf.fit(test_x, test_y)
# print clf.coef_  #2
# print clf.intercept_




Z = np.array(Z)
X = np.reshape(np.array(X), (2,1))


print X.shape
print Z.shape

clf = LassoLars(fit_intercept=False)
# clf = LinearRegression(fit_intercept=False)
clf.fit(Z.T, X)

print 'coef' + str(clf.coef_)
print 'alphas' + str(clf.alphas_)
print 'active' + str(clf.active_)
print 'path' + str(clf.coef_path_)


afsdasdf


lambda1 = .001

w = np.dot(np.dot(X.T, Z.T), pinv(np.dot(Z, Z.T) + lambda1*np.identity(len(Z))))

print w.shape

print w

print np.dot(w, Z)

asfdsas

plt.scatter(Xqwerwq.T[0], Xqwerwq.T[1])

pred = np.dot(w, Z)


# a = len(pred)
# print a

# colors = np.random.rand(a)

plt.scatter(pred.T[0], pred.T[1], s=[120]*len(pred), c=np.random.rand(len(pred)), alpha=0.5)


plt.savefig('../plots/reg_test')
print 'Saved plot'




