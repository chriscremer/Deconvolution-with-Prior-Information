
import numpy as np
import sys
sys.path.insert(0,'../../TCGA_data/')
sys.path.insert(0,'../../kaplan_meier/')



import tcga_data_assembler as tda

X, y = tda.main()

#DIVIDE INTO TRAIN AND TEST
split_at = len(y)*.7
train_x = X[:split_at]
train_y = y[:split_at]
test_x = X[split_at:]
test_y = y[split_at:]
print 'Train len', len(train_y)
print 'Test len', len(test_y)

#PREPROCESS
print 'Preprocessing ...'
from sklearn.preprocessing import MinMaxScaler
scaler = MinMaxScaler()
scaler.fit(train_x)
train_x = scaler.transform(train_x)
test_x = scaler.transform(test_x)

#TRANSFORM
print 'Decomposing ...'
from sklearn import decomposition
decomposer = decomposition.NMF(n_components=20, max_iter=1500, tol=.01)
# from deconvol_model import Deconvol
# decomposer = Deconvol(n_components=20)
decomposer.fit(train_x)
train_x = decomposer.transform(train_x)
test_x = decomposer.transform(test_x)


#PREDICT WHAT GROUP EACH SAMPLE GOES INTO
print 'Training ...'
from sklearn.linear_model import Lasso
model = Lasso()
model.fit(train_x, train_y)
prediction = model.predict(test_x)

group1 = []
group2 = []
dividing_line = np.mean(train_y)
print 'Mean', dividing_line
for i in range(len(prediction)):

	if prediction[i] < dividing_line:
		group1.append(i)
	else:
		group2.append(i)
print 'Group 1 len', len(group1)
print 'Group 2 len', len(group2)

#MAKE A KM PLOT
import km_test
km_test.km_plot2(y[group1], [True]*len(y[group1]), y[group2], [True]*len(y[group2]), '../plots/KM_BRCA_NMF.png')







#CLASSIFICATION

# from sklearn.linear_model import LogisticRegression


# model = LogisticRegression(penalty='l1', C=1.)

# model.fit(train_x, train_y)
# print model.score(test_x, test_y)
# prediction = model.predict_proba(test_x)
# for i in range(len(test_y)):
# 	print 'Predict: ' + str(prediction[i][1]) + ' Actual: ' + str(test_y[i])


