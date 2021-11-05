"""
Load Simon's data and fit estimate their density using a mixture of Gaussians.

"""


import sklearn
import sklearn.mixture
from sklearn.mixture import BayesianGaussianMixture
import numpy as np
import cPickle as pickle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


data = np.loadtxt('ergodic_data.txt')
X = data[:, 1:-2]
y = data[:, -1]

# Print something
print 'Total number of observations:', X.shape[0]
print 'Dimensionality of X:', X.shape[1]

# Extract points that are labeled by a 1:
X = X[y==1.0]
print 'Total number of observations (excluding y=0):', X.shape[0]

# Start building the mixture of Gaussians to estimate the density
bgm = BayesianGaussianMixture(n_components=20, max_iter=1000, verbose=1)
bgm.fit(X)

# We started with 20 components, but not all of them will be important.
# Some of them are given, very low probabilities.
# We visualize this in this plot
fig, ax = plt.subplots()
ax.bar(np.arange(bgm.n_components) + 0.5, bgm.weights_)
fig.savefig('bgm_probability_of_each_component.pdf')

# Save it to a file
with open('density.pcl', 'wb') as fd:
    pickle.dump(bgm, fd, protocol=pickle.HIGHEST_PROTOCOL)
