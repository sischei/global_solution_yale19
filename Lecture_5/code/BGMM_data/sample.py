"""
Shows how you can sample from the mixture of Gaussians.


"""


import sklearn.mixture
import cPickle as pickle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


# Load the model
with open('density.pcl', 'rb') as fd:
    gm = pickle.load(fd)

# Start sampling
X = gm.sample(n_samples=10000)[0]
print 'Sampled X ', X.shape
print 'Here you go:'

# Let's compare them with the original data
data = np.loadtxt('ergodic_data.txt')
X_orig = data[:, 1:-2]
X_orig = X_orig[data[:, -1] == 1, :]

# I can only plot projections of the data in 2D
# So, let me do all possible scatter plots
d = X.shape[1]
for i in range(d):
    for j in range(i+1, d):
        fig, ax = plt.subplots()
        ax.plot(X[:, i], X[:, j], '.')
        ax.plot(X_orig[:, i], X_orig[:, j], 'x')
        ax.set_xlabel('$x_{%d}$' % i)
        ax.set_ylabel('$x_{%d}$' % j)
        ax.legend(['Samples from estimated density', 'Observations'])
        fig.savefig('compare_i=%d_j=%d.pdf' % (i, j))
        del fig
