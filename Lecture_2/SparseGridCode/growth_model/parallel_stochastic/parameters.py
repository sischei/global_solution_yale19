#======================================================================
#
#     sets the parameters for the model
#     "Growth Model"
#
#     The model is described in Scheidegger & Bilionis (2017)
#     https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2927400
#
#     Simon Scheidegger; 07/17
#======================================================================

import numpy as np

#======================================================================

# Depth of "Classical" Sparse grid
iDepth = 3     # how many levels
iOut=1         # how many outputs
which_basis = 1 #linear basis function (2: quadratic local basis)

# control of iterations
numstart = 0
numits   = 1

# How many random points for computing the errors
No_samples = 1000

#====================================================================== 

# Comm Switch --parallelize over the discrete states 1: on, 0: off 
comm_switch = 0 

#====================================================================== 

# Model Parameters
n_agents = 4

beta=0.8
rho=0.95
zeta=0.5
psi=0.36
gamma=2.0
delta=0.025
eta=1
big_A=(1.0-beta)/(psi*beta)

# Ranges For States
range_cube=1 # range of [0..1]^d in 1D
k_bar=0.2
k_up=3.0

ntheta=5
theta_range=np.array([0.9, 0.95, 1, 1.05, 1.1])

# Ranges for Controls
c_bar=1e-2
c_up=10000.0

l_bar=1e-2
l_up=1.0

inv_bar=1e-2
inv_up=10000.0

#======================================================================
#======================================================================
