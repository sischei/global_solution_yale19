#======================================================================
# 
#     sets the parameters for the model
#     "Growth Model"
#
#     The model is described in Scheidegger & Bilionis (2017)
#     https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2927400
#
#     Simon Scheidegger, 01/19
#====================================================================== 

import numpy as np

#====================================================================== 

# How many training points for GPR
n_agents= 2  # number of continuous dimensions of the model
No_samples = 10*n_agents

# control of iterations
numstart = 1   # which is iteration to start (numstart = 1: start from scratch, number=/0: restart)
numits = 7    # which is the iteration to end

filename = "restart/restart_file_step_"  #folder with the restart/result files

#====================================================================== 

# Model Paramters

beta=0.8
rho=0.95
zeta=0.5
psi=0.36
gamma=2.0
delta=0.025
eta=1
big_A=(1.0-beta)/(psi*beta)

# Ranges For States
k_bar=0.2
k_up=3.0
range_cube = k_up - k_bar # range of [0..1]^d in 1D


# Ranges for Controls
c_bar=1e-2
c_up=10.0

l_bar=1e-2
l_up=10.0

inv_bar=1e-2
inv_up=10.0

#====================================================================== 

# Number of test points to compute the error in the postprocessing
No_samples_postprocess = 20





