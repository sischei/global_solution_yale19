#======================================================================
#
#     This routine solves an infinite horizon growth model 
#     with dynamic programming and sparse grids
#
#     The model is described in Scheidegger & Bilionis (2017)
#     https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2927400
#
#     external libraries needed:
#     - IPOPT (https://projects.coin-or.org/Ipopt)
#     - PYIPOPT (https://github.com/xuy/pyipopt)
#     - scikit-learn GPR (https://scikit-learn.org)
#
#     Simon Scheidegger, 01/19 
#
#     Cameron Gordon, 11/21 - conversion to Python3    
#======================================================================

import nonlinear_solver_initial as solver     #solves opt. problems for terminal VF
import nonlinear_solver_iterate as solviter   #solves opt. problems during VFI
from parameters import *                      #parameters of model
import interpolation as interpol              #interface to sparse grid library/terminal VF
import interpolation_iter as interpol_iter    #interface to sparse grid library/iteration
import postprocessing as post                 #computes the L2 and Linfinity error of the model
import numpy as np
import matplotlib.pyplot as plt 

#======================================================================
# Start with Value Function Iteration


# Set up a container to save the k/obj pairs as they are sampled & optimised 

sample_container = {"k_init":[],
                    "value": []} 


for i in range(numstart, numits):
# terminal value function
    if (i==1):
        print("start with Value Function Iteration")
        iter_container = interpol.GPR_init(i)

    
    else:     
        print("Now, we are in Value Function Iteration step", i)
        iter_container = interpol_iter.GPR_iter(i)
    
    for j in range(len(iter_container)): 
        sample_container['k_init'].append(iter_container[j][0]) 
        sample_container['value'].append(iter_container[j][1]) 
    
#======================================================================
print("===============================================================")
print(" ")
print(" Computation of a growth model of dimension ", n_agents ," finished after ", numits, " steps")
print(" ")
print("===============================================================")
#======================================================================

# compute errors   
avg_err=post.ls_error(n_agents, numstart, numits, No_samples_postprocess)

#======================================================================
print("===============================================================")
print(" ")
#print " Errors are computed -- see error.txt"
print(" ")
print("===============================================================")
#======================================================================

sample_container['k_init']=np.array(sample_container['k_init'])
sample_container['vallue']=np.array(sample_container['value'])
import matplotlib

from mpl_toolkits.mplot3d import Axes3D

matplotlib.use('tkagg')


print("k_init",sample_container['k_init'][:,0])
print('value', sample_container['value'])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')#(projection='3d')

ax.scatter(sample_container['k_init'][:,0], sample_container['k_init'][:,1] , sample_container['value'])
#plt.colorbar()
#plt.show()
plt.show()
