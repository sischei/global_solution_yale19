#======================================================================
#
#     This routine interfaces with Gaussian Process Regression
#     The crucial part is 
#
#     y[iI] = solver.initial(Xtraining[iI], n_agents,gp_old)[0]  
#     => at every training point, we solve an optimization problem
#
#     check kernels here: https://scikit-learn.org/stable/auto_examples/gaussian_process
#       /plot_gpr_prior_posterior.html#sphx-glr-auto-examples-gaussian-process-plot-gpr-prior-posterior-py
#
#
#     Simon Scheidegger, 01/19
#     Cameron Gordon 11/21 print statements to Python3
#======================================================================

import numpy as np
from parameters import *
import nonlinear_solver_iterate as solver
#import cPickle as pickle
import pickle 

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel, Matern

#======================================================================

def GPR_iter(iteration,save_data=True):

    # iteration container to save the k / obj pairs 
    # note we have multiple sectors so this may be interesting 

    iter_container = []  
    


    
    # Load the model from the previous iteration step
    restart_data = filename + str(iteration-1) + ".pcl"
    with open(restart_data, 'rb') as fd_old:
        gp_old = pickle.load(fd_old)
        print("data from iteration step ", iteration -1 , "loaded from disk")
    fd_old.close()
    
    ##generate sample aPoints
    np.random.seed(666)   #fix seed
    dim = n_agents
    Xtraining = np.random.uniform(k_bar, k_up, (No_samples, dim))
    y = np.zeros(No_samples, float) # training targets    
    
    # solve bellman equations at training points
    for iI in range(len(Xtraining)):
        y[iI], consumption, investment, labor = solver.iterate(Xtraining[iI], n_agents,gp_old) 

        if iI == len(Xtraining-1) and iteration == numits-1: 
            y[iI], consumption, investment, labor = solver.iterate(Xtraining[iI], n_agents,gp_old,final=True) # this pulls out the objective value (i.e. value) 
        iter_container.append([Xtraining[iI],y[iI],iteration, consumption, investment, labor])
    
    #print data for debugging purposes
    #for iI in range(len(Xtraining)):
        #print Xtraining[iI], y[iI]
  
    # Instantiate a Gaussian Process model  
    kernel = RBF() 



    # Instantiate a Gaussian Process model
    #kernel = 1.0 * RBF(length_scale=1.0, length_scale_bounds=(1e-1, 10.0))
    
    #kernel = 1.0 * RBF(length_scale=100.0, length_scale_bounds=(1e-1, 2e2)) \
    #+ WhiteKernel(noise_level=1, noise_level_bounds=(1e-3, 1e+0))   

    #kernel = 1.0 * RBF(length_scale=100.0, length_scale_bounds=(1e-1, 2e2)) 
    #kernel = 1.0 * Matern(length_scale=1.0, length_scale_bounds=(1e-1, 10.0),nu=1.5)
    
    gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=10)

    # Fit to data using Maximum Likelihood Estimation of the parameters
    gp.fit(Xtraining, y)    
     
    ##save the model to a file
    output_file = filename + str(iteration) + ".pcl"
    print(output_file)
    with open(output_file, 'wb') as fd:
        pickle.dump(gp, fd, protocol=pickle.HIGHEST_PROTOCOL)
        print("data of step ", iteration ,"  written to disk")
        print(" -------------------------------------------")
    fd.close()    
    
    if save_data==True: 
        return iter_container 

#======================================================================
