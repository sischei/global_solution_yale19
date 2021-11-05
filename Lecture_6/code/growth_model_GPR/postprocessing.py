#======================================================================
#
#     This module contains routines to postprocess the VFI 
#     solutions.
#
#     Simon Scheidegger, 01/19
#======================================================================

import numpy as np
from parameters import *
import cPickle as pickle
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import (RBF, Matern, RationalQuadratic,
                                              ExpSineSquared, DotProduct,
                                              ConstantKernel)


#======================================================================    
# Routine compute the errors
def ls_error(n_agents, t1, t2, num_points):
    file=open('errors.txt', 'w')
    
    np.random.seed(0)
    unif=np.random.rand(num_points, n_agents)
    k_sample=k_bar+(unif)*(k_up-k_bar)
    to_print=np.empty((1,3))
        
    for i in range(t1, t2-1):
        sum_diffs=0
        diff = 0
      
        # Load the model from the previous iteration step
        restart_data = filename + str(i) + ".pcl"
        with open(restart_data, 'rb') as fd_old:
            gp_old = pickle.load(fd_old)
            print "data from iteration step ", i , "loaded from disk"
        fd_old.close()      
      
        # Load the model from the previous iteration step
        restart_data = filename + str(i+1) + ".pcl"
        with open(restart_data, 'rb') as fd_new:
            gp_new = pickle.load(fd_new)
            print "data from iteration step ", i+1 , "loaded from disk"
        fd_new.close()        
      
        y_pred_old, sigma_old = gp_old.predict(k_sample, return_std=True)
        y_pred_new, sigma_new = gp_new.predict(k_sample, return_std=True)

        # plot predictive mean and 95% quantiles
        #for j in range(num_points):
            #print k_sample[j], " ",y_pred_new[j], " ",y_pred_new[j] + 1.96*sigma_new[j]," ",y_pred_new[j] - 1.96*sigma_new[j]

        diff = y_pred_old-y_pred_new
        max_abs_diff=np.amax(np.fabs(diff))
        average = np.average(np.fabs(diff))
        
        to_print[0,0]= i+1
        to_print[0,1]= max_abs_diff
        to_print[0,2]= average
        
        np.savetxt(file, to_print, fmt='%2.16f')
        print "==================================="

        
    file.close()
    
    return 
        
#======================================================================
        
        
        
        
    