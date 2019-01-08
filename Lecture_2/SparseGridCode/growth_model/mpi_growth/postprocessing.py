#======================================================================
#
#     This module contains routines to postprocess the VFI 
#     solutions.
#
#     Simon Scheidegger, 11/16 ; 07/17
#======================================================================

import TasmanianSG
import numpy as np

from parameters import *

#======================================================================
# Routine to plot sparse grid solutions
def plot_routine(n_agents, grid, dim, num_points):
    f=open("forplot.txt", 'w')
    k_points=0.5*(k_bar + k_up)*np.ones((num_points, n_agents))
    
    k_dim=np.linspace(k_bar, k_up, num_points)
    
    k_points[:,[dim]]=k_dim.reshape(num_points,1)
    
    vals=grid.evaluateBatch(k_points)
    
    to_print=np.hstack((k_points, vals))
    np.savetxt(f, to_print, fmt= '% 2.5f')
    
    f.close()
    return
    
#======================================================================    
# Routine compute the errors
def ls_error(n_agents, t1, t2, num_points):
    file=open('errors.txt', 'w')
    
    np.random.seed(0)
    unif=np.random.rand(num_points, n_agents)
    k_sample=k_bar+(unif)*(k_up-k_bar)
    to_print=np.empty((1,3))
        
    for i in range(t1, t2):
        sum_diffs=0
        diff = 0
      
        v_prev=TasmanianSG.TasmanianSparseGrid()
        v_next=TasmanianSG.TasmanianSparseGrid()
        
        v_prev.read("valnew_1." + str(i) + ".txt")
        v_next.read("valnew_1." + str(i+1) + ".txt")
        
        diff=v_next.evaluateBatch(k_sample) - v_prev.evaluateBatch(k_sample)
        max_abs_diff=np.amax(np.fabs(diff))
        average = np.average(np.fabs(diff))
        
        to_print[0,0]=i+1
        to_print[0,1]=max_abs_diff
        to_print[0,2]= average
        
        np.savetxt(file, to_print, fmt='%2.16f')

        
    file.close()
    
    return 
        
#======================================================================
        
        
        
        
    