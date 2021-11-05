#======================================================================
#
#     This module contains a routine to check the accuracy of the
#     interpolation in the first iteration of the VFI
#     solutions.
#
#     Simon Scheidegger, 11/16 ; 07/17
#======================================================================

import TasmanianSG
import numpy as np

import interpolation as interpol
from parameters import *
import nonlinear_solver_initial as solver

#======================================================================

def test_sg(n_agents, iDepth, num_points=[]):
    #unif=np.random.rand(num_points, n_agents)
    #k_sample=k_bar+(unif)*(k_up-k_bar)
    
    grid=interpol.sparse_grid(n_agents, iDepth)
    
    k_sample=grid.getPoints()
    
    v_interpol=grid.evaluateBatch(k_sample)
    
    v_func=np.empty((len(k_sample),1))
    
    for i in range(len(k_sample)):
        v_func[i][0]=solver.initial(k_sample[i], n_agents)[0]
                
    res=np.fabs(v_interpol-v_func)
    
    return res
  
#======================================================================
