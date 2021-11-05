#======================================================================
#
#     This routine interfaces with the TASMANIAN Sparse grid
#     The crucial part is 
#
#     aVals[iI]=solver.initial(aPoints[iI], n_agents)[0]  
#     => at every gridpoint, we solve an optimization problem
#
#     Simon Scheidegger, 11/16 ; 07/17
#======================================================================

import TasmanianSG
import numpy as np
from parameters import *
import nonlinear_solver_initial as solver

#======================================================================

def sparse_grid(n_agents, iDepth):
    
    k_range=np.array([k_bar, k_up])

    ranges=np.empty((n_agents, 2))

    for i in range(n_agents):
        ranges[i]=k_range

    iDim=n_agents
    iOut=1
    
    grid_list=[]
    points_list=[]
    num_points=[]
    
    for i in range(ntheta):
        grid_list.append(TasmanianSG.TasmanianSparseGrid())
        grid_list[i].makeLocalPolynomialGrid(iDim, iOut, iDepth, which_basis, "localp")
        grid_list[i].setDomainTransform(ranges)
        
        points_list.append(grid_list[i].getPoints())
        num_points.append(points_list[i].shape[0])

    for itheta in range(ntheta):
        theta_init=theta_range[itheta]
        aPoints=points_list[itheta]
        iNumP1=num_points[itheta]
        aVals=np.empty([iNumP1, 1])
        file=open("comparison0_"+str(theta_range[itheta])+"_.txt", 'w')
        for iI in range(iNumP1):
            aVals[iI]=solver.initial(aPoints[iI], theta_init, n_agents)[0]
            v=aVals[iI]*np.ones((1,1))
            to_print=np.hstack((aPoints[iI].reshape(1,n_agents), v))
            np.savetxt(file, to_print, fmt='%2.16f')    
        file.close()
        grid_list[itheta].loadNeededPoints(aVals)
        
    
    return grid_list

#======================================================================

