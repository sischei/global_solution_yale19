#======================================================================
#
#     This routine interfaces with the TASMANIAN Sparse grid
#     The crucial part is
#
#     aVals[iI]=solver.initial(aPoints[iI], n_agents)[0]
#     => at every gridpoint, we solve an optimization problem
#
#     Simon Scheidegger; 07/17
#======================================================================

import TasmanianSG
import numpy as np
from parameters import *
import nonlinear_solver_initial as solver
from mpi4py import MPI

#======================================================================

def sparse_grid(n_agents):
    
    comm=MPI.COMM_WORLD
    rank=comm.Get_rank()
    size = comm.Get_size()
    
    grid  = TasmanianSG.TasmanianSparseGrid()
    
    aPoints=0
    iNumP1_buf=np.zeros(1, int)
    iNumP1=iNumP1_buf[0]
    aVals_gathered=0
    
    
    if rank==0:
        k_range=np.array([k_bar, k_up])

        ranges=np.empty((n_agents, 2))


        for i in range(n_agents):
            ranges[i]=k_range

        iDim=n_agents

        grid.makeLocalPolynomialGrid(iDim, iOut, iDepth, which_basis, "localp")
        grid.setDomainTransform(ranges)

        aPoints=grid.getPoints()
        
        f=open("grid.txt", 'w')
        np.savetxt(f, aPoints, fmt='% 2.5f')
        f.close()
        
        iNumP1=aPoints.shape[0]
        iNumP1_buf[0]=iNumP1
        aVals_gathered=np.empty((iNumP1, 1))
    
    
    comm.Barrier()
    comm.Bcast(iNumP1_buf, root=0)
    iNumP1=iNumP1_buf[0]
    
    
    nump=iNumP1//size
    r=iNumP1 % size
    
    if rank<r:
        nump+=1
    
    displs_scat=np.empty(size)
    sendcounts_scat=np.empty(size)
    
    displs_gath=np.empty(size)
    sendcounts_gath=np.empty(size)
    
    for i in range(r):
        displs_scat[i]=i*(1+iNumP1//size)*n_agents
        sendcounts_scat[i]=(1+iNumP1//size)*n_agents
        
        displs_gath[i]=i*(1+iNumP1//size)
        sendcounts_gath[i]=(1+iNumP1//size)
        
    for i in range(r, size):
        displs_scat[i]=(r+i*(iNumP1//size))*n_agents
        sendcounts_scat[i]=(iNumP1//size)*n_agents
        
        displs_gath[i]=(r+i*(iNumP1//size))
        sendcounts_gath[i]=(iNumP1//size)
        
    local_aPoints=np.zeros((nump, n_agents))
    
    gridp=open("grid"+str(rank)+".txt", 'w')
    comm.Scatterv([aPoints, sendcounts_scat, displs_scat, MPI.DOUBLE], local_aPoints)
    np.savetxt(gridp, local_aPoints)
    gridp.close()
    
    
    local_aVals=np.empty([nump, 1])
    
     
    file=open("comparison0.txt", 'w') 
    for iI in range(nump):
        local_aVals[iI]=solver.initial(local_aPoints[iI], n_agents)[0]
        v_and_rank=np.array([[local_aVals[iI], rank]])
        to_print=np.hstack((local_aPoints[iI].reshape(1,n_agents), v_and_rank))
        np.savetxt(file, to_print, fmt='%2.16f')

    file.close()
    comm.Gatherv(local_aVals, [aVals_gathered, sendcounts_gath, displs_gath, MPI.DOUBLE])
    
    if rank==0:
        grid.loadNeededPoints(aVals_gathered)
        
    return grid

#======================================================================

    