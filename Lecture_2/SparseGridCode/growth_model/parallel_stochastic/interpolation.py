#======================================================================
#
#     This routine interfaces with the TASMANIAN Sparse grid
#     The crucial part is
#
#     local_aVals[iI]=solver.initial(local_aPoints[iI], theta_init, n_agents)[0]
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
    
    grid_list=[]
    points_list=[]
    num_points=[]
    
    iDim=n_agents
    
    iNumP1_buf=np.zeros(1,int)
    aVals_gathered=0
    aPoints=0
    
    k_range=np.array([k_bar, k_up])

    ranges=np.empty((n_agents, 2))

    for i in range(n_agents):
        ranges[i]=k_range
        
    
    if (comm_switch!=0 or size<ntheta):
        
        if rank==0:
        
            for i in range(ntheta):
                grid_list.append(TasmanianSG.TasmanianSparseGrid())
                grid_list[i].makeLocalPolynomialGrid(iDim, iOut, iDepth, which_basis, "localp")
                grid_list[i].setDomainTransform(ranges)
        
                points_list.append(grid_list[i].getPoints())
                num_points.append(points_list[i].shape[0])
        
        for itheta in range(ntheta):
            theta_init=theta_range[itheta]
            
            if rank==0:
                aPoints=points_list[itheta]
                iNumP1_buf[0]=num_points[itheta]
                aVals_gathered=np.empty((iNumP1_buf[0], 1))
            
            comm.Barrier()
            comm.Bcast(iNumP1_buf, root=0)
            comm.Barrier()
            
            iNumP1=iNumP1_buf[0]
            nump=iNumP1//size
            r=iNumP1 % size
            
            if rank<r:
                nump+=1
                
            displs_scat=np.empty(size, int)
            sendcounts_scat=np.empty(size, int)
            
            displs_gath=np.empty(size, int)
            sendcounts_gath=np.empty(size, int)
            
            
            for i in range(r):
                
                displs_scat[i]=i*(1+iNumP1//size)*n_agents
                sendcounts_scat[i]=(1+iNumP1//size)*n_agents
                
                
                
                displs_gath[i]=i*(1+iNumP1//size)
                sendcounts_gath[i]=(1+iNumP1//size)
                
            for i in range(r,size):
                displs_scat[i]=(r+i*(iNumP1//size))*n_agents
                sendcounts_scat[i]=(iNumP1//size)*n_agents
                
                displs_gath[i]=(r+i*(iNumP1//size))
                sendcounts_gath[i]=(iNumP1//size)
                

            local_aPoints=np.zeros((nump, n_agents))
            comm.Barrier()
            gridp=open("grid_"+str(theta_init)+"_"+str(rank)+".txt", 'w')
            comm.Scatterv([aPoints, sendcounts_scat, displs_scat, MPI.DOUBLE], local_aPoints)
            
            np.savetxt(gridp, local_aPoints)
            gridp.close()
            
            local_aVals=np.empty([nump, 1])
            
            file=open("comparison0_"+str(theta_init)+"_.txt", 'a')
            for iI in range(nump):
                local_aVals[iI]=solver.initial(local_aPoints[iI], theta_init, n_agents)[0]
                v_and_rank=np.array([[local_aVals[iI], rank]])
                to_print=np.hstack((local_aPoints[iI].reshape(1,n_agents), v_and_rank))
                np.savetxt(file, to_print, fmt='%2.16f')
                
            file.close()
            comm.Gatherv(local_aVals, [aVals_gathered, sendcounts_gath, displs_gath, MPI.DOUBLE])
            
            if rank==0:
                grid_list[itheta].loadNeededPoints(aVals_gathered)
                       
            
    else:
        color=rank % ntheta
        
        theta_comm=comm.Split(color, rank)
        theta_rank=theta_comm.Get_rank()
        theta_size=theta_comm.Get_size()
        
        theta_aVals=0
        
        
        theta_init=theta_range[color]
        
        theta_grid=TasmanianSG.TasmanianSparseGrid()
        
        if theta_rank==0:
            theta_grid.makeLocalPolynomialGrid(iDim, iOut, iDepth, which_basis, "localp")
            theta_grid.setDomainTransform(ranges)
            aPoints=theta_grid.getPoints()
            
            iNumP1_buf[0]=aPoints.shape[0] 
            theta_aVals=np.empty((iNumP1_buf[0], 1))
        
        theta_comm.Barrier()
        theta_comm.Bcast(iNumP1_buf, root=0)
        
        iNumP1=iNumP1_buf[0]
        
        nump=iNumP1//theta_size
        r=iNumP1 % theta_size
        
        
        if theta_rank<r:
            nump+=1
        
        displs_scat=np.empty(theta_size)
        sendcounts_scat=np.empty(theta_size)
        
        displs_gath=np.empty(theta_size)
        sendcounts_gath=np.empty(theta_size)
        
        
        for i in range(r):
            displs_scat[i]=i*(1+iNumP1//theta_size)*n_agents
            sendcounts_scat[i]=(1+iNumP1//theta_size)*n_agents
            
            displs_gath[i]=i*(1+iNumP1//theta_size)
            sendcounts_gath[i]=(1+iNumP1//theta_size)
            
        for i in range(r, theta_size):
            displs_scat[i]=(r+i*(iNumP1//theta_size))*n_agents
            sendcounts_scat[i]=(iNumP1//theta_size)*n_agents
            
            displs_gath[i]=(r+i*(iNumP1//theta_size))
            sendcounts_gath[i]=(iNumP1//theta_size)
            
        theta_aPoints=np.zeros((nump, n_agents))
        
        gridp=open("grid_"+str(theta_init)+"_"+str(theta_rank)+".txt", 'a')
        theta_comm.Scatterv([aPoints, sendcounts_scat, displs_scat, MPI.DOUBLE], theta_aPoints)
        np.savetxt(gridp, theta_aPoints)
        gridp.close()
        
        local_aVals=np.empty([nump, 1])
        
        file=open("comparison0_"+str(theta_init)+"_.txt", 'a')
        for iI in range(nump):
            local_aVals[iI]=solver.initial(theta_aPoints[iI], theta_init, n_agents)[0]
            v_and_rank=np.array([[local_aVals[iI], theta_rank]])
            to_print=np.hstack((theta_aPoints[iI].reshape(1,n_agents), v_and_rank))
            np.savetxt(file, to_print, fmt='%2.16f')
        file.close()
        theta_comm.Gatherv(local_aVals, [theta_aVals, sendcounts_gath, displs_gath, MPI.DOUBLE])
        
        if theta_rank==0:
            theta_grid.loadNeededPoints(theta_aVals)
            theta_grid.write("theta_grid0_"+str(color)+".txt")
        
        comm.Barrier()
        if rank==0:
            for th in range(ntheta):
                grid_list.append(TasmanianSG.TasmanianSparseGrid())
                grid_list[th].read("theta_grid0_"+str(th)+".txt")        
    return grid_list

#======================================================================
