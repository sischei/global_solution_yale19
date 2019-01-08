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
#     - TASMANIAN (http://tasmanian.ornl.gov/)
#
#     Simon Scheidegger; 07/17
#======================================================================



import nonlinear_solver_initial as solver     #solves opt. problems for terminal VF
import nonlinear_solver_iterate as solviter   #solves opt. problems during VFI
from parameters import *                      #parameters of model
import interpolation as interpol              #interface to sparse grid library/terminal VF
import interpolation_iter as interpol_iter    #interface to sparse grid library/iteration
import postprocessing as post                 #computes the L2 and Linfinity error of the model

import TasmanianSG                            #sparse grid library
import numpy as np
from mpi4py import MPI

#======================================================================

comm=MPI.COMM_WORLD
rank=comm.Get_rank()


t1=MPI.Wtime()
# Start with Value Function Iteration
valnew=[]
if (numstart==0):
    valnew=interpol.sparse_grid(n_agents)
     
    if rank==0:
        for itheta in range(ntheta):
            valnew[itheta].write("valnew_" + str(numstart) + "." + str(itheta) + ".txt") #write file to disk for restart

    comm.Barrier()
    
    if rank!=0:
        for itheta in range(ntheta):
            valnew.append(TasmanianSG.TasmanianSparseGrid())
            valnew[itheta].read("valnew_" + str(numstart) + "." + str(itheta) + ".txt")
# value function during iteration
else:
    for itheta in range(ntheta):
        valnew.append(TasmanianSG.TasmanianSparseGrid())
        valnew[itheta].read("valnew_" + str(numstart) + "." + str(itheta) + ".txt")

  
valold=[]
valold=valnew

for i in range(numstart, numits):
    valnew=[]
    valnew=interpol_iter.sparse_grid_iter(n_agents, valold)
    
    if rank==0:
        for itheta in range(ntheta):
            valnew[itheta].write("valnew_" + str(i+1) + "." + str(itheta) + ".txt")
    
    comm.Barrier()
    
    if rank!=0:
        
        for itheta in range(ntheta):
            valnew.append(TasmanianSG.TasmanianSparseGrid())
            valnew[itheta].read("valnew_" + str(i+1) + "." + str(itheta) + ".txt")
    
    valold=[]
    valold=valnew

t2=MPI.Wtime()

#======================================================================
if rank==0:
    print "time: ", t2-t1, "s"
    print "==============================================================="
    print " "
    print " Computation of a growth model of dimension ", n_agents ," finished after ", numits, " steps"
    print " "
    print "==============================================================="
#======================================================================

# compute errors
if rank==0:
    avg_err=post.ls_error(n_agents, numstart, numits, No_samples)
#======================================================================
if rank==0:
    print "==============================================================="
    print " "
    print " Errors are computed -- see error.txt"
    print " "
    print "==============================================================="
#======================================================================

