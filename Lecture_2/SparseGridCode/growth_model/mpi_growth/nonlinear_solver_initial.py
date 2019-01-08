#======================================================================
#
#     This routine interfaces with IPOPT
#     It sets the optimization problem for every gridpoint
#     at the beginning of the VFI.
#
#     Simon Scheidegger, 11/16 ; 07/17
#======================================================================

from parameters import *
from ipopt_wrapper import EV_F, EV_GRAD_F, EV_G, EV_JAC_G
import numpy as np
import pyipopt

#======================================================================

def initial(k_init, n_agents):
    # IPOPT PARAMETERS below 
    nvars=3*n_agents
    N=nvars         # number of vars
    M=3*n_agents+1  # number of constraints
    NELE_JAC=N*M
    NELE_HESS=(N**2-N)/2 + N    # number of non-zero entries of Hess matrix

    # check that number of nonlinear equations is consistent 
    if (N!=3*n_agents):
        print("there is an error with the number of non-lin eqs!")
        quit

    # Vector of variables -> solution of non-linear equation system 
    X=np.empty(N)

    LAM=np.empty(M) # multipliers
    G=np.empty(M)   # (in-)equality constraints

    # Vector of lower and upper bounds
    G_L=np.empty(M)
    G_U=np.empty(M)

    X_L=np.empty(N)
    X_U=np.empty(N)

    Z_L=np.empty(N)
    Z_U=np.empty(N)

    # get coords of an individual grid points 
    grid_pt_box=k_init
    X_L[:n_agents]=c_bar
    X_U[:n_agents]=c_up

    X_L[n_agents:2*n_agents]=l_bar
    X_U[n_agents:2*n_agents]=l_up

    X_L[2*n_agents:3*n_agents]=inv_bar
    X_U[2*n_agents:3*n_agents]=inv_up

    # Set bounds for the constraints 
    G_L[:n_agents]=c_bar
    G_U[:n_agents]=c_up

    G_L[n_agents:2*n_agents]=l_bar
    G_U[n_agents:2*n_agents]=l_up

    G_L[2*n_agents:3*n_agents]=inv_bar
    G_U[2*n_agents:3*n_agents]=inv_up

    G_L[3*n_agents]=0.0 # both values set to 0 for equality contraints
    G_U[3*n_agents]=0.0

    # initial guesses for first iteration
    cons_init=0.5*(X_U[:n_agents] - X_L[:n_agents]) + X_L[:n_agents]
    lab_init=0.5*(X_U[n_agents:2*n_agents] - X_L[n_agents:2*n_agents]) + X_L[n_agents:2*n_agents]
    inv_init=0.5*(X_U[2*n_agents:3*n_agents] - X_L[2*n_agents:3*n_agents]) + X_L[2*n_agents:3*n_agents]

    X[:n_agents]=cons_init
    X[n_agents:2*n_agents]=lab_init
    X[2*n_agents:3*n_agents]=inv_init
    #X=np.ones(nvars)
    
    # Create ev_f, eval_f, eval_grad_f, eval_g, eval_jac_g for given k_init and n_agent 
    def eval_f(X):
        return EV_F(X, k_init, n_agents)
    
    def eval_grad_f(X):
        return EV_GRAD_F(X,k_init, n_agents)
    
    def eval_g(X):
        return EV_G(X, k_init, n_agents)
        
    def eval_jac_g(X, flag):
        return EV_JAC_G(X, flag, k_init, n_agents)
        
    # First create a handle for the Ipopt problem 
    nlp=pyipopt.create(nvars, X_L, X_U, M, G_L, G_U, NELE_JAC, NELE_HESS, eval_f, eval_grad_f, eval_g, eval_jac_g)
    nlp.num_option("obj_scaling_factor", -1.00)
    nlp.num_option("tol", 1e-6)
    nlp.num_option("acceptable_tol", 1e-5)
    nlp.str_option("derivative_test", "first-order")
    nlp.str_option("hessian_approximation", "limited-memory")
    nlp.int_option("print_level", 0)
    
    x, z_l, z_u, constraint_multipliers, obj, status=nlp.solve(X)
    nlp.close()
    # x: Solution of the primal variables
    # z_l, z_u: Solution of the bound multipliers
    # constraint_multipliers: Solution of the constraint multipliers
    # obj: Objective value
    # status: Exit Status

    # Unpack Consumption, Labor, and Investment 
    c=x[:n_agents]
    l=x[n_agents:2*n_agents]
    inv=x[2*n_agents:3*n_agents]
    
    to_print=np.hstack((obj,x))
    
    # == debug ==
    #f=open("results.txt", 'a')
    #np.savetxt(f, np.transpose(to_print) #, fmt=len(x)*'%10.10f ')
    #for num in to_print:
    #    f.write(str(num)+"\t")
    #f.write("\n")
    #f.close()
    
    return obj, c, l, inv
    