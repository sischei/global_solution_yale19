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

import matplotlib
from mpl_toolkits.mplot3d import Axes3D 
import pickle 
import time 


#======================================================================
# Start with Value Function Iteration


# Set up a container to save the k/obj pairs as they are sampled & optimised 

start = time.time() 

sample_container = {"k_init":[],
                    "value": [], 
                    "iteration": [], 
                    "consumption": [], 
                    "investment": [], 
                    "labor": []} 


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
        sample_container['iteration'].append(iter_container[j][2])
        sample_container['consumption'].append(iter_container[j][3])
        sample_container['investment'].append(iter_container[j][4])
        sample_container['labor'].append(iter_container[j][5])
    
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
end = time.time() 

def plot_scatterplot(): 

    # for all sampled points (not all will have converged, but will give an approximate view of the surface) 
    sample_container['k_init']=np.array(sample_container['k_init'])
    sample_container['value']=np.array(sample_container['value'])



    matplotlib.use('tkagg')

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')#(projection='3d') 
    ax.set_xlabel("k (sector 1)")
    ax.set_ylabel("k (sector 2)")
    ax.set_zlabel("Value")

    #colormap = matplotlib.cm(sample_container['iteration'])

    img = ax.scatter(sample_container['k_init'][:,0], sample_container['k_init'][:,1] , sample_container['value'],c=sample_container['iteration'])
    plt.colorbar(img)

    #plt.show()
    plt.show()

def get_gaussian_process(): 
    with open("./restart/restart_file_step_"+str(numits-1)+".pcl", 'rb') as fd:
        gp_old = pickle.load(fd)
    
    fd.close() 
    return gp_old 

def get_values(k_inits): 
    Gaussian_Process = get_gaussian_process()
    values = Gaussian_Process.predict(k_inits, return_std=False) 
    return values 



def convergence_check(): 
    # tests for convergence by checking the predicted values at the sampled points of the final 
    # iterate and then testing on the optimized value #v_old - value_test

    k_test = [] 
    value_test = [] 

    # load the final instance of Gaussian Process 

    with open("./restart/restart_file_step_"+str(numits-1)+".pcl", 'rb') as fd:
        gp_old = pickle.load(fd)
    
    fd.close() 
    for i in iter_container: 
        k_test.append(i[0])
        value_test.append(i[1])

    k_test = np.array(k_test) 
    value_test = np.array(value_test) 

    v_old = gp_old.predict(k_test, return_std=False) 

    print("=================== Convergence Check ===================")
    print(" ") 
    print("Should be close to zero for all values")

    print(np.around(v_old-value_test,5)) 

def extract_variables(default=True, k_vals = None): 
    # extract the consumption, investment, labour variables (from the final iteration if default=True) 
    # if false, specify random points and calculate 


    k_test = [] 
    value_test = []
    consumption = [] 
    investment = [] 
    labor = [] 

    if default: 
        for i in iter_container: 
            k_test.append(i[0])
            value_test.append(i[1])
            consumption.append(i[3])
            investment.append(i[4])
            labor.append(i[5])


    k_test = np.array(k_test) 
    value_test = np.array(value_test) 
    consumption = np.array(consumption) 
    investment = np.array(investment) 
    labor = np.array(labor)


    return k_test, value_test, consumption, investment, labor  

convergence_check()
k_test, value_test, consumption, investment, labor = extract_variables() 
#print(consumption)
#print(investment)
#print(labor)
def help(): 
    print(" ========== Finished ==========") 
    print("Time elapsed: ",round(end-start,2)) 

    print("Call variables k_test, value_test, consumption, investment, labor")
    print("Use plot_scatterplot() to visualise") 
    print("Use get_gaussian_process() to get the Gaussian Process")
    print("Predict values for given k_values with get_values(k_values), e.g. values = get_values(k_test)")
    print("For options / prompts type help()")

help()

#plot_scatterplot() 


"""

"""

# test and debug 


