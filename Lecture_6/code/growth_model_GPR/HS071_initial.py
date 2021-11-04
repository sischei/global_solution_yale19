

from ipopt_wrapper_A import EV_F, EV_GRAD_F, EV_G, EV_JAC_G
import numpy as np


class HS071(): 
    """
    Class for the optimization problem to be passed to cyipopt 
    Further optimisations may be possible here by including a hessian (optional param) 
    """

    def __init__(self, X, n_agents, k_init, NELE_JAC, NELE_HESS): 
        self.x = X 
        self.n_agents = n_agents 
        self.k_init = k_init 
        self.NELE_JAC = NELE_JAC 
        self.NELE_HESS = NELE_HESS

    # Create ev_f, eval_f, eval_grad_f, eval_g, eval_jac_g for given k_init and n_agent 
    def eval_f(self, x): 
        return EV_F(x, self.k_init, self.n_agents) 

    def eval_grad_f(self, x): 
        return EV_GRAD_F(x, self.k_init, self.n_agents) 

    def eval_g(self, x): 
        return EV_G(x, self.k_init, self.n_agents) 

    def eval_jac_g(self, x, flag): 
        return EV_JAC_G(x, flag, self.k_init, self.n_agents) 

    def objective(self, x): 
        # Returns the scalar value of the objective given x. 
        return self.eval_f(x) 

    def gradient(self, x): 
        # Returns the gradient fo the objective with respect to x.""" 
        return self.eval_grad_f(x) 

    def constraints(self, x): 
        # Returns the constraints 
        return self.eval_g(x) 

    def jacobian(self, x): 
        # Returns the Jacobian of the constraints with respect to x. 
        return self.eval_jac_g(x, False) 

    def intermediate(self, alg_mod, iter_count, obj_value, inf_pr, inf_du, mu,
                     d_norm, regularization_size, alpha_du, alpha_pr,
                     ls_trials):
        """Prints information at every Ipopt iteration."""

        msg = "Objective value at iteration #{:d} is - {:g}"

        print(msg.format(iter_count, obj_value))