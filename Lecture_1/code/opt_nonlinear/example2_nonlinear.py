import numpy as np
from scipy.optimize import root

def func2(x):
  f = [x[0] * np.cos(x[1]) - 4, x[1]*x[0] - x[1] - 5]
  df = np.array([[np.cos(x[1]), -x[0] * np.sin(x[1])],[x[1], x[0] - 1]])
  return f, df


sol = root(func2, [1, 1], jac=True, method='lm')
solution = sol.x

print "the solution of this nonlinear set of equations is: ", solution