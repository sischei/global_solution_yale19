import numpy as np
mat=np.matrix([[1, 2, 3],[4, 5, 6],[7, 8, 9]])
print mat
np.savetxt('matrix.txt',mat,fmt='%.2f')