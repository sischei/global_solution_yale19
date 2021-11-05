x = 1.2 # assign some value
N = 25
 # maximum power in sum
k = 1
s = x
sign = 1.0

import math

while k < N:
  sign = - sign
  k = k + 2
  term = sign*x**k/math.factorial(k)
  s = s + term
  
print( "sin(%g) = %g (approximation with %d terms)" % (x, s, N))
