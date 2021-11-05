Cdegrees = []
n = 21
C_min = -10
C_max = 40
dC = (C_max - C_min)/float(n-1)
 # increment in C
for i in range(0, n):
  C = -10 + i*dC
  Cdegrees.append(C)
  Fdegrees = []

for C in Cdegrees:
  F = (9.0/5)*C + 32
  Fdegrees.append(F)

for i in range(len(Cdegrees)):
  C = Cdegrees[i]
  F = Fdegrees[i]
  print "%5.1f %5.1f" % (C, F)
