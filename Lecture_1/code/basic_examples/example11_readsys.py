import sys
print "This is the name of the script: ", sys.argv[0]

print "please write the degrees celcius outside:"
C = sys.argv[1]

F= 9*float(C)/5 + 32
print "it is ", F , "  degrees F"


print "Number of arguments: ", len(sys.argv)
print "The arguments are: " , str(sys.argv)


