##### read from file 


a1 = []
a2 = []
a3 = []
a4 = []

with open('data.txt') as f:
    for line in f:
        data = line.split()
        a1.append(int(data[0]))
        a2.append(int(data[1]))
        a3.append(int(data[2]))
        a4.append(int(data[3]))
        
print a1, a2, a3, a4

f.close()





