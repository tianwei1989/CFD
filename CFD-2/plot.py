import matplotlib.pyplot as plt
import numpy as np
X1_1=[]
Y1_1=[]
X1_2=[]
Y1_2=[]
X2_1=[]
Y2_1=[]
X2_2=[]
Y2_2=[]
for index, line in enumerate (open ('./Explicit_0.002/results.txt','r')):
    if index <1:
       continue
    else:
        X1_1.append(line.split()[1])
        Y1_1.append(line.split()[2])
for index, line in enumerate (open ('./Explicit_0.00232/results.txt','r')):
    if index <1:
       continue
    else:
        X1_2.append(line.split()[1])
        Y1_2.append(line.split()[2])
for index, line in enumerate (open ('./Implicit_0.002/results.txt','r')):
    if index <1:
       continue
    else:
        X2_1.append(line.split()[1])
        Y2_1.append(line.split()[2])
for index, line in enumerate (open ('./Implicit_0.01/results.txt','r')):
    if index <1:
       continue
    else:
        X2_2.append(line.split()[1])
        Y2_2.append(line.split()[2])
plt.plot(X1_1,Y1_1,color='black', linewidth=2.0,linestyle='-', label='explicit scheme, dt=0.002 s')
plt.plot(X1_2,Y1_2,color='red', linewidth=2.0,linestyle='--', label='explicit scheme, dt=0.00232 s')
plt.plot(X2_1,Y2_1,color='black', linewidth=2.0,linestyle='-.', label='implicit scheme, dt=0.002 s')
plt.plot(X2_2,Y2_2,color='black', linewidth=2.0,linestyle=':', label='implicit scheme, dt=0.01 s')
plt.xlabel('X')
plt.ylabel('Velocity')
plt.ylim(0.00,40.00)
plt.xlim(0.00,0.04)
plt.grid('on')
plt.legend(loc='upper right')
plt.savefig('results.png')
plt.show()
plt.close()

