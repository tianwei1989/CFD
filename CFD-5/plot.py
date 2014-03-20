import matplotlib.pyplot as plt
from matplotlib import pyplot
import numpy as np
X1_1=[]
Y1_1=[]
Y1_2=[]
Y1_3=[]

X2_1=[]
Y2_1=[]
Y2_2=[]

X3_1=[]
Y3_1=[]
Y3_2=[]

for index, line in enumerate (open ('result3.dat','r')):
    if index <1:
       continue
    else:
        X1_1.append(line.split()[0])
        Y1_1.append(line.split()[1])
        Y1_2.append(line.split()[10])
        Y1_3.append(line.split()[20])

# for index, line in enumerate (open ('gsl_rs.dat','r')):
    # if index <1:
       # continue
    # else:
        # X2_1.append(line.split()[0])
        # Y2_1.append(line.split()[1])
        # Y2_2.append(line.split()[2])

# for index, line in enumerate (open ('gsl_sor_rs.dat','r')):
    # if index <1:
       # continue
    # else:
        # X3_1.append(line.split()[0])
        # Y3_1.append(line.split()[1])
        # Y3_2.append(line.split()[2])


plt.plot(X1_1,Y1_1,color='black', linestyle='--',linewidth=2.0)
plt.plot(X1_1,Y1_2,color='black', linewidth=2.0,marker='x', markevery=5)
plt.plot(X1_1,Y1_3,color='black', linewidth=2.0,marker='s', markevery=5)

plt.xlabel('X')
plt.ylabel('Velocity')
plt.ylim(-1, 2)
# pyplot.yscale('log')
# plt.grid('on')
# plt.legend(loc='upper right')
# plt.title('comparison of convergence speed of GSP,GSL and GSL_SOR')
plt.savefig('3.png')
plt.show()
plt.close()


