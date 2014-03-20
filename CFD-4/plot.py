import matplotlib.pyplot as plt
from matplotlib import pyplot
import numpy as np
X1_1=[]
Y1_1=[]
Y1_2=[]

X2_1=[]
Y2_1=[]
Y2_2=[]

X3_1=[]
Y3_1=[]
Y3_2=[]

for index, line in enumerate (open ('gsp_rs.dat','r')):
    if index <1:
       continue
    else:
        X1_1.append(line.split()[0])
        Y1_1.append(line.split()[1])
        Y1_2.append(line.split()[2])

for index, line in enumerate (open ('gsl_rs.dat','r')):
    if index <1:
       continue
    else:
        X2_1.append(line.split()[0])
        Y2_1.append(line.split()[1])
        Y2_2.append(line.split()[2])

for index, line in enumerate (open ('gsl_sor_rs.dat','r')):
    if index <1:
       continue
    else:
        X3_1.append(line.split()[0])
        Y3_1.append(line.split()[1])
        Y3_2.append(line.split()[2])


plt.plot(X1_1,Y1_1,color='black', linewidth=2.0,marker='o', markevery=200,label='Residual of X by GSP')
plt.plot(X2_1,Y2_1,color='black', linewidth=2.0,marker='x',markevery=20, label='Residual of X by GSL')
plt.plot(X3_1,Y3_1,color='black', linewidth=2.0,marker='s',markevery=10, label='Residual of X by GSL_SOR')

plt.xlabel('Iteration Steps')
plt.ylabel('Residual')
plt.xlim(-100)
pyplot.yscale('log')
plt.grid('on')
plt.legend(loc='upper right')
plt.title('comparison of convergence speed of GSP,GSL and GSL_SOR')
plt.savefig('rs_x.png')
plt.show()
plt.close()
plt.plot(X1_1,Y1_2,color='black', linewidth=2.0, marker='o', markevery=200, label='Residual of Y by GSP')
plt.plot(X2_1,Y2_2,color='black', linewidth=2.0, marker='x', markevery=20, label='Residual of Y by GSL')
plt.plot(X3_1,Y3_2,color='black', linewidth=2.0, marker='s', markevery=10, label='Residual of Y by GSL_SOR')

plt.xlabel('Iteration Steps')
plt.ylabel('Residual')
plt.xlim(-100)
pyplot.yscale('log')
plt.grid('on')
plt.legend(loc='upper right')
plt.title('comparison of convergence speed of GSP,GSL and GSL_SOR')
plt.savefig('rs_y.png')
plt.show()

