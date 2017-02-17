import math
import numpy as np
import matplotlib.pyplot as plt

############################################
########## Specify the parameters ##########
############################################
a = -1.0
b = 3.0
l = 2.0
alpha = 0.05

t_0 = 0.0
t_end = 1.0
dt = 0.001
x_0 = 0.0
x_end = l
dx = 0.01

#######################################
########## Set up the arrays ##########
#######################################
t_steps = int((t_end - t_0)/dt + 1)
t = [0]*t_steps
    
x_steps = int((x_end - x_0)/dx + 1)
x = [0]*x_steps

uAna = [0]*x_steps    
uBD2 = np.zeros([x_steps, t_steps])
A = np.identity(x_steps)
RHS = np.zeros([x_steps,1])

########################################
########## Initial Conditions ##########
########################################
t[0] = t_0
for n in range(1, t_steps - 1):  
    t[n+1] = t[n] + dt

x[0] = x_0
for n in range(1, x_steps - 1):  
    x[n+1] = x[n] + dx

for i in range(0, x_steps - 1): 
    uBD2[i,0] = 0
    
#######################################
########## Set up the matrix ##########
#######################################
A[0,0] = -1.0
A[0,1] = 1.0
A[x_steps - 1, x_steps - 2] = -1.0
A[x_steps - 1, x_steps - 1] = 1.0

for i in range(1, x_steps - 1):
    
    #############################################
    ########## Fills the main diagonal ##########
    #############################################
    A[i,i] = 1.0 + (4.0*alpha*dt)/(3.0*(dx**2))

    ############################################
    ########## Fills the sub-diagonal ##########
    ############################################
    if (i > 0):
        A[i, i-1] = -(2.0*alpha*dt) / (3.0*(dx**2))

    ##############################################   
    ########## Fills the super-diagonal ##########
    ##############################################
    if (i < x_steps - 1):
        A[i, i+1] = -(2.0*alpha*dt) / (3.0*(dx**2))

#####################################
########## Sets up the RHS ##########
#####################################
uBD2[0, 1] = uBD2[1, 1] - b*dx
uBD2[x_steps - 1, 1] = uBD2[x_steps - 2, 1] + a*dx

for i in range(1, x_steps - 1):
        uBD2[i, 1] = uBD2[i, 0] + ((alpha * dt)/ dx**2) * (uBD2[i+1, 0] - (2.0*uBD2[i, 0]) + uBD2[i-1, 0])
        
for n in range(0, t_steps - 1):
    RHS[0] = b * dx
    RHS[x_steps - 1] = a * dx

    for i in range(1, x_steps - 1):
        RHS[i] = ((4.0 / 3.0) * uBD2[i, n]) - (1.0 / 3.0) * uBD2[i, n-1]
        
    uBD2[:,range(n+1,n+2)] = np.linalg.solve(A, RHS)

#########################################
########## Analytcial Solution ##########
#########################################
def ana(xf, tf):
    summation = 0.0

    for m in range(1, 10000):
        summation = summation + (((-a*l*2.0 * ((-1)**m))/(math.pi**2 * m**2)) +
            ((b*l*2.0)/(math.pi**2 * m**2))) * math.exp(-t[tf]*((m**2 * (math.pi)**2)/l**2)*alpha) * math.cos((m*(math.pi)*xf)/l)
    return ((a-b)/(2.0*l))*(xf**2) + b*xf + ((alpha*(a-b))/l)*t[tf] -(((a+2*b)*l)/6.0) + summation

for i in range(0, x_steps):
    uAna[i] = ana(x[i], t_steps-1)

#######################################
########## Creates the plots ##########
#######################################
plt.title('BD2 and Analytical', fontsize = 20)
plt.plot(x, uBD2[:, t_steps - 1], linewidth = 3, color = 'green', label = 'BD2')
plt.plot(x, uAna, linewidth = 3, color = 'blue', label = 'Analytical')
plt.xlabel('x', fontsize = 12)
plt.ylabel('BD2 and Analytical', fontsize = 12)
plt.legend(loc = 4)
plt.show()

plt.title('BD2 and Analytical Difference', fontsize = 20)
plt.plot(x,abs(uBD2[:,t_steps-1] - uAna),linewidth = 3,color = 'green')
plt.xlabel('x', fontsize = 12)
plt.ylabel('|BD2 - Analytical|', fontsize = 12)
plt.show()



