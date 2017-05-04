# -*- coding: utf-8 -*-
import math
import numpy as np
import matplotlib.pyplot as plt

############################################
########## Specify the parameters ########## 
############################################
t_0 = 0.0
t_end = 1.5
dt = 0.0025

x_0 = 0.0
x_end = 1.0
dx = 0.005

alpha = 0.05
 
#######################################
########## Set up the arrays ##########
#######################################
t_steps = int((t_end - t_0)/dt) + 1
t = [0]*t_steps

x_steps = int((x_end - x_0)/dx) + 1
x = [0]*x_steps

BD2 = np.zeros([t_steps, x_steps])
BD2val = np.zeros([t_steps, x_steps])
A = np.identity(x_steps)
RHS = np.zeros([x_steps, 1])
U = [0]*x_steps

########################################
########## Initial Conditions ##########
########################################
t[0] = t_0
for n in range(0, t_steps - 1):  
    t[n+1] = t[n] + dt

x[0] = x_0
for i in range(0, x_steps - 1):  
    x[i+1] = x[i] + dx

for i in range(0, x_steps - 1):
    BD2[0, i] = (x[i]**2 - x[i])
  
#######################################
########## Set up the matrix ##########
#######################################
A[0,0] = 1.0
A[x_steps - 1, x_steps - 1] = 1.0

for i in range(1, x_steps - 1):
    
    U[i] = 0.5*math.sin(3*math.pi*x[i])
    #############################################
    ########## Fills the main diagonal ##########
    #############################################
    A[i, i] = 3.0/2.0 + ((2.0*alpha*dt)/(dx**2.0))
    
    ############################################
    ########## Fills the sub-diagonal ##########
    ############################################
    if (i > 0):
        A[i, i-1] = -((dt*U[i])/(2*dx)) - ((alpha*dt)/(dx**2))
        
    ##############################################
    ########## Fills the super-diagonal ##########
    ##############################################
    if (i < x_steps - 1):
        A[i, i+1] = ((dt*U[i])/(2*dx)) - ((alpha*dt)/(dx**2))

#####################################
########## Sets up the RHS ##########
#####################################

for n in range(1,2):
    for i in range(1, x_steps - 1): 
        U[i] = 0.5*math.sin(3*math.pi*x[i])
        f = -math.sin(t[n])*(x[i]**2 - x[i]) + U[i]*math.cos(t[n])* \
                (2.0*x[i] - 1.0) - 2.0*alpha*math.cos(t[n])

        BD2[1, i] = BD2[0, i] + ((alpha * dt)/ dx**2) * (BD2[0, i+1] \
            - (2.0*BD2[0, i]) + BD2[0, i-1]) - ((U[i]*dt/(2.0*dx))* \
            (BD2[0, i+1] - BD2[0, i-1])) + f*dt
             
for n in range(1, t_steps - 1):
    RHS[0] = 0
    RHS[x_steps - 1] = 0

    for i in range(1, x_steps - 1):
        U = 0.5*math.sin(3*math.pi*x[i])
        f = -math.sin(t[n+1])*(x[i]**2 - x[i]) + U*math.cos(t[n+1])* \
            (2.0*x[i] - 1.0) - 2.0*alpha*math.cos(t[n+1])
        RHS[i] = 2.0*BD2[n, i] - 0.5*BD2[n-1, i] + f*dt
        
    BD2[range(n+1,n+2),:] = np.transpose(np.linalg.solve(A, RHS))
 
for n in range(0, t_steps):
    for i in range(0, x_steps): 
        BD2val[n, i] = math.cos(t[n])*(x[i]**2 - x[i])

#######################################
########## Creates the plots ##########
####################################### 
T, X = np.meshgrid(t, x, indexing = 'ij')
cp = plt.contourf(T, X, np.absolute(BD2-BD2val), 50)
cbar = plt.colorbar(cp)
#plt.title('Convection Diffusion Equation', fontsize = 20)
#plt.title('Validation', fontsize=20)
#plt.xlabel('t', fontsize = 12)
#plt.ylabel('x', fontsize = 12)
#plt.plot(BD2[t_steps-1,:],x, linewidth = 3, color='blue', label = 'BD2')
#plt.plot(BD2val[t_steps-1,:],x, linewidth = 3, color='green',label = 'BD2val')
#plt.legend(loc = 4)
#plt.ylabel('|BD2 - BD2val|', fontsize = 12)
plt.title('Validation Error', fontsize=20)
#plt.plot(np.absolute(BD2[t_steps - 1,:] - BD2val[t_steps - 1,:]),x,
#    linewidth = 3,color = 'green')
plt.show()