import math
import numpy as np
import matplotlib.pyplot as plt

############################################
########## Specify the parameters ########## 
############################################
t_0 = 0
t_end = 1.0
dt = 0.001

x_0 = 0
x_end = 2.0*math.pi
dx = 0.002
 
#######################################
########## Set up the arrays ##########
#######################################
t_steps = int((t_end - t_0)/dt) + 1
t = [0]*t_steps

x_steps = int((x_end - x_0)/dx) + 1
x = [0]*x_steps

u = np.zeros([x_steps,t_steps]) # u(x,t)-->u[i,n]
uval = np.zeros([x_steps,t_steps]) # u(x,t)-->u[i,n]
########################################
########## Initial Conditions ##########
########################################
t[0] = t_0
for n in range(0, t_steps - 1):  
    t[n+1] = t[n] + dt

x[0] = x_0
for i in range(0, x_steps - 1):  
    x[i+1] = x[i] + dx

for i in range(0,x_steps):
    u[i,0] = math.sin(3.0*x[i])**2
    uval[i,0] = math.sin(3.0*x[i])**2
    
for n in range(0,t_steps):
    u[0,n] = 0.0
    uval[0,n] = 0.0

for n in range(0,t_steps-1):
    for i in range(1,x_steps-1):
        f = 6*math.cos(3*x[i])*(math.sin(3*x[i])**3)*(math.cos(t[n])**2) \
            - math.sin(t[n])*(math.sin(3*x[i])**2)
        u[i,n+1] = u[i,n] - (((u[i,n]*dt)/dx)*(u[i,n] - u[i-1,n])) + f*dt

for n in range(0,t_steps-1):
    for i in range(1,x_steps-1):
        uval[i,n+1] = math.cos(t[n])*(math.sin(3*x[i])**2)

#######################################
########## Creates the plots ##########
####################################### 
#plt.title('Validation', fontsize = 20)
#plt.title('Nonlinear Convection Equation', fontsize = 20)       
#plt.plot(x,u[:,t_steps-1], linewidth = 3, color='green', label = 'u')
#plt.plot(x,uval[:,t_steps-1], linewidth = 3, color='blue', label = 'uval')
#plt.xlabel('x', fontsize = 12)
#plt.ylabel('u and uval', fontsize = 12)
#plt.legend(loc = 4)
#plt.show()

plt.title('Validation Error', fontsize = 20)
#plt.plot(x,np.absolute(u[:,t_steps - 1] - uval[:,t_steps - 1]),
#    linewidth = 3,color = 'green')
#plt.xlabel('x', fontsize = 12)
#plt.ylabel('|u - uval|', fontsize = 12)
#plt.show()

X, T = np.meshgrid(x, t, indexing = 'ij')
cp = plt.contourf(X, T, np.absolute(u-uval), 50)
#plt.title('Nonlinear Convection Equation', fontsize = 20)
cbar = plt.colorbar(cp)
plt.show()

#cp = plt.contourf(X, T, u, 50)
#plt.title('Nonlinear Convection Validation', fontsize = 18)
#cbar = plt.colorbar(cp)
#plt.show()