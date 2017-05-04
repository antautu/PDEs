import math
import numpy as np
import matplotlib.pyplot as plt

############################################
########## Specify the parameters ########## 
############################################
t_0 = 0.0
t_end = 0.6
dt = 0.001

x_0 = 0.0
x_end = 1.0
dx = 0.002

v = 0.005
 
#######################################
########## Set up the arrays ##########
#######################################
t_steps = int((t_end - t_0)/dt) + 1
t = [0]*t_steps

x_steps = int((x_end - x_0)/dx) + 1
x = [0]*x_steps

BE = np.zeros([x_steps, t_steps])
BEval = np.zeros([x_steps, t_steps])
A = np.identity(x_steps)
RHS = np.zeros([x_steps, 1])

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
    BE[i, 0] = math.sin(6*x[i]*math.pi)
    BEval[i, 0] = math.sin(6*x[i]*math.pi)
      
#######################################
########## Set up the matrix ##########
#######################################
A[0,0] = 1.0
A[x_steps - 1, x_steps - 1] = 1.0

for i in range(1, x_steps - 1):
    
    #############################################
    ########## Fills the main diagonal ##########
    #############################################
    A[i, i] = (2.0*v*dt)/(dx**2) + 1
    
    ############################################
    ########## Fills the sub-diagonal ##########
    ############################################
    if (i > 0):
        A[i, i-1] = -(v*dt)/(dx**2)
        
    ##############################################
    ########## Fills the super-diagonal ##########
    ##############################################
    if (i < x_steps - 1):
        A[i, i+1] = -(v*dt)/(dx**2)

#####################################
########## Sets up the RHS ##########
#####################################    
for n in range(0, t_steps - 1):
    RHS[0] = 0
    RHS[x_steps - 1] = 0

    for i in range(1, x_steps - 1):
        f = -math.sin(t[n])*math.sin(6*x[i]*math.pi) + 6*math.pi*\
           (math.cos(t[n])**2)*math.sin(6*x[i]*math.pi)*\
           math.cos(6*x[i]*math.pi) + v*36*(math.pi**2)*math.cos(t[n])*\
            math.sin(6*x[i]*math.pi)
        RHS[i] = BE[i, n] - (dt/(2.0*dx))*BE[i, n]*(BE[i+1, n] - BE[i-1, n]) \
           # + f*dt
        
    BE[:, range(n+1,n+2)] = np.linalg.solve(A, RHS)

for n in range(0, t_steps):
    for i in range(0, x_steps):
        BEval[i,n] = math.cos(t[n])*math.sin(6*x[i]*math.pi)    

#######################################
########## Creates the plots ##########
####################################### 
X, T = np.meshgrid(x, t, indexing = 'ij')
cp = plt.contourf(X, T, BE, 50)
cbar = plt.colorbar(cp)
plt.title('Nonlinear Burgers Equation', fontsize = 20)
#plt.show()

#cp = plt.contourf(X, T, BEval, 50)
#cbar = plt.colorbar(cp)
#plt.show()

#plt.title('Nonlinear Burgers Validation', fontsize = 20)
plt.plot(x,BE[:,t_steps-1], linewidth = 3, color='black', label = 'BE')
#plt.plot(x,BEval[:,t_steps-1], linewidth = 3, color='blue', label = 'BEval')
#plt.legend(loc = 4)
#plt.show()

#plt.title('Validation Error', fontsize = 20)
#plt.plot(x,np.absolute(BE[:,t_steps - 1] - BEval[:,t_steps - 1])
#    ,linewidth = 3,color = 'green')
#plt.xlabel('x', fontsize = 12)
#plt.ylabel('|BE - BEval|', fontsize = 12)
plt.show()


