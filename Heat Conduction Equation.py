import math
import numpy as np
import matplotlib.pyplot as plt

############################################
########## Specify the parameters ########## 
############################################
t_0 = 0.0
t_end = 2400
dt = 0.01

a = 0.00925

r_0 = 0.0
r_end = a
dr = 0.00004625

B = 25000
h = 10
rho = 2750
cp = 800

klin = 0.2

#######################################
########## Set up the arrays ##########
#######################################
t_steps = int((t_end - t_0)/dt) + 1
t = [0]*t_steps

r_steps = int((r_end - r_0)/dr) + 1
r = [0]*r_steps

BE = np.zeros([r_steps, t_steps])
BEval = np.zeros([r_steps, t_steps])
A = np.identity(r_steps)
RHS = np.zeros([r_steps, 1])

########################################
########## Initial Conditions ##########
########################################
t[0] = t_0
for n in range(0, t_steps - 1):  
    t[n+1] = t[n] + dt

r[0] = r_0
for i in range(0, r_steps - 1):  
    r[i+1] = r[i] + dr

for i in range(0, r_steps - 1):
    BE[i, 0] = 0
    BEval[i,0] = 0
    
#######################################
########## Set up the matrix ##########
#######################################
A[0,0] = 1.0 
A[0,1] = -1.0 
A[r_steps - 1, r_steps - 2] = -klin
A[r_steps - 1, r_steps - 1] = klin + h*dr

for i in range(1, r_steps - 1):
    
    #############################################
    ########## Fills the main diagonal ##########
    #############################################
    A[i, i] = 1.0 + (2.0*klin*dt)/(dr**2 * rho*cp)
    
    ############################################
    ########## Fills the sub-diagonal ##########
    ############################################
    if (i > 0):
        A[i, i-1] = (klin*dt)/(2.0*r[i]*dr*rho*cp) - (klin*dt)/(dr**2 * rho*cp)
        
    ##############################################
    ########## Fills the super-diagonal ##########
    ##############################################
    if (i < r_steps - 1):
        A[i, i+1] = (-klin*dt)/(2.0*r[i]*dr*rho*cp) - (klin*dt)/(dr**2 * rho*cp)
        
#####################################
########## Sets up the RHS ##########
#####################################    
for n in range(0, t_steps-1):
    RHS[0] = 0
    g = -200.0*((klin + 0.1*math.exp(-BE[r_steps-1,n]))*math.sin(t[n])* \
        math.sin(200*(r_steps-1))) +(h*math.sin(t[n])*math.cos(200*(r_steps-1)))
    RHS[r_steps - 1] = -0.1*math.exp(-BE[r_steps-1,n])*(BE[r_steps-1,n] \
        - BE[r_steps-2,n]) + g*dr

    for i in range(1, r_steps - 1):
        knon = 0.1*math.exp(-BE[i,n])
        f = rho*cp*math.cos(t[n])*math.cos(200*r[i]) + math.sin(t[n])* \
            math.sin(200*r[i])*((200*klin/r[i]) + (200*knon/r[i]) + \
            (40000*knon*math.sin(t[n])*math.sin(200*r[i]))) + \
            math.sin(t[n])*math.cos(200*r[i])*((40000*klin) + (40000*knon)) - B
       
        RHS[i] = BE[i, n] + ((knon*dt)/(2*r[i]*dr*rho*cp))*(BE[i+1,n] \
            - BE[i-1,n]) + ((knon*dt)/(dr**2 * rho*cp))*(BE[i+1,n] - \
            2.0*BE[i,n] + BE[i-1,n]) - ((knon*dt)/(4.0*dr**2 * rho*cp))* \
            ((BE[i+1,n] - BE[i-1,n])**2) + (B*dt/(rho*cp)) + f*dt/(rho*cp)
   
    BE[:, range(n+1,n+2)] = np.linalg.solve(A, RHS)

for n in range(0, t_steps):
    for i in range(0, r_steps):
        BEval[i,n] = math.sin(t[n])*math.cos(200*r[i]) 

#######################################
########## Creates the plots ##########
####################################### 
#plt.plot(r, BE[:,t_steps-1], linewidth = 3, color = 'green', label = 'BE')
#plt.plot(r, BEval[:,t_steps-1], linewidth = 3, color = 'blue', label = 'BEval')
#plt.title('Heat Conduction Validation', fontsize = 20)
#plt.legend(loc = 4)
#plt.show() 

#plt.title('Validation Error', fontsize = 20)
#plt.plot(r,np.absolute(BE[:,t_steps - 1] - BEval[:,t_steps - 1])
#    ,linewidth = 3,color = 'green')
#plt.xlabel('r', fontsize = 12)
#plt.ylabel('|BE - BEval|', fontsize = 12)
#plt.show()

R, T = np.meshgrid(r, t, indexing = 'ij')
#cp = plt.contourf(R, T, BE, 50)
#cbar = plt.colorbar(cp)
#plt.title('Heat Conduction Validation', fontsize = 20)
#plt.show()

cp = plt.contourf(R, T, np.absolute(BE-BEval), 50)
cbar = plt.colorbar(cp)
plt.title('Validation Error', fontsize = 20)
plt.show()