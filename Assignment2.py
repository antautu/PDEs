import math
import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.special

############################################
########## Specify the parameters ##########
############################################
a = 2.0
b = 3.0
alpha = 0.5

t_0 = 0.0
t_roast = 1.0  
t_end = 1.2   ### Need to change to 1.0 to get BD2 and Analytical to match up 
dt = 0.001
r_0 = 0.0
r_end = a
dr = 0.002

zeros = 100
#######################################
########## Set up the arrays ##########
#######################################
t_steps = int((t_end - t_0)/dt) + 1
t = [0]*t_steps
    
r_steps = int((r_end - r_0)/dr) + 1
r = [0]*r_steps

Y = [0]*zeros
uAna = [0]*r_steps  

Y = scipy.special.jn_zeros(1, zeros)
uBD2 = np.zeros([r_steps, t_steps])
A = np.identity(r_steps)
RHS = np.zeros([r_steps,1])

########################################
########## Initial Conditions ##########
########################################
t[0] = t_0
for n in range(0, t_steps - 1):  
    t[n+1] = t[n] + dt

r[0] = r_0
for n in range(0, r_steps - 1):  
    r[n+1] = r[n] + dr

for i in range(0, r_steps - 1): 
    uBD2[i,0] = 0
    
#######################################
########## Set up the matrix ##########
#######################################
A[0,0] = -3.0/2.0
A[0,1] = 2.0
A[0,2] = -1.0/2.0
A[r_steps - 1, r_steps - 3] = 1.0/2.0
A[r_steps - 1, r_steps - 2] = -2.0
A[r_steps - 1, r_steps - 1] = 3.0/2.0


for i in range(1, r_steps - 1):
    
    #############################################
    ########## Fills the main diagonal ##########
    #############################################
    A[i,i] =1.0 + (4.0/3.0)*(dt*alpha/(dr**2))

    ############################################
    ########## Fills the sub-diagonal ##########    
    ############################################
    if (i > 0):
        A[i, i-1] = ((dt*alpha)/(3.0*dr*r[i])) - ((2.0/3.0)*(dt*alpha/(dr**2)))
         
    ##############################################   
    ########## Fills the super-diagonal ##########
    ##############################################
    if (i < r_steps - 1):
        A[i, i+1] = -((dt*alpha/(3.0*dr*r[i])) + ((2.0/3.0)*(dt*alpha/(dr**2))))
        
#####################################
########## Sets up the RHS ##########
#####################################
uBD2[0, 1] = uBD2[1, 1]
uBD2[r_steps - 1, 1] = uBD2[r_steps - 2, 1] + b*dr

for i in range(1, r_steps - 1): 
    uBD2[i, 1] = uBD2[i,0]+(alpha*dt)*(((uBD2[i+1,0]-2.0*uBD2[i,0]+uBD2[i-1,0])/dr**2) + ((uBD2[i+1,0]-uBD2[i-1,0])/(dr*2.0*r[i])))

for n in range(1, t_steps - 1):
    RHS[0] = 0.0
    
    if (t[n] <= t_roast):
        RHS[r_steps - 1] = b*dr
    else:
        RHS[r_steps - 1] = 0.0
        A[r_steps - 1, r_steps - 3] = 0.0
        A[r_steps - 1, r_steps - 2] = 0.0
        A[r_steps - 1, r_steps - 1] = 1.0
 

    for i in range(1, r_steps - 1):
        RHS[i] = (4.0/3.0)*uBD2[i, n] - (1.0/3.0)*uBD2[i, n - 1]
        
    uBD2[:,range(n+1,n+2)] = np.linalg.solve(A, RHS)


#########################################
########## Analytcial Solution ##########
#########################################
def ana(rf, tf):
    summation = 0.0

    for m in range(0, zeros): 
        summation = summation + (2.0*a*b/(Y[m]**2)*(scipy.special.jv(2, Y[m]))/(scipy.special.jv(0, Y[m])**2))*math.exp((-alpha*(Y[m]**2)*tf/(a**2)))*scipy.special.jv(0, (Y[m]*rf/a))
   
    return (b*(rf**2)/(2.0*a)) + ((2.0*b*alpha)/a)*tf - ((b*a)/4.0) + summation

for i in range(0, r_steps):
    uAna[i] = ana(r[i], t[t_steps-1])

#######################################
########## Creates the plots ##########
#######################################
#plt.title('BD2 and Analytical', fontsize = 20)
#plt.plot(r, uBD2[:, t_steps - 1], linewidth = 3, color = 'green', label = 'BD2')
#plt.plot(r, uAna, linewidth = 3, color = 'blue', label = 'Analytical')
#plt.xlabel('r', fontsize = 12)
#plt.ylabel('BD2 and Analytical', fontsize = 12)
#plt.legend(loc = 2)
#plt.show()

#plt.title('BD2 and Analytical Difference', fontsize = 20)
#plt.plot(r,np.absolute(uBD2[:,t_steps - 1] - uAna),linewidth = 3,color = 'green')
#plt.xlabel('r', fontsize = 12)
#plt.ylabel('|BD2 - Analytical|', fontsize = 12)
#plt.show()

plt.title('BD2 Roast and BD2 End', fontsize = 20)
plt.plot(r, uBD2[:, 1000], linewidth = 3, color = 'blue', label = 'BD2 Roast')
plt.plot(r, uBD2[:, t_steps - 2], linewidth = 3, color = 'green', label = 'BD2 End')
plt.xlabel('r', fontsize = 12)
plt.ylabel('BD2', fontsize = 12)
plt.legend(loc = 2)
plt.show()


