import math
import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.special

############################################
########## Specify the parameters ##########
############################################
a = 1.75
b = 2.5

x_0 = 0.0
x_end = a
dx = 0.0175

y_0 = 0.0
y_end = b
dy = 0.025

#######################################
########## Set up the arrays ##########
#######################################
x_steps = int((x_end - x_0)/dx + 1)
x = [0]*x_steps
    
y_steps = int((y_end - y_0)/dy + 1)
y = [0]*y_steps

uhom = np.zeros([x_steps, y_steps])
uDaQ = np.zeros([x_steps, y_steps])

u_old = 0.5 + np.zeros([x_steps, y_steps])
u_new = 0.5 + np.zeros([x_steps, y_steps])

########################################
########## Initial Conditions ##########
########################################
x[0] = x_0
for n in range(0, x_steps - 1):  
    x[n+1] = x[n] + dx

y[0] = y_0
for n in range(0, y_steps - 1):   
    y[n+1] = y[n] + dy

for j in range(0, y_steps - 1):
    u_old[0,j] = 0.0
    u_new[0,j] = 0.0
    
for j in range(0, y_steps - 1):
    u_old[x_steps - 1,j] = 1.0
    u_new[x_steps - 1,j] = 1.0
    
for i in range(0, x_steps - 1):
    u_old[i,0] = 0.0
    u_new[i,0] = 0.0
    
for i in range(0, x_steps - 1):
    u_old[i,y_steps - 1] = 1.0
    u_new[i,y_steps - 1] = 1.0
#########################
########## SOR ##########
#########################
change = 1.0
tolerance = 0.01 
w = (2.0 / (1.0 + math.sin(dx*math.pi)))
while change > tolerance:
    change = 0.0
    for i in range(1, x_steps - 1):
        for j in range(1, y_steps - 1):
            du = (1.0/(2.0*(dx**2) + 2.0*(dy**2)))*(((dy**2)*(u_old[i+1,j] + u_old[i-1,j])) + ((dx**2)*(u_old[i,j+1] + u_old[i,j-1])) - ((2.0*(dx**2) + 2.0*(dy**2))*u_old[i, j]))
            u_new[i, j] = u_old[i, j] + w*du
            change = change + abs(du)
        print change
        u_old = u_new
            
            
####################################################
########## Analytcial From Homogenization ##########
####################################################
def anahom(xf, yf):
    summation = 0.0
    
    for m in range(1, 320): 
        c = ((m*math.pi*a)/b)
        summation = summation + (((2.0/(m*math.pi))*math.sin((m*math.pi*yf)/b)) / (math.exp(-2.0*c) - 1)) * ((math.exp(((m*math.pi*xf)/b) - c)*(math.cos(m*math.pi)*math.exp(-c) - 1)) + (math.exp(((-m*math.pi*xf)/b) - c)*(1 - math.cos(m*math.pi)*math.exp(c))))
    return (1/b)*yf + summation 
    
for j in range(0, y_steps):
    for i in range(0, x_steps):
        uhom[i,j] = anahom(x[i], y[j])
        
########################################################
########## Analytical From Divide-and-Conquer ##########
########################################################
def anaDaQ(xf, yf):
    summation = 0.0
    
    for m in range(1, 320):
        cv = ((m*math.pi*a)/b)
        cw = ((m*math.pi*b)/a)
        v = ((2.0/(m*math.pi))*(1 - math.cos(m*math.pi)) / (1 - math.exp(-2.0*cv))) * math.sin((m*math.pi*yf)/b)*(math.exp(((m*math.pi*xf)/b) - cv) - math.exp(((-m*math.pi*xf)/b) - cv))
        w = ((2.0/(m*math.pi))*(1 - math.cos(m*math.pi)) / (1 - math.exp(-2.0*cw))) * math.sin((m*math.pi*xf)/a)*(math.exp(((m*math.pi*yf)/a) - cw) - math.exp(((-m*math.pi*yf)/a) - cw))
        summation = summation + v + w  
    return summation
    
for j in range(0, y_steps):
    for i in range(0, x_steps):
        uDaQ[i,j] = anaDaQ(x[i], y[j])
        
#######################################
########## Creates the plots ##########
#######################################   
levels = np.linspace(0.0, 1.2, 50)
X, Y = np.meshgrid(x, y, indexing = 'ij')
plt.subplot(131)
plt.title('Analytical From Homogenization', fontsize = 18)
plt.xlabel('x', fontsize = 12)
plt.ylabel('y', fontsize = 12)
cp = plt.contourf(X, Y, uhom, levels)
plt.gca().set_aspect('equal', adjustable='box')
#cbar = plt.colorbar(cp)
#cbar.set_label('Heat', fontsize = 18)
 
plt.subplot(132)
plt.title('Analytical From Divide-and-Conquer', fontsize = 18)
plt.xlabel('x', fontsize = 12)
plt.ylabel('y', fontsize = 12)
cp = plt.contourf(X, Y, uDaQ, levels)
plt.gca().set_aspect('equal', adjustable='box')
#cbar = plt.colorbar(cp)
#cbar.set_label('Heat', fontsize = 18)
  
plt.subplot(133)
plt.title('Numerical SOR', fontsize = 18)
plt.xlabel('x', fontsize = 12)
plt.ylabel('y', fontsize = 12)
cp = plt.contourf(X, Y, u_new, levels)
plt.gca().set_aspect('equal', adjustable='box')
cax = plt.axes([0.125, 0.075, 0.775, 0.03])
cbar = plt.colorbar(cp, cax = cax, orientation = 'horizontal')
cbar.set_label('Heat', fontsize = 18)
plt.show()

#plt.subplot(121)
#plt.title('SOR and Homogenization Difference', fontsize = 18)
#plt.xlabel('x', fontsize = 12)
#plt.ylabel('|SOR - Homogenization|', fontsize = 12)
#cp = plt.contourf(X, Y, abs(u_new - uhom), levels)
#cbar = plt.colorbar(cp)
#cbar.set_label('Error', fontsize = 18)

#plt.subplot(122)
#plt.title('SOR and Divide-and-Conquer Difference', fontsize = 18)
#plt.xlabel('x', fontsize = 12)
#plt.ylabel('|SOR - Divide-and-Conquer|', fontsize = 12)
#cp = plt.contourf(X, Y, abs(u_new - uDaQ), levels)
#cbar = plt.colorbar(cp)
#cbar.set_label('Error', fontsize = 18)
#plt.show()      