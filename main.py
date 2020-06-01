#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""  This script is to simulate Luenberger observer with the model in 
example 6.2.2 in book Sliding Mode Cotrol Theory and Applications by Christopher Edwars """

__author__ = '{Miguel Angel Pimentel Vallejo}'
__email__ = '{miguel.pimentel@umich.mx}'
__date__= '{15/may/2020}'

#Import the modules needed to run the script.
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
from control.matlab import *

#Fuction tha contains the script model
def obs_model(t,x):
    
    #Linear model parameters
    A = np.matrix([[0,1],[-2,0]])
    B = np.matrix([[0],[1]])
    C = np.matrix([[1,1]])

    #Observer parameters
    P = (10**(15))*np.matrix([[0.844,0],[0,1.688]])
    k = 2
    L = acker (A.T,C.T,[-0.5,-0.6])

    #Declare the list to contain the derivative result
    # xdot position -> derivative meaning
    #[0] -> \.x_1, 
    #[1] -> \.x_2, 
    #[2] -> \.\epsilon_1, 
    #[3] -> \.\epsilon_2,
    #[4] -> \.e_1,
    #[5] -> \.e_2.
    xdot = [0,0,0,0,0,0]

    #Input
    u = 0
    
    #Output error
    y_e = C*np.matrix(x[0:2]).T - C*np.matrix(x[2:4]).T 
        
    #Calculus of the dynamics of the system
    xdot[0:2] = (A*np.matrix(x[0:2]).T + B*u).T.tolist()[0]
    xdot[2:4] = (A*np.matrix(x[2:4]).T + B*u + L.T*(y_e[0,0])  ).T.tolist()[0] 
    xdot[4:6] = ((A- L.T*C)*np.matrix(x[4:6]).T  ).T.tolist()[0] 

    return xdot

#Fuction tha contains the script model with uncertainty
def obs_model_uncer(t,x):
    
    #Linear model parameters
    A = np.matrix([[0,1],[-2,0]])
    B = np.matrix([[0],[1]])
    C = np.matrix([[1,1]])

    #Uncertainty matrix
    A_uncer = np.matrix([[0.2,0.8],[-2.1,0.25]])

    #Observer parameters
    P = (10**(15))*np.matrix([[0.844,0],[0,1.688]])
    k = 2
    L = acker (A.T,C.T,[-0.5,-0.6])

    #Declare the list to contain the derivative result
    # xdot position -> derivative meaning
    #[0] -> \.x_1, 
    #[1] -> \.x_2, 
    #[2] -> \.\epsilon_1, 
    #[3] -> \.\epsilon_2,
    xdot = [0,0,0,0]

    #Input
    u = 0
    
    #Output error
    y_e = C*np.matrix(x[0:2]).T - C*np.matrix(x[2:4]).T 
        
    #Calculus of the dynamics of the system
    xdot[0:2] = (A*np.matrix(x[0:2]).T + B*u).T.tolist()[0]
    xdot[2:4] = (A_uncer*np.matrix(x[2:4]).T + B*u + L.T*(y_e[0,0])  ).T.tolist()[0] 

    return xdot


#Initials condictions for the system
x0 = [-1,1]

#Initials condictions for the observer
xo0 = [0,0]

#Intial error
e0 = x0[0]-xo0[0]
e1 = x0[1]-xo0[1]

#Initial condictions vector
x0 = [x0[0],x0[1],xo0[0],xo0[1],e0,e1]

#Initial time
t0 = 0

#ODE with Runge Kutta
r = ode(obs_model).set_integrator('dopri5',atol=1.e-3,rtol=1.e-3)
r.set_initial_value(x0, t0)

#Final time
tf = 15

#Step size
dt = 0.01

#Create list to save the result of the solver
x1 = [x0[0]]
x2 = [x0[1]]
x3 = [x0[2]]
x4 = [x0[3]]
x5 = [x0[4]]
x6 = [x0[5]]
t = [t0]
P = (10**(15))*np.matrix([[0.844,0],[0,1.688]])
k = 2
u = 0

#Loop to solve the ODE
while r.successful() and r.t < tf:
    r.t+dt
    r.integrate(r.t+dt)
    x1.append(r.y[0])
    x2.append(r.y[1])
    x3.append(r.y[2])
    x4.append(r.y[3])
    x5.append(r.y[4])
    x6.append(r.y[5])
    t.append(r.t)
       
#plot results

#Figure 1 plot the system vs observer
labels = ['$x_1$','$x_2$',r'$\^x_1$',r'$\^x_2$']
plt.figure()
plt.title('System and observer')
est1 = plt.plot(t,x1)
est2 = plt.plot(t,x2)
ob1 = plt.plot(t,x3,'--')
ob2 = plt.plot(t,x4,'--')
plt.xlabel('time')
plt.ylabel('magnitude')
plt.grid(True)
plt.legend(est1 + est2 + ob1 + ob2,labels,loc='lower right')

#Figure 2 plot the error bewteen the observer and system
labels = ['$e_1 $','$e_2$']
plt.figure()
plt.title('Error bewteen the system and observer')
est5 = plt.plot(t,x5)
est6 = plt.plot(t,x6)
plt.xlabel('time')
plt.ylabel('magnitude')
plt.grid(True)
plt.legend(est5 + est6 ,labels)

#Initial conditions
x0 = [x0[0],x0[1],xo0[0],xo0[1]]

#ODE with Runge Kutta
r = ode(obs_model_uncer).set_integrator('dopri5',atol=1.e-3,rtol=1.e-3)
r.set_initial_value(x0, t0)

#Final time
tf = 5

#Step size
dt = 0.01

#Create list to save the result of the solver
x1 = [x0[0]]
x2 = [x0[1]]
x3 = [x0[2]]
x4 = [x0[3]]
e1 = [x0[0]-x0[2]]
e2 = [x0[1]-x0[3]]
t = [t0]
P = (10**(15))*np.matrix([[0.844,0],[0,1.688]])
k = 2
u = 0

#Loop to solve the ODE
while r.successful() and r.t < tf:
    r.t+dt
    r.integrate(r.t+dt)
    x1.append(r.y[0])
    x2.append(r.y[1])
    x3.append(r.y[2])
    x4.append(r.y[3])
    e1.append(r.y[0]-r.y[2])
    e2.append(r.y[1]-r.y[3])
    t.append(r.t)
       
#plot results

#Figure 1 plot the system vs observer
labels = ['$x_1$','$x_2$',r'$\^x_1$',r'$\^x_2$']
plt.figure()
plt.title('System and observer')
est1 = plt.plot(t,x1)
est2 = plt.plot(t,x2)
ob1 = plt.plot(t,x3,'--')
ob2 = plt.plot(t,x4,'--')
plt.xlabel('time')
plt.ylabel('magnitude')
plt.grid(True)
plt.legend(est1 + est2 + ob1 + ob2,labels,loc='lower right')

#Figure 2 plot the error bewteen the observer and system
labels = ['$e_1 $','$e_2$']
plt.figure()
plt.title('Error bewteen the system and observer')
est5 = plt.plot(t,e1)
est6 = plt.plot(t,e2)
plt.xlabel('time')
plt.ylabel('magnitude')
plt.grid(True)
plt.legend(est5 + est6 ,labels)

plt.show()