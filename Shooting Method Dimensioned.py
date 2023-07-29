#This code allows you to solve a second-order ODE

#A common situation where you will use this code is to solve the ODE
#for the concetration profile inside a catalyst

#Let's start with a first-order reaction to compare the analytical and numerical solutions

#Importing the requisite libraries
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import fmin
import matplotlib.pyplot as plt

#Let us first consider the analytical solutions

#Analytical solution 1 is for C(r=0) = 0
#Analytical solution 2 is for dCdr(r=0) = 0

R = 0.006#m
Deff = 8E-8 #m2/s
kov = 0.159
#odeint only solves IVPs
#This means we need to use the shooting method to solve the BVP
#The shooting method is simply an optimization problem
#We will solve the problem for a single value of phi to demonstrate the method

 #value for phi that we will use to solve the BVP numerically
r = np.linspace(0.00000001,R,25) #Re-defining r so that we have a clearer plot


#Range of concentration at surface of catalyst was chosen from 0.1 to 0.6 at 0.1 increments.



#The solution array is S = [C, J] where dCdr = J


def bvp1(S, r):
    CAsurf = 0.1
    P = (kov*CAsurf*R**2)/(Deff)
    dCdr = S[1]
    d2Cdr= P/(CAsurf*R**2)*S[0]**2 - 1/r *S[1]
    dSdr = [dCdr, d2Cdr]
    return dSdr

#Shooting method for case of C = 0 at r = 0

def shooting1(J_guess):
    CAsurf = 0.1
    S0 = [0,J_guess]
    S1 = odeint(bvp1, S0, r) 
    C_numerical1 = S1[:,0]
    k = np.size(r)-1
    sol1 = abs(C_numerical1[k] - CAsurf)
    return sol1

def bvp2(S, r):
    CAsurf = 0.2
    P = (kov*CAsurf*R**2)/(Deff)
    dCdr = S[1]
    d2Cdr= P/(CAsurf*R**2)*S[0]**2  - 1/r *S[1]
    dSdr = [dCdr, d2Cdr]
    return dSdr

#Shooting method for case of C = 0 at r = 0

def shooting2(J_guess):
    CAsurf = 0.2
    S0 = [0,J_guess]
    S1 = odeint(bvp2, S0, r) 
    C_numerical1 = S1[:,0]
    k = np.size(r)-1
    sol1 = abs(C_numerical1[k] - CAsurf)
    return sol1

def bvp3(S, r):
    CAsurf = 0.3
    P = (kov*CAsurf*R**2)/(Deff)
    dCdr = S[1]
    d2Cdr= P/(CAsurf*R**2)*S[0]**2  - 1/r *S[1]
    dSdr = [dCdr, d2Cdr]
    return dSdr

#Shooting method for case of C = 0 at r = 0

def shooting3(J_guess):
    CAsurf = 0.3
    S0 = [0,J_guess]
    S1 = odeint(bvp3, S0, r) 
    C_numerical1 = S1[:,0]
    k = np.size(r)-1
    sol1 = abs(C_numerical1[k] - CAsurf)
    return sol1

def bvp4(S, r):
    CAsurf = 0.4
    P = (kov*CAsurf*R**2)/(Deff)
    dCdr = S[1]
    d2Cdr= P/(CAsurf*R**2)*S[0]**2  - 1/r *S[1]
    dSdr = [dCdr, d2Cdr]
    return dSdr

#Shooting method for case of C = 0 at r = 0

def shooting4(J_guess):
    CAsurf = 0.4
    S0 = [0,J_guess]
    S1 = odeint(bvp4, S0, r) 
    C_numerical1 = S1[:,0]
    k = np.size(r)-1
    sol1 = abs(C_numerical1[k] - CAsurf)
    return sol1

def bvp5(S, r):
    CAsurf = 0.5
    P = (kov*CAsurf*R**2)/(Deff)
    dCdr = S[1]
    d2Cdr= P/(CAsurf*R**2)*S[0]**2  - 1/r *S[1]
    dSdr = [dCdr, d2Cdr]
    return dSdr

#Shooting method for case of C = 0 at r = 0

def shooting5(J_guess):
    CAsurf = 0.5
    S0 = [0,J_guess]
    S1 = odeint(bvp5, S0, r) 
    C_numerical1 = S1[:,0]
    k = np.size(r)-1
    sol1 = abs(C_numerical1[k] - CAsurf)
    return sol1

def bvp6(S, r):
    CAsurf = 0.6
    P = (kov*CAsurf*R**2)/(Deff)
    dCdr = S[1]
    d2Cdr= P/(CAsurf*R**2)*S[0]**2  - 1/r *S[1]
    dSdr = [dCdr, d2Cdr]
    return dSdr

#Shooting method for case of C = 0 at r = 0

def shooting6(J_guess):
    CAsurf = 0.6
    S0 = [0,J_guess]
    S1 = odeint(bvp6, S0, r) 
    C_numerical1 = S1[:,0]
    k = np.size(r)-1
    sol1 = abs(C_numerical1[k] - CAsurf)
    return sol1


#Running the optimization

J_ini_correct1 = fmin(shooting1, -20, xtol=1e-8)
J_ini_correct2 = fmin(shooting2, -20, xtol=1e-8)
J_ini_correct3 = fmin(shooting3, -20, xtol=1e-8)
J_ini_correct4 = fmin(shooting4, -20, xtol=1e-8)
J_ini_correct5 = fmin(shooting5, -20, xtol=1e-8)
J_ini_correct6 = fmin(shooting6, -20, xtol=1e-8)

#Resolving the ODE with the correct initial value of J
S_ini1 = [0, J_ini_correct1]
S_ini2 = [0, J_ini_correct2]
S_ini3 = [0, J_ini_correct3]
S_ini4 = [0, J_ini_correct4]
S_ini5 = [0, J_ini_correct5]
S_ini6 = [0, J_ini_correct6]

S_correct1 = odeint(bvp1, S_ini1, r)
S_correct2 = odeint(bvp2, S_ini2, r)
S_correct3 = odeint(bvp3, S_ini3, r)
S_correct4 = odeint(bvp4, S_ini4, r)
S_correct5 = odeint(bvp5, S_ini5, r)
S_correct6 = odeint(bvp6, S_ini6, r)


#Plotting the result
plt.plot(r, S_correct1[:,0], "g--", label = 'Numerical solution CAsurf = 0.1')
plt.plot(r, S_correct2[:,0], "r--", label = 'Numerical solution CAsurf = 0.2')
plt.plot(r, S_correct3[:,0], "b--", label = 'Numerical solution CAsurf = 0.3')
plt.plot(r, S_correct4[:,0], "o--", label = 'Numerical solution CAsurf = 0.4')
plt.plot(r, S_correct5[:,0], "k--", label = 'Numerical solution CAsurf = 0.5')
plt.plot(r, S_correct6[:,0], "p--", label = 'Numerical solution CAsurf = 0.6')
plt.xlabel("Radius $r$ (m)")
plt.xlim(0.001, 0.006)
plt.ylabel("Concentration " r'$C_A$')
#plt.title('$C_A$ = 0 at r = 0 and ' r'$\phi$' " = 2")
plt.legend()
plt.show()

#Plotting the result
plt.plot(r, S_correct1[:,1], "g--", label = 'Numerical solution CAsurf = 0.1')
plt.plot(r, S_correct2[:,1], "r--", label = 'Numerical solution CAsurf = 0.2')
plt.plot(r, S_correct3[:,1], "b--", label = 'Numerical solution CAsurf = 0.3')
plt.plot(r, S_correct4[:,1], "o--", label = 'Numerical solution CAsurf = 0.4')
plt.plot(r, S_correct5[:,1], "k--", label = 'Numerical solution CAsurf = 0.5')
plt.plot(r, S_correct6[:,1], "p--", label = 'Numerical solution CAsurf = 0.6')
plt.xlabel("Radius $r$ (m)")
plt.ylabel("Derivative of Concentration " r'$\frac{dC_A}{dr}$')
plt.xlim(0.001, 0.006)
plt.ylim(0.001, 600)
#plt.title('$C_A$ = 0 at r = 0 and ' r'$\phi$' " = 2")
plt.legend()
plt.show()
