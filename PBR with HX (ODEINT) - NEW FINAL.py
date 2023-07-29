#Importing the requisite libraries
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

#Assuming that we have a first order reaction, A > 2B
#Ideal gases
#Specific heats are assumed to be independent of temperature
#No B in the feed
#The reactor is equipped with a heat exchanger
#There is a pressure drop in the system

#Stoichiometry
a = 1
b = 2
i = 1
delta = b - a

#Inlet conditions
yAo = 0.1 #pure A in inlet
CAo = 10 #mol/m3
vo = 2 # m3/s
FAo = CAo*vo
W = 200 #kg
thetaA = 1
thetaB = 0#No species B in inlet
thetaI = 0.9/0.1 #90% of inlet is inert species
eps = (yAo*delta)/a
Tcin = 273#K

#Reaction rate data
A = 96554.6349 # m3/mol.s
EA = 46815.50118 #J/mol.K
R = 8.314 #J/mol.K

#Thermodynamic data
Hrxno = -250000 #J/mol
CpA = 2000 #J/mol.K
CpB = 3000 #J/mol.K
CpI = 5000 #J/mol.K
deltaCp = ((b*CpB + i*CpI - a*CpA - i*CpI))/a
To = 423 #K
Tref = 298 #K

#Heat exchange and coolant data
U = 1000#W/sqm.K
alpha_c = 5000 #surface area per unit volume m-1
mc = 500*1E0#coolant mass flow rate (kg/s)
CpCool = 200000 #J/kg.K

#Ergun parameters
alpha = 5.0E-4 #1/kg

#omega = 1
Rad = 0.006
rhoB = 2#kg/m3
DA = 1E-7
Deff = 8E-8
deltaBL = 0.1*Rad
kconv = DA/deltaBL
CAso = 0.01
#dCAdr was calculated using shooting method for an initial CA_surf value of 0.1.
dCAdr = 1.32362593e+02
porosity = 0.6
#Define the functions
# S = [XA, T, y, Tc]

def pbr(S, W):
    k = A*np.exp(-EA/(R*S[1]))
    CA = (CAo*(1-S[0])*S[2]*To)/(S[1]*(1+eps*S[0]))
    CAs = (kconv*S[0]+Deff*dCAdr*(1-porosity))/(kconv)
    eta = (2*Deff*dCAdr)/(Rad*k*CAs)
    omega = (kconv*eta)/(kconv + k*eta*(1-porosity)*(Rad/3))
    rA = -k*omega*CA**2
    dXAdW = (1/rhoB)*(-rA)/(FAo)
    dTdW = (1/rhoB)*(-U*alpha_c*(S[1]-S[3]) - (-rA)*(Hrxno + deltaCp*(S[1] - Tref)))/(FAo*((thetaA*CpA + thetaI*CpI + thetaB*CpB) + deltaCp*S[0]))
    dydW = (-alpha/2)*(1/S[2])*(S[1]/To)*(1+eps*S[0])
    dTcdW = (1/rhoB)*(U*alpha_c*(S[1] - S[3]))/(mc*CpCool)
    dSdW = [dXAdW, dTdW, dydW, dTcdW]
    return dSdW                                                                                                          
#Solving the system of ODEs
W = np.linspace(0,200,201)
sol = odeint(pbr,[0,To,1,Tcin],W)
#Extracting the values for x and y from the solution array

XA = sol[:,0]
T = sol[:,1]
y = sol[:,2]
Tc = sol[:,3]
W

#Preparing the plots

plt.plot(W, XA, "k", label = '$X_{A}$')
plt.plot(W, y, "g", label = 'y')
plt.xlabel("Bed weight (kg)")
plt.ylabel("Conversion or dimensionless pressure")
plt.legend()
plt.show()

plt.plot(W, T, "r", label = "Reactor")
plt.plot(W, Tc, "c", label = "Coolant")
plt.xlabel("Bed weight (kg)")
plt.ylabel("Temperature (K)")
plt.ylim(250, 450)
plt.legend()
plt.show()

#XA plots for Varying To
To_values = np.linspace(To, To+500, 501)
XA_values = np.linspace(0,0,501)
total = len(To_values)
i = 0
while(i<total):
    To = To_values[i]
    sol = odeint(pbr,[0,To,1,Tcin],W)
    XA_values[i] = sol[200,0]
    i += 1
plt.plot(To_values, XA_values,"blue")
plt.xlabel("Reactor Temperature (K)")
plt.ylabel("Conversiont $X_A$")
plt.show()