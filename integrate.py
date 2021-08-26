import sys
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from parameters import *

#Define set of equations
def diffeqs(y,x,p):

    h1,h2,s = y
    r,t,epst,epss,delt,dels = p

    D00 = -epst*h1/(2*(1-h1))
    D01 = -delt*(1-h1)/(2*t*h1**2)
    D02 = -epss*s*h2**2/(2*t*h1)
    D03 = -dels/(2*t*h1*h2)

    D12 = epss*h2
    D13 = dels/(s*h2**2)
    D14 = -r/(2*s**2*h2)

    D22 = -epss*s
    D23 = -dels/(h2**3)
    D24 = r/(s*h2**2)

    D30 = epst*t*h1**2/(1-h1)
    D31 = delt*(1-h1)/h1

    D42 = epss*s*h2**2
    D43 = dels/h2

    D = D00,D01,D02,D03,D12,D13,D14,D22,D23,D24
    Dpsi = D30,D31,D42,D43

    return D,Dpsi

def fun(y,x,p):

    h1,h2,s = y

    r,t,epst,epss,delt,dels = p
    D,Dpsi = diffeqs(y,x,p)

    D00,D01,D02,D03,D12,D13,D14,D22,D23,D24 = D

    f = [D00+D01+D02+D03,
         D12+D13+D14,
         D22+D23+D24,
        ]

    return f



def integrate(parms):
    #Read parameters
    for i in parms.keys():
        exec(i+"=parms['"+i+"']['basicval']")

    #Scaling terms
    H = Ht
    T = T1-Ta
    S = a*T/b
    X = L
    Psi = g*a*T*H**2/(2*f)
    A0 = 16*Af**2*L**2/(9*np.pi)
    Gamma = g*a*T*Cp*rho*H**2/(2*f*A0)
    R = g*a**2*T**2*H**2/(2*f*L*b*S1)

    #Nondim params
    gamma = Gam/Gamma
    r = Run/R
    h1_0 = Hin/H
    h2_0 = H2_0/H
    s_0 = (S1-S2_0)/S

    epst = ctheta*(Ht-Hs)*L/(Ht*dc)
    delt = 2*k*f*dc*L/(ntheta*g*a*T*(Ht-Hs)*Ht**2)
    dels = 2*k*f*ds*L/(nsigma*g*a*T*Ht**3)
    epss = csigma*L/ds

    #Integration params
    t = .5 #Initial guess for t
    t_err = .1
    t_stop = .0001
    Niter = 30 #Number of rounds to converge t
    x = np.arange(0,1,.0001)

    #Start integration
    y0 = [h1_0,h2_0,s_0]

    p = r,t,epst,epss,delt,dels

    while t_err > t_stop:
        out = odeint(fun,y0,x,args=(p,))
        [h1_1,h2_1,s2_1] = out[-1,:]
        t_new = ((s2_1*h2_1**2-gamma)+np.sqrt((s2_1*h2_1**2-gamma)**2+4*gamma*(h1_0**2-h1_1**2)))/(2*(h1_0**2-h1_1**2))
        t_err = np.abs(t-t_new)
        t = t_new
        if t > 1 or t < 0:
            sys.exit('Error: temperature exceeded range')

    #Get full series
    h1 = out[:,0]
    h2 = out[:,1]
    s = out[:,2]

    #Get scaled diagnostics
    H1s = -H*h1
    H2s = -H*h2
    Ss = S*s
    Psi1s = Psi*t*h1**2*1.e-6
    Psi2s = Psi*s*h2**2*1.e-6
    Psi0s = Psi*(t*(h1_0**2-h1**2)-s*h2**2)*1.e-6

    y = h1,h2,s
    p = r,t,epst,epss,delt,dels
    D,Dpsi = diffeqs(y,x,p)
    D30,D31,D42,D43 = Dpsi
    Dvt = D30*Psi/L*1.e-6
    Dwt = D31*Psi/L*1.e-6
    Dvs = D42*Psi/L*1.e-6
    Dws = D43*Psi/L*1.e-6
    xx = x*L

    output = {}
    output['H1']   = H1s
    output['H2']   = H2s
    output['S']    = S
    output['Psi1'] = Psi1s
    output['Psi2'] = Psi2s
    output['Psi0'] = Psi0s
    output['Dvt']  = Dvt
    output['Dwt']  = Dwt
    output['Dvs']  = Dvs
    output['Dws']  = Dws
    output['xx']   = xx
    output['T'] = t*T

    return output

