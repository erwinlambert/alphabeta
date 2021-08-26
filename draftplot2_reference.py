from utils2 import *
import numpy as np
import matplotlib.pyplot as plt

#Define input parameters
dx = 1.
L = 5.e6
x = np.arange(0.,L,dx)

R = .02+0.*x

#Solve equations
Hsg,S2,Msg = solvePW(x,R)
Hth,T1,FTint,T0 = solveAW(x,Msg)

#Get plotcolors
red,blu,pur = getcols()

#Make variable plot
prettyplot(6,8)
fig,ax = plt.subplots(3,1)

x = x/1000.

ax[0].plot(x,T1,color=red,label=r'$T_1$')
ax[0].plot(x,T0+0.*x,color=pur,label=r'$T_0$')
ax[1].plot(x,S2,color=blu,label=r'$S_2$')
ax[2].plot(x,Hth,color=red,label=r'$H_\theta$')
ax[2].plot(x,Hsg,color=blu,label=r'$H_\sigma$')
ax[2].invert_yaxis()

for AX in ax:
   AX.legend(ncol=2)

for AX in ax[0:2]:
   AX.set_xticklabels([])

ax[2].set_xlabel('Distance [km]')

ax[0].set_ylabel('Temperature [C]')
ax[1].set_ylabel('Salinity [g/kg]')
ax[2].set_ylabel('Pycnocline depth [m]')

ax[0].set_yticks([0,2,4,6,8,10])
ax[1].set_yticks([33,34,35])


saveshow('reference')

#Make transport plot
Psi1,Psi2,Psi0 = transport(Hth,Hsg,T1,S2,T0)

prettyplot(8,6)

fig,ax = plt.subplots(2,1)

ax[0].plot(x,Psi0*1.e-6,color=pur,label=r'DW')
ax[0].plot(x,Psi1*1.e-6,color=red,label=r'AW')
ax[0].plot(x,Psi2*1.e-6,color=blu,label=r'PW')

#ax[1].fill_between(dim(x,'x'),0.*x,dim(dftdx,'ft'),color='.5',alpha=.5)
#ax[1].plot(dim(x,'x'),dim(dftdx,'ft'),color='.5')

ax[0].legend()

ax[1].set_xlabel('Distance [km]')
ax[0].set_ylabel('Volume transport [Sv]')
ax[1].set_ylabel('Lateral heat flux [TW/1000km]')

saveshow('reftransport')


