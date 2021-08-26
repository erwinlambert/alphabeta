from utils2 import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys

#Define input parameters
dx = 1.
L = 5.e6
x = np.arange(0.,L,dx)
rr = np.arange(.001,.06,.002)

Psi1 = np.zeros((len(x),len(rr)))
Psi2 = np.zeros((len(x),len(rr)))
Psi0 = np.zeros((len(x),len(rr)))

for r,RR in enumerate(rr):
   R = RR+0.*x

   #Solve equations
   Hsg,S2,Msg = solvePW(x,R)
   Hth,T1,FTint,T0 = solveAW(x,Msg)

   #Get transports
   Psi1[:,r],Psi2[:,r],Psi0[:,r] = transport(Hth,Hsg,T1,S2,T0)


#Get plotcolors
red,blu,pur = getcols()

#Make variable plot
prettyplot(10,3)
mpl.rcParams['figure.subplot.left']   = .07
mpl.rcParams['figure.subplot.right']  = .95

fig,ax = plt.subplots(1,3)

X,RR = np.meshgrid(x/1000.,rr)

cax = ax[0].contourf(X,RR,Psi1.T*1.e-6,10,cmap=plt.get_cmap('Reds'))
cbar = plt.colorbar(cax,ax=ax[0])
cbar.ax.set_title('[Sv]')

cax = ax[1].contourf(X,RR,Psi0.T*1.e-6,10,cmap=plt.get_cmap('Purples'))
cbar = plt.colorbar(cax,ax=ax[1])
cbar.ax.set_title('[Sv]')

cax = ax[2].contourf(X,RR,Psi2.T*1.e-6,10,cmap=plt.get_cmap('Blues'))
cbar = plt.colorbar(cax,ax=ax[2])
cbar.ax.set_title('[Sv]')

for AX in ax[1:]:
   AX.set_yticklabels([])

for AX in ax:
   AX.set_xlabel('Distance $x$ [km]')

ax[0].set_ylabel(r'Runoff $R$ [m$^2$s$^{-1}$]')
ax[0].set_title(r'AW transport')
ax[1].set_title(r'DW transport')
ax[2].set_title(r'PW transport')

for d,dd in enumerate([r'$\textbf{a)}$',r'$\textbf{b)}$',r'$\textbf{c)}$']):
   ax[d].text(0,1.05,dd,transform=ax[d].transAxes)

saveshow('construnoff2')

os.system('./copy_figs')
