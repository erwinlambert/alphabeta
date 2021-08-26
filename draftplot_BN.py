from utils2 import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys

#Define input parameters
dx = 10.
L = 5.e6
x = np.arange(0.,L,dx)


FSmean = np.arange(0.,.04,.005)
Rmean = np.arange(.01,.04,.005)
width = 5.e5
xmax = 0.

Psi1 = np.zeros((len(Rmean),len(FSmean)))
Psi2 = np.zeros((len(Rmean),len(FSmean)))
Psi0 = np.zeros((len(Rmean),len(FSmean)))

for i,ii in enumerate(Rmean):
   for j,jj in enumerate(FSmean):
      print ii,jj
      R = np.exp(-((x-xmax)/width)**2.)
      R = R*jj/np.mean(R) + ii

      if R.min()>0:
         #Solve equations
         Hsg,S2,Msg = solvePW(x,R)
         Hth,T1,FTint,T0 = solveAW(x,Msg)

         #Get transports
         psi1,psi2,psi0 = transport(Hth,Hsg,T1,S2,T0)
         Psi1[i,j] = psi1[0]
         Psi0[i,j] = psi0[-1]
         Psi2[i,j] = psi2[-1]

      else:
         Psi1[i,j] = np.nan
         Psi0[i,j] = np.nan
         Psi2[i,j] = np.nan

#Get plotcolors
red,blu,pur = getcols()

#Make variable plot
prettyplot(10,3)
mpl.rcParams['figure.subplot.left']   = .07
mpl.rcParams['figure.subplot.right']  = .95

fig,ax = plt.subplots(1,3)

I,J = np.meshgrid(Rmean*L*1.e-3,FSmean*L*1.e-3)

cax = ax[0].contourf(I,J,Psi1.T*1.e-6,30,cmap=plt.get_cmap('Reds'))
plt.colorbar(cax,ax=ax[0])

cax = ax[1].contourf(I,J,Psi0.T*1.e-6,30,cmap=plt.get_cmap('Purples'))
plt.colorbar(cax,ax=ax[1])

cax = ax[2].contourf(I,J,Psi2.T*1.e-6,30,cmap=plt.get_cmap('Blues'))
plt.colorbar(cax,ax=ax[2])


for AX in ax[1:]:
   AX.set_yticklabels([])

for AX in ax:
   AX.set_xlabel('Background FW flux [mSv]')

ax[0].set_ylabel(r'FW flux Baltic/North Sea [mSv]')
ax[0].set_title(r'AW inflow [Sv]')
ax[1].set_title(r'DW outflow [Sv]')
ax[2].set_title(r'PW outflow [Sv]')

for d,dd in enumerate([r'$\textbf{a)}$',r'$\textbf{b)}$',r'$\textbf{c)}$']):
   ax[d].text(0,1.05,dd,transform=ax[d].transAxes)

saveshow('rnfBN')

os.system('./copy_figs')
