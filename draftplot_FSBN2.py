from utils2 import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys

#Define input parameters
dx = 10.
L = 5.e6
x = np.arange(0.,L,dx)

#step = .005
step = .0005

FSmean = np.logspace(np.log10(.0005),np.log10(.02),num=10)
BNmean = np.logspace(np.log10(.0005),np.log10(.01),num=10)


#FSmean = np.log(np.arange(.0005,.02+step,step))
#BNmean = np.log(np.arange(.0005,.01+step,step))
rmean = .004
width = 5.e5
width2 = 2.5e5
xmax = 3.*L/5.

Psi1 = np.zeros((len(BNmean),len(FSmean)))
Psi2 = np.zeros((len(BNmean),len(FSmean)))
Psi0 = np.zeros((len(BNmean),len(FSmean)))

for i,ii in enumerate(BNmean):
   for j,jj in enumerate(FSmean):
      print ii,jj
      Rfs = np.exp(-((x-xmax)/width2)**2.)
      Rfs = Rfs*jj/np.mean(Rfs)
      Rbn = np.exp(-((x-0.)/width)**2.) 
      Rbn = Rbn*ii/np.mean(Rbn)

      R = rmean+Rfs+Rbn

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

I,J = np.meshgrid(BNmean*L*1.e-3,FSmean*L*1.e-3)

cax = ax[0].contourf(I,J,Psi1.T*1.e-6,30,cmap=plt.get_cmap('Reds'))
cbar = plt.colorbar(cax,ax=ax[0])
cbar.ax.set_title('[Sv]')

cax = ax[1].contourf(I,J,Psi0.T*1.e-6,30,cmap=plt.get_cmap('Purples'))
cbar = plt.colorbar(cax,ax=ax[1])
cbar.ax.set_title('[Sv]')

cax = ax[2].contourf(I,J,Psi2.T*1.e-6,30,cmap=plt.get_cmap('Blues'))
cbar = plt.colorbar(cax,ax=ax[2])
cbar.ax.set_title('[Sv]')

for AX in ax:
   AX.scatter(20.,65.,60,marker='X',color='.5')
   AX.set_xscale('log')
   AX.set_yscale('log')

Apd = 65.
Bpd = 20.
apd = np.argmin(FSmean*L*1.e-3-Apd)
bpd = np.argmin(BNmean*L*1.e-3-Bpd)

dp1da = (Psi1[bpd,apd+1]-Psi1[bpd,apd])*1.e-6
dp1db = (Psi1[bpd+1,apd]-Psi1[bpd,apd])*1.e-6

dp0da = (Psi0[bpd,apd+1]-Psi0[bpd,apd])*1.e-6
dp0db = (Psi0[bpd+1,apd]-Psi0[bpd,apd])*1.e-6

dp2da = (Psi2[bpd,apd+1]-Psi2[bpd,apd])*1.e-6
dp2db = (Psi2[bpd+1,apd]-Psi2[bpd,apd])*1.e-6

print 'psi1',dp1db/dp1da
print 'psi0',dp0db/dp0da
print 'psi2',dp2db/dp2da

for AX in ax[1:]:
   AX.set_yticklabels([])

for AX in ax:
   AX.set_xlabel('Baltic FW outflow [mSv]')

ax[0].set_ylabel(r'Arctic FW outflow [mSv]')
ax[0].set_title(r'AW inflow')
ax[1].set_title(r'DW outflow')
ax[2].set_title(r'PW outflow')

for d,dd in enumerate([r'$\textbf{a)}$',r'$\textbf{b)}$',r'$\textbf{c)}$']):
   ax[d].text(0,1.05,dd,transform=ax[d].transAxes)

saveshow('rnfFSBN')

os.system('./copy_figs')
