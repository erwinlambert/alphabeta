from utils2 import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

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
prettyplot(4,8)
mpl.rcParams['figure.subplot.top']    = .95
mpl.rcParams['figure.subplot.bottom'] = .1

fig,ax = plt.subplots(3,1)

x = x/1000.

ax[0].plot(x,T1,color=red,label=r'$T_1$')
ax[0].plot(x,T0+0.*x,color=pur,label=r'$T_0$')
ax[1].plot(x,S2,color=blu,label=r'$S_2$')
ax[2].plot(x,Hth,color=red,label=r'$H_1$')
ax[2].plot(x,Hsg,color=blu,label=r'$H_2$')
ax[2].invert_yaxis()

for AX in ax:
   AX.legend(ncol=2)

for AX in ax[0:2]:
   AX.set_xticklabels([])

ax[2].set_xlabel('Distance $x$ [km]')

ax[0].set_ylabel(r'Temperature [$^\circ$C]')
ax[1].set_ylabel('Salinity [g/kg]')
ax[2].set_ylabel('Layer thickness [m]')

ax[0].set_yticks([0,2,4,6,8,10])
ax[1].set_yticks([33,34,35])

for d,dd in enumerate([r'$\textbf{a)}$',r'$\textbf{b)}$',r'$\textbf{c)}$']):
   ax[d].text(0,1.05,dd,transform=ax[d].transAxes)

saveshow('reference')

#Make transport plot
Psi1,Psi2,Psi0 = transport(Hth,Hsg,T1,S2,T0)
#Msg,Mth = diaptransport(Hth,Hsg,T1,S2,T0)
Ksg,Kth = getdiffs(Hth,Hsg,T1,S2,T0)
Ftth = latheatflux(Hth,T1,T0)

print np.mean(Ksg)

prettyplot(4,8)
mpl.rcParams['figure.subplot.top']    = .95
mpl.rcParams['figure.subplot.bottom'] = .1

fig,ax = plt.subplots(3,1)

ax[0].plot(x,Psi0*1.e-6,color=pur,label=r'$\Psi_0$')
ax[0].plot(x,Psi1*1.e-6,color=red,label=r'$\Psi_1$')
ax[0].plot(x,Psi2*1.e-6,color=blu,label=r'$\Psi_2$')

ax[0].legend(ncol=3)

ax[1].plot(x,Kth,color=red,label=r'$\kappa_\theta$')
ax[1].plot(x,Ksg,color=blu,label=r'$\kappa_\sigma$')

ax[1].legend(ncol=2)

ax[2].fill_between(x,Ftth*1.e-6,color='.5',alpha=.5)
ax[2].plot(x,Ftth*1.e-6,color='.5',label=r'$F^T_{int}$')
ax[2].legend()

ax[2].set_xlabel('Distance $x$ [km]')
ax[0].set_ylabel('Volume transport [Sv]')
ax[1].set_ylabel('Diffusivity [m$^2$s$^{-1}$]')
ax[2].set_ylabel('Heat flux [TW/1000km]')

for AX in ax:
   AX.set_ylim(bottom=0.)

for AX in ax[:2]:
   AX.set_xticklabels([])


for d,dd in enumerate([r'$\textbf{a)}$',r'$\textbf{b)}$',r'$\textbf{c)}$']):
   ax[d].text(0,1.05,dd,transform=ax[d].transAxes)


saveshow('reftransport')

#Make lateral heat flux plot
#Ftth = latheatflux(Hth,T1,T0)

#prettyplot(5,4)

#fig,ax = plt.subplots(1,1)

#ax.fill_between(x,Ftth*1.e-6,color='.5',alpha=.5)
#ax.plot(x,Ftth*1.e-6,color='.5')

#ax.set_xlabel('Distance $x$ [km]')
#ax.set_ylabel('Lateral heat flux $F^T_{int}$ [TW/1000km]')

#saveshow('latheatflux')

os.system('./copy_figs')
