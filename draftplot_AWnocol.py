from utils2 import *
from scipy.integrate import odeint
import os
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import sys

#Choose Msg
Msg = .4

#Get T0
dx = 1.
L = 5.e6
x = np.arange(0.,L,dx)
R = .02+0.*x
dots = [0,1000000,2000000,3000000,4000000,-1]
Hsg,S2,MMsg = solvePW(x,R)
Hth,T1,FTint,T0 = solveAW(x,MMsg)

#Define H,S grid
hstep = .05
tstep = .05
hmin = 5.
hmax = 800.
tmin = T0
tmin2 = -2.
tmax = 8.

h = np.arange(hmin+hstep*3/2.,hmax+hstep*1/2.,hstep)
t = np.arange(tmin+tstep/2.,tmax+tstep*1/2.,tstep)

T,H = np.meshgrid(t,h)

#Get fields
dHdx0,dHdx1,dHdx2,dTdx0,dTdx1,dTdx2,Psi1 = AWproc(H,T,T0,Msg)

#Prepare plot
prettyplot(10,3)
fig,ax = plt.subplots(1,3)

#cmap = plt.get_cmap('gist_stern_r')
red,blu,pur = getcols()
#cmap = plt.get_cmap('gist_heat_r')
#cmap = plt.get_cmap('gist_earth_r')
#cmap2 = plt.get_cmap('inferno_r')

cmap = plt.get_cmap('Reds')

n0 = 6
n1 = .005
n2 = 1.
n3 = 2.
n4 = .5
n5 = 4.

#Make theta-diffusion plot
dax = ax[0]
#dax.contour(t,h,Psi1**n4,n0,colors='.5',linestyles='solid',linewidths=1.)
dax.contourf(t,h,Psi1*1.e-6,cmap=cmap)
mag = (dHdx0**2+dTdx0**2)**n1
lw = n2+n3*mag/mag.max()
dax.streamplot(t,h,dTdx0,dHdx0,color='k',density=.4,arrowsize=1.)
#dax.plot(T0+0.*h,h,color=pur,linewidth=n5,label=r'$T_0 = $'+str(.1*int(10.*T0)))
#dax.legend(loc='lower left')

#Make sigma-diffusion plot
dax = ax[1]
#dax.contour(t,h,Psi1**n4,n0,colors='.5',linestyles='solid',linewidths=1.)
dax.contourf(t,h,Psi1*1.e-6,cmap=cmap)
mag = (dHdx1**2+dTdx1**2)**n1
lw = n2+n3*mag/mag.max()
dax.streamplot(t,h,dTdx1,dHdx1,color='k',density=.4,arrowsize=1.)
#dax.plot(T0+0.*h,h,color=pur,linewidth=n5)

#Make surf heat loss plot
dax = ax[2]
#dax.contour(t,h,Psi1**n4,n0,colors='.5',linestyles='solid',linewidths=1.)
im = dax.contourf(t,h,Psi1*1.e-6,cmap=cmap)
mag = (dHdx2**2+dTdx2**2)**n1
lw = n2+n3*mag/mag.max()
dax.streamplot(t,h,dTdx2,dHdx2,color='k',density=.8,arrowsize=1.)
#dax.plot(T0+0.*h,h,color=pur,linewidth=n5)

fig.subplots_adjust(right=.9)
cax = fig.add_axes([.95,.2,.02,.7])
cbar = plt.colorbar(im,cax=cax)
cbar.ax.set_title('[Sv]')

#Make up plot
for AX in ax:
   AX.set_xlim([tmin2,tmax])
   AX.set_ylim([hmin,hmax])
   AX.set_xticks(np.arange(tmin2,tmax+1,1))
   AX.set_xlabel(r'AW temperature $T_1$ [$^\circ$C]')
   AX.invert_yaxis()

ax[1].set_yticklabels([])
ax[2].set_yticklabels([])
ax[0].set_ylabel(r'AW thickness $H_1$ [m]')
ax[0].set_title(r'Cross-thermocline diffusion')
ax[1].set_title(r'Cross-halocline diffusion')
ax[2].set_title(r'Surface heat loss')


for d,dd in enumerate([r'$\textbf{a)}$',r'$\textbf{b)}$',r'$\textbf{c)}$']):
   ax[d].text(0,1.05,dd,transform=ax[d].transAxes)


saveshow('AWproc')

#Prepare plot
prettyplot(4,4)
fig,ax = plt.subplots(1,1)

#Plot transformation field
#ax.contour(t,h,Psi1**n4,n0,colors='.5',linestyles='solid',linewidths=1.)
ax.contourf(t,h,Psi1*1.e-6,cmap=cmap)
dHdx = dHdx0+dHdx1+dHdx2
dTdx = dTdx0+dTdx1+dTdx2
mag = (dHdx**2+dTdx**2)**n1
lw = n2+n3*mag/mag.max()
ax.streamplot(t,h,dTdx,dHdx,color='k',density=.8,arrowsize=2.)

#Plot integration
ax.plot(T1,Hth,color='y',linewidth=4.,zorder=9)
ax.scatter(T1[0],Hth[0],150,color='y',label='Inflow',zorder=9,clip_on=False)
ax.scatter(T1[-1],Hth[-1],150,marker='s',color='y',label='Outflow',zorder=9,clip_on=False)
ax.scatter(T1[dots],Hth[dots],10,color='w',zorder=10,clip_on=False)
leg = ax.legend(loc='lower left',scatterpoints=1)
leg.set_alpha(1.)

#Make up plot
ax.set_xlim([tmin2,tmax])
ax.set_ylim([hmin,hmax])
ax.set_xticks(np.arange(tmin2,tmax+1,1))
ax.set_xlabel('AW temperature $T_1$ [$^\circ$C]')
ax.invert_yaxis()

ax.set_ylabel(r'AW thickness $H_1$ [m]')

saveshow('AWtot')

#Prepare plot
prettyplot(4,3.5)
fig,ax = plt.subplots(1,1)

Feddy,Fvert = AWdiff(H,T,T0)
cax = ax.contourf(t,h,Feddy/(Feddy+Fvert),20,cmap=plt.get_cmap('BrBG'),vmin=0.,vmax=1.)
cbar = plt.colorbar(cax)

ax.plot(T1,Hth,color='y',linewidth=4.,zorder=9)
ax.scatter(T1[0],Hth[0],150,color='y',label='Inflow',zorder=9,clip_on=False)
ax.scatter(T1[-1],Hth[-1],150,marker='s',color='y',label='Outflow',zorder=9,clip_on=False)
ax.scatter(T1[dots],Hth[dots],10,color='w',zorder=10,clip_on=False)
cbar.set_ticks(np.arange(0,1.1,.1))
cbar.set_ticklabels(['vert','90','80','70','60','50/50','60','70','80','90','eddy'])

#Make up plot
ax.set_xlim([tmin2,tmax])
ax.set_ylim([hmin,hmax])
ax.set_xticks(np.arange(tmin2,tmax+1,1))
ax.set_xlabel('AW temperature $T_1$ [$^\circ$C]')
ax.invert_yaxis()

ax.set_ylabel(r'AW thickness $H_1$ [m]')

saveshow('AWdiff')


#Prepare plot
prettyplot(8,4)
#mpl.rcParams['figure.subplot.wspace'] = .3

fig,ax = plt.subplots(1,2)

#Plot transformation field
dax = ax[1]
#dax.contour(t,h,Psi1**n4,n0,colors='.5',linestyles='solid',linewidths=1.)
dax.contourf(t,h,Psi1*1.e-6,cmap=cmap)
dHdx = dHdx0+dHdx1+dHdx2
dTdx = dTdx0+dTdx1+dTdx2
mag = (dHdx**2+dTdx**2)**n1
lw = n2+n3*mag/mag.max()
dax.streamplot(t,h,dTdx,dHdx,color='k',density=.8,arrowsize=2.)

#Plot integration
dax.plot(T1,Hth,color='y',linewidth=4.,zorder=9)
dax.scatter(T1[0],Hth[0],150,color='y',label='Inflow',zorder=9,clip_on=False)
dax.scatter(T1[-1],Hth[-1],150,marker='s',color='y',label='Outflow',zorder=9,clip_on=False)
dax.scatter(T1[dots],Hth[dots],10,color='w',zorder=10,clip_on=False)
leg = dax.legend(loc='lower left',scatterpoints=1)
leg.set_alpha(1.)



Hth,T1,FTint,T0 = solveAW(x,0.*x)

tmin = T0
h = np.arange(hmin+hstep*3/2.,hmax+hstep*1/2.,hstep)
t = np.arange(tmin+tstep/2.,tmax+tstep*1/2.,tstep)
T,H = np.meshgrid(t,h)
dHdx0,dHdx1,dHdx2,dTdx0,dTdx1,dTdx2,Psi1 = AWproc(H,T,T0,0.)

#Plot transformation field
dax = ax[0]
#dax.contour(t,h,Psi1**n4,n0,colors='.5',linestyles='solid',linewidths=1.)
im = dax.contourf(t,h,Psi1*1.e-6,cmap=cmap)
dHdx = dHdx0+dHdx1+dHdx2
dTdx = dTdx0+dTdx1+dTdx2
mag = (dHdx**2+dTdx**2)**n1
lw = n2+n3*mag/mag.max()
dax.streamplot(t,h,dTdx,dHdx,color='k',density=.8,arrowsize=2.)

#Plot integration
dax.plot(T1,Hth,color='y',linewidth=4.,zorder=9)
dax.scatter(T1[0],Hth[0],150,color='y',label='Inflow',zorder=9,clip_on=False)
dax.scatter(T1[-1],Hth[-1],150,marker='s',color='y',label='Outflow',zorder=9,clip_on=False)
dax.scatter(T1[dots],Hth[dots],10,color='w',zorder=10,clip_on=False)

leg = dax.legend(loc='lower left',scatterpoints=1)
leg.set_alpha(1.)

fig.subplots_adjust(right=.9)
cax = fig.add_axes([.95,.2,.02,.7])
cbar = plt.colorbar(im,cax=cax)
cbar.ax.set_title('[Sv]')

ax[0].text(-1.8,100,r'$\kappa_\sigma = 0 $')
ax[1].text(-1.8,100,r'$\kappa_\sigma = 0.4 $')

for d,dd in enumerate([r'$\textbf{a)}$',r'$\textbf{b)}$']):
   ax[d].text(0,1.05,dd,transform=ax[d].transAxes)

#Make up plot

for AX in ax:
   AX.set_xlim([tmin2,tmax])
   AX.set_ylim([hmin,hmax])
   AX.set_xticks(np.arange(tmin2,tmax+1,1))
   AX.set_xlabel('AW temperature $T_1$ [$^\circ$C]')
   AX.invert_yaxis()

ax[0].set_ylabel(r'AW thickness $H_1$ [m]')
ax[1].set_yticklabels([])

saveshow('AWtot2')




os.system('./copy_figs')
