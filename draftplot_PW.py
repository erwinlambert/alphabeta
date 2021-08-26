from utils2 import *
from scipy.integrate import odeint
import os
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import sys

#Choose runoff
R = .02 #m2/s

#Define H,S grid
hstep = .05
sstep = .05
hmin = 5.
hmax = 300.
smin = 32.
smax = 35.1

h = np.arange(hmin+hstep*3/2.,hmax+hstep*1/2.,hstep)
s = np.arange(smin+sstep/2.,smax+sstep*1/2.,sstep)

S,H = np.meshgrid(s,h)

#Get fields
dHdx0,dHdx1,dSdx0,dSdx1,Psi2 = PWproc(H,S,R)

#Prepare plot
prettyplot(8,4)
fig,ax = plt.subplots(1,2)

#cmap = plt.get_cmap('gist_stern_r')
red,blu,pur = getcols()
#cmap = plt.get_cmap('gist_heat_r')
cmap = plt.get_cmap('gist_earth_r')
cmap2 = plt.get_cmap('inferno_r')
#cmap = plt.get_cmap('Greys')

n0 = 6
n1 = .005
n2 = 1.
n3 = 2.
n4 = .5

#Make diffusion plot
dax = ax[0]
dax.contour(s,h,Psi2**n4,n0,colors='.5',linestyles='solid',linewidths=1.)
mag = (dHdx0**2+dSdx0**2)**n1
lw = n2+n3*mag/mag.max()
dax.streamplot(s,h,dSdx0,dHdx0,color=lw,cmap=cmap,density=.8,arrowsize=10.)

#Make runoff plot
dax = ax[1]
dax.contour(s,h,Psi2**n4,n0,colors='.5',linestyles='solid',linewidths=1.)
mag = (dHdx1**2+dSdx1**2)**n1
lw = n2+n3*mag/mag.max()
dax.streamplot(s,h,dSdx1,dHdx1,color=lw,cmap=cmap,density=.8,arrowsize=10.)

#Make up plot
for AX in ax:
   AX.set_xlim([smin,smax])
   AX.set_ylim([hmin,hmax])
   AX.set_xticks(np.arange(smin,smax,1))
   AX.set_xlabel(r'PW salinity $S_2$ [g/kg]')
   AX.invert_yaxis()

ax[1].set_yticklabels([])
ax[0].set_ylabel(r'PW thickness $H_2$ [m]')
ax[0].set_title('Cross-halocline diffusion')
ax[1].set_title('Runoff')

for d,dd in enumerate([r'$\textbf{a)}$',r'$\textbf{b)}$']):
   ax[d].text(0,1.05,dd,transform=ax[d].transAxes)


saveshow('PWproc')



#Prepare plot
prettyplot(4,4)
fig,ax = plt.subplots(1,1)

#Plot transformation field
ax.contour(s,h,Psi2**n4,n0,colors='.5',linestyles='solid',linewidths=1.)
dHdx = dHdx0+dHdx1
dSdx = dSdx0+dSdx1
mag = (dHdx**2+dSdx**2)**n1
lw = n2+n3*mag/mag.max()
ax.streamplot(s,h,dSdx,dHdx,color=lw,cmap=cmap,density=.8,arrowsize=10.)

#Integrate PW
dx = 1.
L = 5.e6
x = np.arange(0.,L,dx)
R = .02+0.*x
dots = [0,1000000,2000000,3000000,4000000,-1]
Hsg,S2,Msg = solvePW(x,R)

ax.plot(S2,Hsg,color=blu,linewidth=4.,zorder=9)
ax.scatter(S2[0],Hsg[0],150,color=blu,label='Inflow',zorder=9,clip_on=False)
ax.scatter(S2[-1],Hsg[-1],150,marker='s',color=blu,label='Outflow',zorder=9,clip_on=False)

ax.scatter(S2[dots],Hsg[dots],10,color='w',zorder=10,clip_on=False)
leg = ax.legend(loc='lower left',scatterpoints=1)
leg.set_alpha(1.)

#Make up plot
ax.set_xlim([smin,smax])
ax.set_ylim([hmin,hmax])
ax.set_xticks(np.arange(smin,smax,1))
ax.set_xlabel('PW salinity $S_2$ [g/kg]')
ax.invert_yaxis()

ax.set_ylabel('PW thickness $H_2$ [m]')

saveshow('PWtot')



#Prepare plot
prettyplot(4,3.5)
fig,ax = plt.subplots(1,1)

Feddy,Fvert = PWdiff(H,S)
cax = ax.contourf(s,h,Feddy/(Feddy+Fvert),10,cmap=plt.get_cmap('BrBG'),vmin=0.,vmax=1.)
cbar = plt.colorbar(cax)

ax.plot(S2,Hsg,color=blu,linewidth=4.,zorder=9)
ax.scatter(S2[0],Hsg[0],150,color=blu,label='Inflow',zorder=9,clip_on=False)
ax.scatter(S2[-1],Hsg[-1],150,marker='s',color=blu,label='Outflow',zorder=9,clip_on=False)

ax.scatter(S2[dots],Hsg[dots],10,color='w',zorder=10,clip_on=False)
cbar.set_ticks(np.arange(0,1.1,.1))
cbar.set_ticklabels(['vert','90','80','70','60','50/50','60','70','80','90','eddy'])


#Make up plot
ax.set_xlim([smin,smax])
ax.set_ylim([hmin,hmax])
ax.set_xticks(np.arange(smin,smax,1))
ax.set_xlabel('PW salinity $S_2$ [g/kg]')
ax.invert_yaxis()

ax.set_ylabel('PW thickness $H_2$ [m]')

saveshow('PWdiff')










os.system('./copy_figs')
