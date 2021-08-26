from utils import *
from scipy.integrate import odeint
import os
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import sys


t0 = .3
dd = 5.

hstep = .1
tstep = .01
hmin = 0.
hmax = dd
tmin = 0.
tmax = 1.-t0
#t == t1-t0

h = np.arange(hmin+hstep/2.,hmax-hstep/2.,hstep)
t = np.arange(tmin+tstep/2.,tmax+tstep*3/2.,tstep)

t1 = t+t0

T,H = np.meshgrid(t,h)
T1 = T+t0

psi = T*H**2.
ht = T**2.*H**2.

#Prepare plot
prettyplot(12,4)
fig,ax = plt.subplots(1,3)

cmap = plt.get_cmap('gist_stern_r')
cmap2 = plt.get_cmap('inferno_r')


for d in [0,1,2]:
#   ax[d].contourf(t1,h,psi**.3,cmap=plt.get_cmap('Reds'),alpha=.7,zorder=1)
#   ax[d].contour(t1,h,ht**.3,5,colors='b',linestyles='solid')


#   ax[d].contour(t1,h,psi**.3,7,colors='r',linestyles='solid')
#   ax[d].contour(t1,h,ht**.3,7,colors='b',linestyles='solid')

   if d==0:
      dt = 0.*T
      dh = -H/(dd-H)
   elif d==1:
      dt = 0.*T
      dh = -1/(T*H)
   elif d==2:
      dt = -T1/(T*H**2.)
      dh = T1/(2.*T**2.*H)

   mag = np.log(dt**2+dh**2)
   lw = 1.5+3*mag/mag.max()

#   ax[d].streamplot(t1,h,dt,dh,color='k',linewidth=lw,density=[.3,.5],arrowsize=2.)
   ax[d].streamplot(t1,h,dt,dh,color=lw,cmap=cmap,density=.8,arrowsize=2.)

for AX in ax:
   AX.set_xlim([0,1.])
   AX.set_ylim([hmin,hmax])
   AX.set_xticks([0,t0,1.])
   AX.set_yticks([hmin,1,hmax])

   AX.set_xticklabels([0,'$t_0$',1])
   AX.set_yticklabels([0,1,'$\delta$'])

   AX.set_xlabel('AW temperature')
   AX.invert_yaxis()
ax[0].set_ylabel('Thermocline depth')

ax[0].set_title('Cross-thermocline diffusion')
ax[1].set_title('Cross-halocline diffusion')
ax[2].set_title('Surface heat loss')

saveshow('AWproc')

