from utils import *
from scipy.integrate import odeint
import os
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import sys


hstep = .01
sstep = .01
hmin = 0.
hmax = 1.
smin = 0.
smax = 3.


h = np.arange(hmin+hstep/2.,hmax+hstep*3/2.,hstep)
s = np.arange(smin+sstep/2.,smax+sstep*3/2.,sstep)

#h = dim(h,'h')
#s = dim(s,'s')

S,H = np.meshgrid(s,h)

psi = S*H**2.
fw = S**2.*H**2.

#Prepare plot
prettyplot(12,6)
fig,ax = plt.subplots(1,2)

cmap = plt.get_cmap('gist_stern_r')
cmap2 = plt.get_cmap('inferno_r')

for d in [0,1]:
#   ax[d].contourf(s,h,psi,15,cmap=cmap)
#   if d in [0,1]:
#   ax[d].contourf(s,h,fw**.3,cmap=plt.get_cmap('Blues'),alpha=1)
#   ax[d].contour(s,h,psi**.3,5,colors='r',linestyles='solid')
#   else:   
#   ax[d].contourf(s,h,psi**.3,cmap=plt.get_cmap('Reds'),alpha=.7)
#   ax[d].contour(s,h,fw**.3,5,colors='b',linestyles='solid')

   if d==0:
      ds = -S -1/(H**3.)
      dh = H + 1/(S*H**2.)
      fac = 2.
   elif d==1:
      ds = 1/(S*H**2.)
      dh = -.5/(S**2*H)
      fac = 5.

   mag = (ds**2+dh**2)**.03
   lw = 0.5+3.*mag/mag.max()
#   ax[d].streamplot(s,h,ds,dh,color='k',linewidth=lw,density=[.2,.4],arrowsize=2)
   ax[d].streamplot(s,h,ds,dh,color=lw,cmap=cmap,density=.8,arrowsize=2.)

for AX in ax:
#   AX.set_xlim([smin,smax])
#   AX.set_ylim([hmin,hmax])
#   AX.set_xticks([smin,smax])
#   AX.set_yticks([hmin,hmax])

   AX.set_xlabel('PW salinity')
   AX.invert_yaxis()
   AX.invert_xaxis()

for AX in ax[1:]:
   AX.set_yticks([])

ax[0].set_ylabel('Halocline depth')

ax[0].set_title('Cross-halocline diffusion')
ax[1].set_title('Runoff')

saveshow('PWproc')

