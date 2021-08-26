from scipy.integrate import odeint
import os
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import sys


hstep = .01
sstep = .01
hmin = 0.
hmax = 2.
smin = 0.
smax = 3.

rmin = .01
rmax = 3.
rstep = .1

rr = np.arange(rmin,rmax+rstep,rstep)

LL3 = 1.2
LL4 = 0.22

dx = .01
s0 = .01
h0 = .01

y0 = np.array([h0,s0])
x = np.arange(0,1+dx,dx)

PSI = np.zeros((len(rr),len(x)))
S   = PSI.copy()
Msg = PSI.copy()

def dydx(y,x,r):
   h,s = y

   dhdx = LL3*h +LL4/(s*h**2.) -r/(2.*s**2.*h)
   dsdx = -LL3*s -LL4/(h**3.) + r/(s*h**2.)

   return np.array([dhdx,dsdx])

#Define plot
mpl.rcParams['xtick.labelsize']       = 12
mpl.rcParams['ytick.labelsize']       = 12
mpl.rcParams['lines.linewidth']       = 1.
mpl.rcParams['axes.labelsize']        = 16
mpl.rcParams['axes.titlesize']        = 16
mpl.rcParams['legend.fontsize']       = 8
mpl.rcParams['figure.subplot.hspace'] = .2
mpl.rcParams['figure.subplot.wspace'] = .2
mpl.rcParams['figure.subplot.top'] = .9
mpl.rcParams['figure.subplot.bottom'] = .15
mpl.rcParams['figure.subplot.left'] = .05
mpl.rcParams['figure.subplot.right'] = .99
mpl.rcParams['figure.figsize']        = 12,4
mpl.rcParams['font.family']           ='serif'
mpl.rcParams['font.serif']            = 'palatino'
mpl.rcParams['text.usetex']          = True
mpl.rcParams['patch.linewidth']       = 0

#red = np.array([215,45,38])/255.     #Red color
blu = np.array([66,118,180])/255.    #Blue color
#pur = np.array([119,43,133])/255. #Purple

cmap = plt.get_cmap('gist_stern_r')
cmap2 = plt.get_cmap('inferno_r')

fig,ax = plt.subplots(1,3)

for i,r in enumerate(rr):
   y = odeint(dydx,y0,x,args=(r,))

   h = y[:,0]
   s = y[:,1]

   PSI[i,:] = s*h**2
   S[i,:] = s
   Msg[i,:] = LL3*s*h**2+LL4/h

N = 30
ax[0].contourf(x,rr,PSI,N,cmap=cmap)
ax[1].contourf(x,rr,S,N,cmap=cmap)
ax[2].contourf(x,rr,np.log(Msg),N,cmap=cmap)

ax[0].set_title('PW Volume transport')
ax[1].set_title('Salinity contrast')
ax[2].set_title('PW formation (log)')

for AX in ax:
   AX.set_xlim([0,1])
   AX.set_ylim([0,rmax])

   AX.set_xlabel('Distance along boundary')
ax[0].set_ylabel('Runoff')

 
fname = '../figures/draftplot_PWx.png'
plt.savefig(fname,dpi=300)
os.system('eog '+fname)

