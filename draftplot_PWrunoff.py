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

rr = np.array([.3,1.,3.])
LL3 = 1.2
LL4 = 0.22

dx = .01
s0 = .01
h0 = .01

y0 = np.array([h0,s0])
x = np.arange(0,1,dx)

def dydx(y,x,r):
   h,s = y

   dhdx = LL3*h +LL4/(s*h**2.) -r/(2.*s**2.*h)
   dsdx = -LL3*s -LL4/(h**3.) + r/(s*h**2.)

   return np.array([dhdx,dsdx])

hh = np.arange(hmin+hstep/2.,hmax+hstep*3/2.,hstep)
ss = np.arange(smin+sstep/2.,smax+sstep*3/2.,sstep)

S,H = np.meshgrid(ss,hh)

psi = S*H**2.
fw = S**2.*H**2.

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

for d in [0,1,2]:
   r = rr[d]
#   ax[d].contourf(ss,hh,psi,15,cmap=cmap)
#   if d in [0,1]:
#   ax[d].contourf(ss,hh,fw**.3,cmap=plt.get_cmap('Blues'),alpha=1)
#   ax[d].contour(ss,hh,psi**.3,5,colors='r',linestyles='solid')
#   else:   
   ax[d].contourf(ss,hh,psi**.3,cmap=plt.get_cmap('Reds'),alpha=.7)
#   ax[d].contour(ss,hh,fw**.3,5,colors='b',linestyles='solid')


   dse = -S
   dhe = H
   dsv = -1/H**3.
   dhv = 1/(S*H**2.)
   dsr = 1./(S*H**2)
   dhr = -1./(2*S**2*H)
  
   ds = LL3*dse+LL4*dsv+r*dsr
   dh = LL3*dhe+LL4*dhv+r*dhr

   mag = np.log(ds**2+dh**2)
   lw = 1.5+3.*mag/mag.max()
   ax[d].streamplot(ss,hh,ds,dh,color='k',linewidth=lw,density=[.3,.5],arrowsize=2)

#  Compute jacobian and mark region with positive feedback!

   ax[d].set_title('Runoff = '+str(rr[d]))

   y = odeint(dydx,y0,x,args=(r,))

   h = y[:,0]
   s = y[:,1]

   col = blu
   ax[d].plot(s,h,color=col,linewidth=5.,zorder=9)
   ax[d].scatter(s[0],h[0],150,color=col,marker='o',zorder=9,clip_on=False)
   ax[d].scatter(s[-1],h[-1],150,color=col,marker='s',zorder=9)

   ax[d].scatter(s[::len(x)/9],h[::len(x)/9],5,color='1',marker='o',zorder=10,clip_on=False)

for AX in ax:
   AX.set_xlim([smin,smax])
   AX.set_ylim([hmin,hmax])
   AX.set_xticks([smin,smax])
   AX.set_yticks([hmin,hmax])

   AX.set_xlabel('Salinity contrast')
   AX.invert_yaxis()
ax[0].set_ylabel('Halocline depth')

 
fname = '../figures/draftplot_PWrunoff.png'
plt.savefig(fname,dpi=300)
os.system('eog '+fname)

