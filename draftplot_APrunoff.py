from scipy.integrate import odeint
from scipy.interpolate import UnivariateSpline
import os
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import sys


dd = 5.

hstep = .1
tstep = .01
hmin = 0.
hmax = dd
tmin = 0.

rr = np.array([.3,1.,3.])

LL1 = 1.08
LL3 = 1.2
LL4 = 0.22
LL5 = 1.06
wc = 1.
gamma = 3.8

dx = .001
t10 = 1.
t00 = .4
h0 = 3.
hh0 = .01
s0 = .01


terr = .001
Nit = 30

y0 = np.array([h0,t10,0.])
yy0 = np.array([hh0,s0])
x = np.arange(0,1,dx)
vert = np.arange(0,dd,.1)

def dydx(y,x,mmsg):
   h,t1,ftth = y

   dhdx = -LL1*h/(2*wc*(dd-h)) -mmsg(x)/(2*(t1-t0)*h) +wc*LL5*t1/(2*(t1-t0)**2*h)
   dtdx = -wc*LL5*t1/((t1-t0)*h**2)

   ftth = LL1*(t1-t0)**2*h**2/(wc*(dd-h))

   return np.array([dhdx,dtdx,ftth])

def dydx2(y,x,r):
   hh,s = y

   dhdx = LL3*hh +LL4/(s*hh**2.) -r/(2.*s**2.*hh)
   dsdx = -LL3*s -LL4/(hh**3.) + r/(s*hh**2.)

   return np.array([dhdx,dsdx])

#Define plot
mpl.rcParams['xtick.labelsize']       = 12
mpl.rcParams['ytick.labelsize']       = 12
mpl.rcParams['lines.linewidth']       = 1.
mpl.rcParams['axes.labelsize']        = 16
mpl.rcParams['axes.titlesize']        = 16
mpl.rcParams['legend.fontsize']       = 12
mpl.rcParams['figure.subplot.hspace'] = .2
mpl.rcParams['figure.subplot.wspace'] = .2
mpl.rcParams['figure.subplot.top'] = .9
mpl.rcParams['figure.subplot.bottom'] = .15
mpl.rcParams['figure.subplot.left'] = .1
mpl.rcParams['figure.subplot.right'] = .95
mpl.rcParams['figure.figsize']        = 8,4
mpl.rcParams['font.family']           ='serif'
mpl.rcParams['font.serif']            = 'palatino'
mpl.rcParams['text.usetex']          = True
mpl.rcParams['patch.linewidth']       = 0

red = np.array([215,45,38])/255.     #Red color
blu = np.array([66,118,180])/255.    #Blue color
pur = np.array([119,43,133])/255. #Purple

cmap = plt.get_cmap('gist_stern_r')
cmap2 = plt.get_cmap('inferno_r')

cols=[red,blu,pur]
col = '.5'

fig,ax = plt.subplots(1,2)

mmsg = UnivariateSpline(x,0.*x)

t0 = t00
for n in range(0,Nit):
   y = odeint(dydx,y0,x,args=(mmsg,))
   t0new = y[-1,-1]/gamma
   t0 = (t0+t0new)/2.
if abs(t0-t0new) > terr:
   sys.exit('No convergence in t0')
else:
   print t0

h = y[:,0]
t1 = y[:,1]

ax[1].plot(t1,h,color='k',linewidth=5.,zorder=1,label='r = 0')
ax[1].scatter(t1[0],h[0],100,color=col,marker='o',zorder=9,clip_on=False)
ax[1].scatter(t1[-1],h[-1],100,color=col,marker='s',zorder=9)


for d in [0,1,2]:
   r = rr[d]
   yy = odeint(dydx2,yy0,x,args=(r,))

   hh = yy[:,0]
   s = yy[:,1]

   ax[0].plot(s,hh,color=cols[d],linewidth=5.,zorder=1)
   ax[0].scatter(s[0],hh[0],100,color=col,marker='o',zorder=9,clip_on=False)
   ax[0].scatter(s[-1],hh[-1],100,color=col,marker='s',zorder=9)

   msg = LL3*s*hh**2.+LL4/hh
   mmsg = UnivariateSpline(x,msg)

   t0 = t00
   for n in range(0,Nit):
      y = odeint(dydx,y0,x,args=(mmsg,))
      t0new = y[-1,-1]/gamma
      t0 = (t0+t0new)/2.
   if abs(t0-t0new) > terr:
      sys.exit('No convergence in t0')
   else:
      print t0

   h = y[:,0]
   t1 = y[:,1]

   ax[1].plot(t1,h,color=cols[d],linewidth=5.,zorder=1,label='r = '+str(r))
   ax[1].scatter(t1[0],h[0],100,color=col,marker='o',zorder=9,clip_on=False)
   ax[1].scatter(t1[-1],h[-1],100,color=col,marker='s',zorder=9)
#   ax[1].plot(t0+0.*vert,vert,color=cols[d],linewidth=2.)

ax[1].legend()

ax[0].set_xlim([0,3])
ax[0].set_ylim([0,2])
ax[0].set_xlabel('Salinity contrast')
ax[0].set_ylabel('Halocline depth')
ax[0].set_title('Polar Water transformation')
ax[0].set_xticks([0,3])
ax[0].set_yticks([0,2])

ax[1].set_xlim([.5,1])
ax[1].set_ylim([2,3])
ax[1].set_xlabel('AW temperature')
ax[1].set_ylabel('Thermocline depth')
ax[1].set_title('Atlantic Water transformation')

ax[1].set_yticks([2,3])
ax[1].set_xticks([.5,1])

for AX in ax:
   AX.invert_yaxis()





fname = '../figures/draftplot_APrunoff.png'
plt.savefig(fname,dpi=300)
os.system('eog '+fname)

