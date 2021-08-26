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

rmin = .1
rmax = 3.
wcmin = .3
wcmax = 3.
step = .1

rr = np.arange(rmin,rmax+step,step)
wcc = np.arange(wcmin,wcmax+step,step)

LL1 = 1.08
LL3 = 1.2
LL4 = 0.22
LL5 = 1.06
gamma = 3.8

dx = .001
t10 = 1.
t00 = .3
h0 = 3.
hh0 = .01
s0 = .01


terr = .001
Nit = 50

y0 = np.array([h0,t10,0.])
yy0 = np.array([hh0,s0])
x = np.arange(0,1,dx)

PSIT = np.zeros((len(rr),len(wcc)))
PSI0 = PSIT.copy()
PSI1 = PSIT.copy()
PSI2 = PSIT.copy()

RR,WCC = np.meshgrid(rr,wcc)

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
mpl.rcParams['figure.figsize']        = 12,4
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

fig,ax = plt.subplots(1,3)

for i,r in enumerate(rr):
   for j,wc in enumerate(wcc):
      print r,wc

      yy = odeint(dydx2,yy0,x,args=(r,))

      hh = yy[:,0]
      s = yy[:,1]

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

      psit = (1-t0)*h0**2.
      psi1 = (t1-t0)*h**2.
      psi2 = s*hh**2.
      psi0 = psit-psi1-psi2

      PSIT[i,j] = psit 
      PSI2[i,j] = psi2[-1]
      PSI0[i,j] = psi0[-1]

      h = y[:,0]
      t1 = y[:,1]

ax[0].contourf(RR,WCC,PSIT.T,10,cmap=cmap)
ax[1].contourf(RR,WCC,PSI0.T,10,cmap=cmap)
ax[2].contourf(RR,WCC,PSI2.T,10,cmap=cmap)

ax[0].set_title('Atlantic Inflow')
ax[1].set_title('Deep Water outflow')
ax[2].set_title('Polar Water outflow')

ax[0].set_ylabel('Slope width')

for AX in ax:
   AX.set_xlabel('Runoff')
   AX.set_xlim([0,rmax])
   AX.set_ylim([0,wcmax])

fname = '../figures/draftplot_APwcr.png'
plt.savefig(fname,dpi=300)
os.system('eog '+fname)

