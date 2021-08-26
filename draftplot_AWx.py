from scipy.integrate import odeint
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

wcmin = .1
wcmax = 3.
wcstep = .02

wcc = np.arange(wcmin,wcmax+wcstep,wcstep)

LL1 = 1.08
LL5 = 1.06
gamma = 3.8

dx = .001
t10 = 1.
t00 = .4
h0 = 3.

terr = .001
Nit = 30

y0 = np.array([h0,t10,0.])
x = np.arange(0,1+dx,dx)

PSI = np.zeros((len(wcc),len(x)))
HT  = PSI.copy()
FT  = PSI.copy()
T1  = PSI.copy()
T0  = np.zeros((len(wcc)))

def dydx(y,x,t0,wc):
   h,t1,ftth = y

   dhdx = -LL1*h/(2*wc*(dd-h)) +wc*LL5*t1/(2*(t1-t0)**2*h)
   dtdx = -wc*LL5*t1/((t1-t0)*h**2)

   ftth = LL1*(t1-t0)**2*h**2/(wc*(dd-h))

   return np.array([dhdx,dtdx,ftth])



#Define plot
mpl.rcParams['xtick.labelsize']       = 12
mpl.rcParams['ytick.labelsize']       = 12
mpl.rcParams['lines.linewidth']       = 3.
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

red = np.array([215,45,38])/255.     #Red color
blu = np.array([66,118,180])/255.    #Blue color
pur = np.array([119,43,133])/255. #Purple

cmap = plt.get_cmap('gist_stern_r')
cmap2 = plt.get_cmap('inferno_r')

fig,ax = plt.subplots(1,4)

for i,wc in enumerate(wcc):
   t0 = t00
   for n in range(0,Nit):
      y = odeint(dydx,y0,x,args=(t0,wc))
      t0new = y[-1,-1]/gamma
      t0 = (t0+t0new)/2.
   if abs(t0-t0new) > terr:
      sys.exit('No convergence in t0')
   else:
      print t0

   h = y[:,0]
   t1 = y[:,1]
   t = t1-t0

   PSI[i,:] = t*h**2
   HT[i,:]  = t**2*h**2
   FT[i,:]  = LL1*t**2*h**2/(2*wc*(dd-h)) 
   T1[i,:]  = t1
   T0[i]    = t0

N = 30
ax[0].contourf(x,wcc,PSI,N,cmap=cmap)
#ax[1].contourf(x,wcc,HT,N,cmap=cmap)
ax[1].contourf(x,wcc,T1,N,cmap=cmap)
ax[2].contourf(x,wcc,FT,N,cmap=cmap)
ax[3].plot(T0,wcc,'k')

ax[0].set_title('AW Volume transport')
ax[1].set_title('AW temperature')
ax[2].set_title('Diffusive heat loss')
#ax[3].set_title('Interior temperature')

for AX in ax:
   AX.set_xlim([0,1.])
   AX.set_ylim([0,wcmax])
#   AX.set_yticks([hmin,1,hmax])

#   AX.set_yticklabels([0,1,'$\delta$'])

   AX.set_xlabel('Distance from inflow')
ax[3].set_xlabel('Interior temperature')
ax[0].set_ylabel('Slope width')

fname = '../figures/draftplot_AWx.png'
plt.savefig(fname,dpi=300)
os.system('eog '+fname)

