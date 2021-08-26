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

mmsg = np.array([0.3,1.,3.])

LL1 = 1.08
LL5 = 1.06
wc = 1.
gamma = 3.8

dx = .001
t10 = 1.
t00 = .4
h0 = 3.

terr = .001
Nit = 30

y0 = np.array([h0,t10,0.])
x = np.arange(0,1,dx)

def dydx(y,x,t0,msg):
   h,t1,ftth = y

   dhdx = -LL1*h/(2*wc*(dd-h)) -msg/(2*(t1-t0)*h) +wc*LL5*t1/(2*(t1-t0)**2*h)
   dtdx = -wc*LL5*t1/((t1-t0)*h**2)

   ftth = LL1*(t1-t0)**2*h**2/(wc*(dd-h))

   return np.array([dhdx,dtdx,ftth])



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

red = np.array([215,45,38])/255.     #Red color
blu = np.array([66,118,180])/255.    #Blue color
pur = np.array([119,43,133])/255. #Purple

cmap = plt.get_cmap('gist_stern_r')
cmap2 = plt.get_cmap('inferno_r')

fig,ax = plt.subplots(1,3)


for d in [0,1,2]:
   t0 = t00
   msg = mmsg[d]
   for n in range(0,Nit):
      y = odeint(dydx,y0,x,args=(t0,msg))
      t0new = y[-1,-1]/gamma
      t0 = (t0+t0new)/2.
   if abs(t0-t0new) > terr:
      sys.exit('No convergence in t0')
   else:
      print t0

   tmax = 1.-t0

   hh = np.arange(hmin+hstep/2.,hmax-hstep/2.,hstep)
   tt = np.arange(tmin+tstep/2.,tmax+tstep*3/2.,tstep)

   tt1 = tt+t0

   T,H = np.meshgrid(tt,hh)
   T1 = T+t0

   psi = T*H**2.
   ht = T**2.*H**2.

   ax[d].contourf(tt1,hh,psi**.3,cmap=plt.get_cmap('Reds'),alpha=.7,zorder=1)
#   ax[d].contourf(tt1,hh,psi**.3,30,cmap=cmap,alpha=.7,zorder=1)
#   ax[d].contour(tt1,hh,ht**.3,5,colors='b',linestyles='solid')

#   ax[d].contour(tt1,hh,psi**.3,7,colors='r',linestyles='solid')
#   ax[d].contour(tt1,hh,ht**.3,7,colors='b',linestyles='solid')

   dh = -LL1*H/(2*wc*(dd-H)) -msg/(2*T*H)  +wc*LL5*T1/(2*T**2*H)
   dt = -wc*LL5*T1/(T*H**2)

   mag = np.log(dt**2+dh**2)
   lw = 1.5+3*mag/mag.max()

   ax[d].streamplot(tt1,hh,dt,dh,color='k',linewidth=lw,density=[.3,.5],arrowsize=2)

   ax[d].set_title('PW formation = '+str(msg))

   h = y[:,0]
   t1 = y[:,1]
   t = t1-t0

   col = blu
   ax[d].plot(t1,h,color=col,linewidth=5.,zorder=9)
   ax[d].scatter(t1[0],h[0],150,color=col,marker='o',zorder=9,clip_on=False)
   ax[d].scatter(t1[-1],h[-1],150,color=col,marker='s',zorder=9)

   ax[d].scatter(t1[::len(x)/9],h[::len(x)/9],5,color='1',marker='o',zorder=10,clip_on=False)

   ax[d].set_xticks([0,t0,1])
   ax[d].set_xticklabels([0,'$t_0=$'+str(.01*int(100*t0)),1])

for AX in ax:
   AX.set_xlim([0,1.])
   AX.set_ylim([0,hmax])
   AX.set_yticks([hmin,1,hmax])

   AX.set_yticklabels([0,1,'$\delta$'])

   AX.set_xlabel('AW temperature')
   AX.invert_yaxis()
ax[0].set_ylabel('Thermocline depth')

fname = '../figures/draftplot_APconst.png'
plt.savefig(fname,dpi=300)
os.system('eog '+fname)

