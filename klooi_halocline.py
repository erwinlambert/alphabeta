from scipy.integrate import odeint
import os
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import sys

#Nondim parameters
LL3 = .8
LL4 = .2
LL5 = .6
r = 1.

hsg0 = .1
s0 = .1

dx = .01

#arrays
y0  = np.array([hsg0,s0])
x   = np.arange(0,1,dx)
########################

#DIFFERENTIAL EQUATIONS
def dydx(y,x):
   #Current values
   hsg,s = y

   #Increments d/dx
   dhsg = LL3*hsg + LL4/(s*hsg**2) + .5*r*LL5/(s**2*hsg)
   ds = -LL3*s - LL4/(hsg**3) + r*LL5/(s*hsg**2)

   return np.array([dhsg,ds])
#######################

y = odeint(dydx,y0,x)

hsg = y[:,0]
s = y[:,1]

#Define plot
mpl.rcParams['xtick.labelsize']       = 12
mpl.rcParams['ytick.labelsize']       = 12
mpl.rcParams['lines.linewidth']       = 3.
mpl.rcParams['axes.labelsize']        = 12
mpl.rcParams['axes.titlesize']        = 12
mpl.rcParams['legend.fontsize']       = 8
mpl.rcParams['figure.subplot.hspace'] = .2
mpl.rcParams['figure.subplot.wspace'] = .2
mpl.rcParams['figure.figsize']        = 6,8
mpl.rcParams['font.family']           = 'serif'
mpl.rcParams['font.serif']            = 'palatino'
#mpl.rcParams['text.usetex']          = True
mpl.rcParams['patch.linewidth']       = 0

red = np.array([215,45,38])/255.     #Red color
blu = np.array([66,118,180])/255.    #Blue color
pur = np.array([119,43,133])/255. #Purple

fig,ax = plt.subplots(3,1)
ax[0].plot(x,hsg)
ax[0].plot(x,s)

ax[1].plot(x,s*hsg**2)
ax[1].plot(x,1/hsg)

ax[2].plot(s,hsg)

fname = '../figures/klooi_halocline.png'
plt.savefig(fname,dpi=300)
os.system('eog '+fname)

