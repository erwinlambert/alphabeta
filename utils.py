from scipy.integrate import odeint
from scipy.interpolate import UnivariateSpline
import numpy as np
import matplotlib as mpl
import os
import matplotlib.pyplot as plt

# Parameters
L    = 5.e6   # m
Ht   = 1000.   # m
Hs   = 200.    # m
Ws   = 100.e3  # m
Wc   = 100.e3  # m
cth  = 0.006   #
csg  = 0.025   #
g    = 9.8     # m/s2
a    = 1.0e-4  # /degC
b    = 8.e-4   # /(g/kg)
Tin  = 8.      # degC
Ta   = -2      # degC
nth  = .5      # 
nsg  = .5      #
f    = 1.4e-4  # /s
k    = 1.e-4   # m2/s
G    = 20.     # W/(m2 degC)
A    = 1.5e12
Cp   = 4.2e6   # J/(m3 degC)
S1   = 35.     # g/kg
R    = 0.02    # m2/s

#Implicit parameters
dT = Tin-Ta
dH = Ht-Hs

#Nondim parameters
l1 = cth*dH*L/(Ws*Hs)
l2 = 2*f*k*Ws*L/(nth*g*a*dT*Hs**2*dH)
l3 = csg*L/Ws
l4 = 2*f*k*Ws*L/(nsg*g*a*dT*Hs**3)
l5 = 2*f*Ws*G*L/(Cp*g*a*dT*Hs**2)

epsilon = a*dT/(b*S1)
delta = Ht/Hs
gamma = 2*f*G*A/(Cp*g*a*dT*Hs**2)

#Integration parameters
t10  = 1.
t00  = .4
hth0 = 3.
hsg0 = .01
s0   = .01
terr = .001
Nit  = 30
eps  = .001

def nondim():
   return l1,l2,l3,l4,l5,epsilon,delta,gamma

def dim(v,type):
   if type == 'psi':
      vv = v*g*a*dT*Hs**2./(2.*f)*1.e-6
   elif type == 's':
      vv = S1-v*a*dT/b
   elif type == 't':
      vv = Ta+v*dT
   elif type == 'r':
      vv = v*g*a**2*dT**2*Hs**2/(2*f*b*S1*L)
   elif type == 'h':
      vv = v*Hs
   elif type == 'wc':
      vv = v*Ws*1.e-3
   elif type == 'fw':
      vv = v*g*a**2*dT**2*Hs**2/(2*f*b*S1)*1.e-3
   elif type == 'x':
      vv = v*L*1.e-3
   elif type == 'ft':
      vv = v*Cp*cth*g*a*dT**2*Hs*(Ht-Hs)/(2*f*Ws)*1.e-6
   elif type == 'dsdx':
      vv = -v*a*dT/(b*L)*1.e6
   elif type == 'dhdx':
      vv = v*Hs/L*1.e6
   elif type == 'dtdx':
      vv = v*dT/L*1.e6
   return vv

def solvePW(x,r):
   r = UnivariateSpline(x,r)
   def dydx(y,x,r):
      y[y<eps]=eps
      h,s = y
      dhdx = l3*h +l4/(s*h**2.) -r(x)/(2.*s**2.*h)
      dsdx = -l3*s -l4/(h**3.) + r(x)/(s*h**2.)
      return np.array([dhdx,dsdx])

   y0 = np.array([hsg0,s0])

   y = odeint(dydx,y0,x,args=(r,))

   hsg = y[:,0]
   s   = y[:,1]
   msg = l3*s*hsg**2.+l4/hsg
   return hsg,s,msg

def solveAW(x,wc,msg):
   msg = UnivariateSpline(x,msg)
   wc  = UnivariateSpline(x,wc)
   def dydx(y,x,wc,msg):
      h,t1,ft = y
      dhdx = -l1*h/(2*wc(x)*(delta-h)) -msg(x)/(2*(t1-t0)*h) +wc(x)*l5*t1/(2*(t1-t0)**2*h)
      dtdx = -wc(x)*l5*t1/((t1-t0)*h**2)

      ft = l1*(t1-t0)**2*h**2/(wc(x)*(delta-h))

      return np.array([dhdx,dtdx,ft])

   y0 = np.array([hth0,t10,0.]) 
   n = 0
   dt0 = 1.
   t0 = t00
   while ((n<Nit) and (dt0>terr)):
      y = odeint(dydx,y0,x,args=(wc,msg,))
      t0new = y[-1,-1]/gamma
      t0 = (t0+t0new)/2.
      dt0 = abs(t0-t0new)
   if abs(t0-t0new) > terr:
      sys.exit('No convergence in t0')
   else:
      print 'T_0 = ',dim(t0,'t'),'t_0 = ',t0

   hth = y[:,0]
   t1  = y[:,1]
   ft  = y[:,2]
   return hth,t1,ft,t0


def prettyplot(w,h):
   mpl.rcParams['xtick.labelsize']       = 12
   mpl.rcParams['ytick.labelsize']       = 12
   mpl.rcParams['lines.linewidth']       = 2.
   mpl.rcParams['axes.labelsize']        = 12
   mpl.rcParams['axes.titlesize']        = 16
   mpl.rcParams['legend.fontsize']       = 16
   mpl.rcParams['legend.loc']            = 'upper right'
   mpl.rcParams['figure.subplot.hspace'] = .2
   mpl.rcParams['figure.subplot.wspace'] = .2
   mpl.rcParams['figure.subplot.top']    = .9
   mpl.rcParams['figure.subplot.bottom'] = .2
   mpl.rcParams['figure.subplot.left']   = .15
   mpl.rcParams['figure.subplot.right']  = .99
   mpl.rcParams['figure.figsize']        = w,h
#   mpl.rcParams['font.family']           = 'sans-serif'
   mpl.rcParams['font.sans-serif']       = 'dejavu sans'
   mpl.rcParams['font.family']           ='serif'
   mpl.rcParams['font.serif']            = 'palatino'
   mpl.rcParams['text.usetex']          = True
   mpl.rcParams['patch.linewidth']       = 1

def getcols():
   red = np.array([215,45,38])/255.     #Red color
   blu = np.array([66,118,180])/255.    #Blue color
   pur = np.array([119,43,133])/255. #Purple
   return red,blu,pur

def saveshow(fname):
   fullname = '../figures/draftplot_'+fname+'.png'
   plt.savefig(fullname,dpi=300)
   os.system('eog '+fullname +' &')

