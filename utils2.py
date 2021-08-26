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
cth  = 0.007   #
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
S1   = 35.2    # g/kg

#Integration parameters
T10  = Tin
T00  = 1.      # degC
Hth0 = 600.    # m
Hsg0 = 10.     # m
S20  = 35.1    # g/kg
Terr = .01     # degC
Nit  = 30
eps  = .001

def solvePW(x,R):
   R = UnivariateSpline(x,R)
   def dydx(y,x,R):
      y[y<eps]=eps
      Hsg,S2 = y

      Fssg = csg*g*b*(S1-S2)**2*Hsg**2/(2*f*Ws) + k*(S1-S2)*Ws/(nsg*Hsg)

      dHsg = 2*f*Fssg/(g*b*(S1-S2)**2*Hsg)-f*R(x)*S2/(g*b*(S1-S2)**2*Hsg)

      dS2 = 2*f*Fssg/(g*b*(S1-S2)*Hsg**2)-2*f*R(x)*S2/(g*b*(S1-S2)*Hsg**2)

      return np.array([dHsg,dS2])

   y0 = np.array([Hsg0,S20])

   y = odeint(dydx,y0,x,args=(R,))

   Hsg = y[:,0]
   S2  = y[:,1]
   if Hsg.min()<0.:
      Msg = 0.*Hsg
      print Hsg.min()
   else:
      Msg = csg*g*b*(S1-S2)*Hsg**2/(2*f*Ws) + k*Ws/(nsg*Hsg)
   return Hsg,S2,Msg

def solveAW(x,Msg):
   Msg = UnivariateSpline(x,Msg)
   def dydx(y,x,Msg):
      Hth,T1,FTint = y

      Ftth = Cp*cth*g*a*(T1-T0)**2*Hth**2*(Ht-Hs)/(2*f*Wc*(Ht-Hth)) + Cp*k*(T1-T0)*Wc*(Ht-Hth)/(nth*Hth*(Ht-Hs))

      dHth = -f*Ftth/(Cp*g*a*(T1-T0)**2*Hth)-f*Msg(x)/(g*a*(T1-T0)*Hth)+G*Wc*f*(T1-Ta)/(Cp*g*a*(T1-T0)**2*Hth)

      dT1 = -2*f*G*Wc*(T1-Ta)/(Cp*g*a*(T1-T0)*Hth**2)

      FTint = Ftth

      return np.array([dHth,dT1,FTint])

   y0 = np.array([Hth0,T10,0.]) 
   n = 0
   dT0 = T00-Ta
   T0 = T00
   while ((n<Nit) and (dT0>Terr)):
      y = odeint(dydx,y0,x,args=(Msg,))
      T0new = Ta + y[-1,2]/(G*A)
      T0 = (T0+T0new)/2.
      dT0 = abs(T0-T0new)
   if abs(T0-T0new) > Terr:
      sys.exit('No convergence in T0')
   else:
      print 'T_0 = ',T0

   Hth = y[:,0]
   T1  = y[:,1]
   FTint  = y[:,2]
   return Hth,T1,FTint,T0

def transport(Hth,Hsg,T1,S2,T0):
   Psi1 = g*a*(T1-T0)*Hth**2/(2*f)
   Psi2 = g*b*(S1-S2)*Hsg**2/(2*f)
   Psi0 = Psi1[0]-Psi1-Psi2
   return Psi1,Psi2,Psi0

def latheatflux(Hth,T1,T0):
   Ftth = Cp*cth*g*a*(T1-T0)**2*Hth**2*(Ht-Hs)/(2*f*Wc*(Ht-Hth)) + Cp*k*(T1-T0)*Wc*(Ht-Hth)/(nth*Hth*(Ht-Hs))
   return Ftth

def diaptransport(Hth,Hsg,T1,S2,T0):
   Msg = csg*g*b*(S1-S2)*Hsg**2/(2*f*Ws) + k*Ws/(nsg*Hsg)
   Mth = cth*g*a*(T1-T0)*Hth**2*(Ht-Hs)/(2*f*Wc*(Ht-Hth)) + k*Wc*(Ht-Hth)/(nth*Hth*(Ht-Hs))
   return Msg,Mth

def PWproc(Hsg,S2,R):
   Fssg = csg*g*b*(S1-S2)**2*Hsg**2/(2*f*Ws) + k*(S1-S2)*Ws/(nsg*Hsg)

   dHdx0 = 2*f*Fssg/(g*b*(S1-S2)**2*Hsg)
   dHdx1 = -f*R*S2/(g*b*(S1-S2)**2*Hsg)

   dSdx0 = 2*f*Fssg/(g*b*(S1-S2)*Hsg**2)
   dSdx1 = -2*f*R*S2/(g*b*(S1-S2)*Hsg**2)

   Psi2 = g*b*(S1-S2)*Hsg**2/(2*f)
   return dHdx0,dHdx1,dSdx0,dSdx1,Psi2

def AWproc(Hth,T1,T0,Msg):
   Ftth = Cp*cth*g*a*(T1-T0)**2*Hth**2*(Ht-Hs)/(2*f*Wc*(Ht-Hth)) + Cp*k*(T1-T0)*Wc*(Ht-Hth)/(nth*Hth*(Ht-Hs))

   dHdx0 = -f*Ftth/(Cp*g*a*(T1-T0)**2*Hth)
   dHdx1 = -f*Msg/(g*a*(T1-T0)*Hth)
   dHdx2 = G*Wc*f*(T1-Ta)/(Cp*g*a*(T1-T0)**2*Hth)

   dTdx0 = 0.*Hth
   dTdx1 = 0.*Hth
   dTdx2 = -2*f*G*Wc*(T1-Ta)/(Cp*g*a*(T1-T0)*Hth**2)

   Psi1 = g*a*(T1-T0)*Hth**2/(2*f)
   return dHdx0,dHdx1,dHdx2,dTdx0,dTdx1,dTdx2,Psi1

def PWdiff(Hsg,S2):
   Feddy = csg*g*b*(S1-S2)**2*Hsg**2/(2*f*Ws)
   Fvert = k*(S1-S2)*Ws/(nsg*Hsg)
   return Feddy,Fvert

def AWdiff(Hth,T1,T0):
   Feddy = Cp*cth*g*a*(T1-T0)**2*Hth**2*(Ht-Hs)/(2*f*Wc*(Ht-Hth))
   Fvert = Cp*k*(T1-T0)*Wc*(Ht-Hth)/(nth*Hth*(Ht-Hs))
   return Feddy,Fvert

def getdiffs(Hth,Hsg,T1,S2,T0):
   Ksg = csg*g*b*(S1-S2)*Hsg**2/(2*f*Ws) + k*Ws/(nsg*Hsg)
   Kth = cth*g*a*(T1-T0)*Hth**2*(Ht-Hs)/(2*f*Wc*(Ht-Hth)) + k*Wc*(Ht-Hth)/(nth*Hth*(Ht-Hs))
   return Ksg,Kth

def prettyplot(w,h):
   mpl.rcParams['xtick.labelsize']       = 12
   mpl.rcParams['ytick.labelsize']       = 12
   mpl.rcParams['lines.linewidth']       = 2.
   mpl.rcParams['axes.labelsize']        = 12
   mpl.rcParams['axes.titlesize']        = 12
   mpl.rcParams['legend.fontsize']       = 12
   mpl.rcParams['legend.loc']            = 'upper right'
   mpl.rcParams['figure.subplot.hspace'] = .2
   mpl.rcParams['figure.subplot.wspace'] = .2
   mpl.rcParams['figure.subplot.top']    = .9
   mpl.rcParams['figure.subplot.bottom'] = .2
   mpl.rcParams['figure.subplot.left']   = .15
   mpl.rcParams['figure.subplot.right']  = .95
   mpl.rcParams['figure.figsize']        = w,h
#   mpl.rcParams['font.family']           = 'sans-serif'
   mpl.rcParams['font.sans-serif']       = 'dejavu sans'
#   mpl.rcParams['font.family']           ='serif'
#   mpl.rcParams['font.serif']            = 'palatino'
   mpl.rcParams['font.size']             = 12
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
