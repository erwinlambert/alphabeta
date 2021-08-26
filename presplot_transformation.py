from scipy.integrate import odeint
import os
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import sys

#PHYSICAL PARAMETERS
L    = 4500.e3   #Length of boundary
Tin  = 10.       #Inflow temperature
Sin  = 35.       #Inflow salinity
Hin  = 600.      #Inflow thermocline depth
cth  = .006      #Eddy efficiency thermocline
csg  = .025      #Eddy efficiency halocline
g    = 9.81      #Gravitational acceleration
Ht   = 1000.     #Total depth
Hs   = 200.      #Shelf depth
f    = 1.4e-4    #Coriolis parameter
wc   = 100.e3     #Width continental slope
ws   = 50.e3     #Width continental shelf
Gint = 12.       #Restoring strength interior 
G1   = 12.       #Restoring strength wm 1 [Wm-2K-1]
G2   = 12.       #Restoring strength wm 2 
c    = 4.184e6   #Specific heat [Jm-3K-1]
a    = 1.2e-4    #Alpha
b    = 7.7e-4    #Beta
k    = 1.e-4     #Vertical diffusivity 
nth  = 0.9       #Slope parameter thermocline
nsg  = 0.5       #Slope parameter halocline
Ta   = -2.       #Atmospheric temperature
R    = 1.e-2     #Runoff .5e-2 or 1.e-2
Aint = .75 *L**2/(4*np.pi)
####################

#INTEGRATION PARAMETERS
dx   = 10.       #Step size
eps  = .01       #Small number
T0   = 5.        #Initial guess T0
Nit  = 30        #Max number of iterations
Terr = .001     #Limit for convergence T0
####################

#PREPARE FOR INTEGRATION
#initial conditions
hth0 = Hin
hsg0 = eps
t10  = Tin
t20  = Tin
s20  = Sin-eps
S1   = Sin
S0   = Sin

#arrays
y0  = np.array([hth0,hsg0,t10,t20,s20,0.])
x   = np.arange(0,L,dx)
########################

#DIFFERENTIAL EQUATIONS
def dydx(y,x,T0):
   #Current values
   Hth,Hsg,T1,T2,S2,FTin = y

   #Diapycnal fluxes
   mtl = cth*g*a*(T1-T0)*Hth**2*(Ht-Hs)/(2*f*wc*(Ht-Hth))
   mtv = k*wc*(Ht-Hth)/(nth*Hth*(Ht-Hs))
   msl = csg*g*b*(S1-S2)*Hsg**2/(2*f*ws)
   msv = k*ws/(nsg*Hsg)

   #Increments d/dx
   dhth = -(mtl+mtv+msl+msv)*f/(g*a*(T1-T0)*Hth) + G1*wc*f*(T1-Ta)/(c*g*a*(T1-T0)**2*Hth)
#   if Hsg > Hs:
#      dhsg = 0.
#      ds2  = -2*f*(msl+msv+R)/(g*b*Hsg**2)
#   else:
   dhsg = 2*f*(msl+msv)/(g*b*(S1-S2)*Hsg) - f*(2*S2-S1)*R/(g*b*(S1-S2)**2*Hsg)
   ds2  = -2*f*R*S2/(g*b*(S1-S2)*Hsg**2) + 2*f*(msl+msv)/(g*b*Hsg**2)
   dt1  = -2*f*G1*wc*(T1-Ta)/(c*g*a*(T1-T0)*Hth**2)
   dt2  = -2*f*G2*ws*(T2-Ta)/(c*g*b*(S1-S2)*Hsg**2) + 2*f*(msl+msv)*(T1-T2)/(g*b*(S1-S2)*Hsg**2)

   #Heat flux boundary-interior
   FTin = (mtl+mtv)*(T1-T0)

   return np.array([dhth,dhsg,dt1,dt2,ds2,FTin])
#######################

for n in range(0,Nit):
   print T0
   y = odeint(dydx,y0,x,args=(T0,))
   T0new = Ta + c/(Gint*Aint)*y[-1,-1]
   T0 = (T0+T0new)/2.
if abs(T0-T0new) > Terr:
   sys.exit('No convergence in T0')

Psi2 = g*b*(S1-y[:,4])*y[:,1]**2./(2.*f)*1.e-6
Psi1 = g*a*(y[:,2]-T0)*y[:,0]**2./(2.*f)*1.e-6
Psi0 = Psi1[0]-Psi1-Psi2

T1 = y[:,2]
T2 = y[:,3]
S2 = y[:,4]

x = x*1.e-3 #km

#Define plot
mpl.rcParams['xtick.labelsize']       = 12
mpl.rcParams['ytick.labelsize']       = 12
mpl.rcParams['lines.linewidth']       = 3.
mpl.rcParams['axes.labelsize']        = 12
mpl.rcParams['axes.titlesize']        = 12
mpl.rcParams['legend.fontsize']       = 8
mpl.rcParams['figure.subplot.hspace'] = .2
mpl.rcParams['figure.subplot.wspace'] = .2
mpl.rcParams['figure.subplot.top'] = .9
mpl.rcParams['figure.subplot.bottom'] = .15
mpl.rcParams['figure.figsize']        = 8,4
mpl.rcParams['font.family']           = 'serif'
mpl.rcParams['font.serif']            = 'palatino'
#mpl.rcParams['text.usetex']          = True
mpl.rcParams['patch.linewidth']       = 0

red = np.array([215,45,38])/255.     #Red color
blu = np.array([66,118,180])/255.    #Blue color
pur = np.array([119,43,133])/255. #Purple

fig,ax = plt.subplots(1,1)
ax.fill_between(x,0.*x,Psi0,facecolor=pur)
ax.fill_between(x,Psi0,Psi1+Psi0,facecolor=red)
ax.fill_between(x,Psi1+Psi0,Psi1+Psi0+Psi2,facecolor=blu)

ax.set_xlabel('Distance along boundary [km]')
ax.set_ylabel('Volume transport [Sv]')

ax.set_xlim([0,L*1.e-3])
ax.set_ylim([0,Psi1[0]+1.])

fname = '../figures/presplot_transformation.png'
plt.savefig(fname,dpi=300)
os.system('eog '+fname)

