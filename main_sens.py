import sys
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from parameters import *
from integrate import *

#Choose variable parameter
#Choices:
#Gam,Run,Hin,H2_0,S2_0,Ta,T1,a,b,dc,ds,ctheta,csigma,L,g,f,A0,Cp,rho,S1,Ht,Hs,k,ntheta,nsigma

#Range of variables and basic choices are defined in parameters.py

varkey = 'L'
Nsteps = 40

#Don't change stuff below here

savename = './Figures/Sensitivity_'+varkey+'.png'

red = np.array([215,45,38])/255.     #Red color
blu = np.array([66,118,180])/255.    #Blue color
pur = np.array([119,43,133])/255. #Purple

#Define figure
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12
mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['axes.labelsize'] = 16
mpl.rcParams['figure.subplot.hspace'] = .2
mpl.rcParams['figure.subplot.wspace'] = .4
mpl.rcParams['figure.figsize'] = 6,6
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'palatino'
mpl.rcParams['font.size'] = 12
mpl.rcParams['text.usetex'] = True
mpl.rcParams['patch.linewidth'] = 0.0

#Load basic parameters
parms = parameters()
for i in parms.keys():
    exec(i+"=parms['"+i+"']['basicval']")

#Get info of variable parameter MAKE THIS LOOP
basicval =   parms[varkey]['basicval']
rangemin =   parms[varkey]['rangemin']
rangemax =   parms[varkey]['rangemax']
conversion = parms[varkey]['conversion']
label =      parms[varkey]['label']

#Calculate
values = np.linspace(rangemin,rangemax,Nsteps)
In = []
Ov = []
Es = []

#Get values for basic parameter
#Integrate
output = integrate(parms)
for i in output.keys():
    exec(i+" = output['"+i+"']")
Inref = Psi1[0]
Ovref = Psi0[-1]
Esref = Psi2[-1]

#Calculate sensitivity 

for i,ii in enumerate(values):
    #Exchange parameter
    parms[varkey]['basicval'] = ii

    #Integrate
    output = integrate(parms)
    for i in output.keys():
        exec(i+" = output['"+i+"']")
    In = np.append(In,Psi1[0])
    Ov = np.append(Ov,Psi0[-1])
    Es = np.append(Es,Psi2[-1])

fig,ax = plt.subplots()

ax.plot(values*conversion,In,color=red,label='Inflow')
ax.plot(values*conversion,Ov,color=pur,label='Overturning')
ax.plot(values*conversion,Es,color=blu,label='Estuarine c')

ax.scatter(basicval*conversion,Inref,color=red,s=100)
ax.scatter(basicval*conversion,Ovref,color=pur,s=100)
ax.scatter(basicval*conversion,Esref,color=blu,s=100)

ax.legend(loc='upper left',prop={'size':12})
ax.set_ylabel('Transport [Sv]')
ax.set_xlabel(label)
ax.set_ylim([0.,15.])

plt.savefig(savename,dpi=300)
os.system('eog ' + savename)

