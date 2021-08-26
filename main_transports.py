import sys
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from parameters import *
from integrate import *

changes = ''

red = np.array([215,45,38])/255.     #Red color
blu = np.array([66,118,180])/255.    #Blue color
pur = np.array([119,43,133])/255. #Purple

#Define figure
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['axes.labelsize'] = 16
mpl.rcParams['figure.subplot.hspace'] = .2
mpl.rcParams['figure.subplot.wspace'] = .4
mpl.rcParams['figure.figsize'] = 10,10
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'palatino'
mpl.rcParams['font.size'] = 12
mpl.rcParams['text.usetex'] = True
mpl.rcParams['patch.linewidth'] = 0.0

#Read basic parameters
parms = parameters()
for i in parms.keys():
    exec(i+"=parms['"+i+"']['basicval']")

#Vary any parameters here
def change(key,newval):
    parms[key]['basicval'] = newval
#    changes = np.append(changes,'_'+key)

change('L',4.e6)
change('Hin',550)

#Integrate
output = integrate(parms)
for i in output.keys():
    exec(i+" = output['"+i+"']")

#Plot stuff
fig,ax = plt.subplots(3,2)
x = xx*1.e-3 #km
x0 = 0.*xx
dx = xx[1]-xx[0]


ax[0,0].plot(x,-Dvt*1.e6,'y',label="v'T'")
ax[0,0].plot(x,-Dwt*1.e6,'m',label="w'T'")
ax[0,0].plot(x,-Dvs*1.e6,'g',label="v'S'")
ax[0,0].plot(x,-Dws*1.e6,'k',label="w'S'")
ax[0,0].plot(x,-(Dvt+Dwt+Dvs+Dws)*1.e6,linewidth=4,color=red,label="Total")
ax[0,0].legend(loc='lower right',prop={'size':12})
ax[0,0].set_ylabel('Layer 1')
ax[0,0].set_xticks([])
ax[0,0].set_title('Volume flux [Sv/1000km]')

ax[1,0].plot(x,Dvt*1.e6,'y',label="v'T'")
ax[1,0].plot(x,Dwt*1.e6,'m',label="w'T'")
ax[1,0].plot(x,(Dvt+Dwt)*1.e6,color=pur,linewidth=4,label="Total")
ax[1,0].legend(loc='upper right',prop={'size':12})
ax[1,0].set_ylabel('Layer 0')
ax[1,0].set_xticks([])

ax[2,0].plot(x,Dvs*1.e6,'g',label="v'S'")
ax[2,0].plot(x,Dws*1.e6,'k',label="w'S'")
ax[2,0].plot(x,(Dvs+Dws)*1.e6,color=blu,linewidth=4,label="Total")
ax[2,0].legend(loc='upper right',prop={'size':12})
ax[2,0].set_ylabel('Layer 2')
ax[2,0].set_xlabel('Distance along boundary [km]')

def tot(d):
    dd = np.cumsum(d*dx)
    return dd

ax[0,1].plot(x,Psi1[0]-tot(Dvt),'y',label="v'T'")
ax[0,1].plot(x,Psi1[0]-tot(Dwt),'m',label="w'T'")
ax[0,1].plot(x,Psi1[0]-tot(Dvs),'g',label="v'S'")
ax[0,1].plot(x,Psi1[0]-tot(Dws),'k',label="w'S'")
ax[0,1].plot(x,Psi1[0]-tot(Dvt+Dwt+Dvs+Dws),linewidth=4,color=red,label="Total")
ax[0,1].legend(loc='lower left',prop={'size':12})
ax[0,1].set_title('Total transport [Sv]')
ax[0,1].set_xticks([])

ax[1,1].plot(x,tot(Dvt),'y',label="v'T'")
ax[1,1].plot(x,tot(Dwt),'m',label="w'T'")
ax[1,1].plot(x,tot(Dvt+Dwt),color=pur,linewidth=4,label="Total")
ax[1,1].legend(loc='upper left',prop={'size':12})
ax[1,1].set_xticks([])

ax[2,1].plot(x,tot(Dvs),'g',label="v'S'")
ax[2,1].plot(x,tot(Dws),'k',label="w'S'")
ax[2,1].plot(x,tot(Dvs+Dws),color=blu,linewidth=4,label="Total")
ax[2,1].legend(loc='upper left',prop={'size':12})
ax[2,1].set_xlabel('Distance along boundary [km]')

savename = './Figures/transports'+changes+'.png'

plt.savefig(savename,dpi=300)
os.system('eog ' + savename)

