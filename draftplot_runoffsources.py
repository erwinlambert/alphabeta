from utils2 import *
import numpy as np
import matplotlib.pyplot as plt
import sys

#Define input parameters
dx = 10.
L = 5.e6
x = np.arange(0.,L,dx)


Amean = 65./5000 #mSv/km = m2/s
Bmean = 20./5000
Rmean = 20./5000
width = 5.e5
width2 = 2.5e5
xmax = 3.*L/5.

Fr = Rmean+0.*x
Fa = np.exp(-((x-xmax)/width2)**2.)
Fa = Fa*Amean/np.mean(Fa)
Fb = np.exp(-((x-0.)/width)**2.)
Fb = Fb*Bmean/np.mean(Fb)

prettyplot(5,3)
fig,ax = plt.subplots(1,1)

x = x/1000.

ax.fill_between(x,0.*x,Fr,color='olive',label='Background')
ax.fill_between(x,Fr,Fr+Fb,color='darkcyan',label='Baltic')
ax.fill_between(x,Fr+Fb,Fr+Fb+Fa,color='peru',label='Arctic')
ax.plot(x,Fr+Fb+Fa,'k')

ax.legend(loc='upper left')

ax.set_ylabel(r'Runoff $R$ [m$^2$s$^{-1}$]')
ax.set_xlabel('Distance $x$ [km]')

ax.set_yticks([0,.04,.08,.12])

saveshow('rnfdistr')

os.system('./copy_figs')
