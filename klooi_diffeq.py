from utils import *
import numpy as np
import matplotlib.pyplot as plt


dx = .001
x = np.arange(0.,1.,dx)
wc = 1.+0.*x
d = .01
D = 0.

wc = wc - D*np.exp(-(x-.5)**2./d)
msg = 0.+0.*x

hth,t1,ft,t0 = solveAW(x,wc,msg)

prettyplot()

fig,ax = plt.subplots(1,4)
ax[0].plot(x,dim(hth,'h'),'r')
ax[1].plot(x,dim(t1,'t'),'b')
ax[2].plot(x,dim((t1-t0)*hth**2.,'psi'),'g')
ax[3].plot(x,dim(wc,'wc'),'k')


ax[0].invert_yaxis()
plt.show()


