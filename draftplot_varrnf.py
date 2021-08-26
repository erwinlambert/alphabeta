from utils import *
import numpy as np
import matplotlib.pyplot as plt

#Define input parameters
dx = .0001
x = np.arange(0.,1.,dx)

rb = 2.
rm = 4.

wc = 1.+0.*x

#Get plotcolors
red,blu,pur = getcols()
col = [red,pur,blu]
tp = [':','-','--']
lab = ['Arctic','uniform','Baltic']

#Get nondimensional parameters
l1,l2,l3,l4,l5,epsilon,delta,gamma = nondim()

#Prepare plot
prettyplot(8,8)
fig,ax = plt.subplots(2,2)

ax[0,0].plot(dim(x,'x'),0.*x,'--k',linewidth=.5)

for i in [0,1,2]:
   if i == 0:
      r = rb+0.*x
      r2 = np.exp(-(x-.7)**2./.02)
      r2 = (rm-rb)*r2/np.mean(r2)
      r = r+r2
   elif i == 1:
      r = rm+0.*x
   elif i == 2:
      r = rb+0.*x
      r2 = np.exp(-(x-.0)**2./.02)
      r2 = (rm-rb)*r2/np.mean(r2)
      r = r+r2

   #Solve equations
   hsg,s,msg = solvePW(x,r)
   hth,t1,ft,t0 = solveAW(x,wc,msg)

   psi2 = s*hsg**2.
   psi1 = (t1-t0)*hth**2.
   psi0 = psi1[0]-psi1-psi2

   dax = ax[0,0]
   dax.plot(dim(x,'x'),dim(r,'r'),tp[i],color='b',label=lab[i])
   dax = ax[1,0]
   dax.plot(dim(x,'x'),dim(psi2,'psi'),tp[i],color=blu)
   dax = ax[0,1]
   dax.plot(dim(x,'x'),dim(psi1,'psi'),tp[i],color=red)   
   dax = ax[1,1]
   dax.plot(dim(x,'x'),dim(psi0,'psi'),tp[i],color=pur)


ax[0,0].legend(loc='upper right')
#for AX in ax[0:-1]:
#   AX.set_xticklabels([])
ax[0,0].set_ylabel('Runoff [m$^2$s$^{-1}$]')
ax[1,0].set_ylabel('PW transport [Sv]')
ax[0,1].set_ylabel('AW transport [Sv]')
ax[1,1].set_ylabel('DW transport [Sv]')
for AX in ax[1,:]:
   AX.set_xlabel('Distance [km]')

saveshow('varrnf')


#Prepare plot
prettyplot(12,4)
fig,ax = plt.subplots(1,3)

#ax[0].plot(dim(x,'x'),0.*x,'--k',linewidth=.5)

for i,rr in enumerate([1.,2.,4.]):
   r = rr+0.*x
   #Solve equations
   hsg,s,msg = solvePW(x,r)
   hth,t1,ft,t0 = solveAW(x,wc,msg)

   psi2 = s*hsg**2.
   psi1 = (t1-t0)*hth**2.
   psi0 = psi1[0]-psi1-psi2

   dax = ax[0]
   dax.plot(dim(x,'x'),dim(psi2,'psi'),tp[i],color=blu,label='R='+str(dim(rr,'r')))
   dax = ax[1]
   dax.plot(dim(x,'x'),dim(psi1,'psi'),tp[i],color=red)
   dax = ax[2]
   dax.plot(dim(x,'x'),dim(psi0,'psi'),tp[i],color=pur)


ax[0].legend(loc='upper left')
#for AX in ax[0:-1]:
#   AX.set_xticklabels([])
ax[0].set_ylabel('PW transport [Sv]')
ax[1].set_ylabel('AW transport [Sv]')
ax[2].set_ylabel('DW transport [Sv]')
for AX in ax:
   AX.set_xlabel('Distance [km]')

saveshow('constrnf')

