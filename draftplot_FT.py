from scipy.integrate import odeint
import os
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import sys


t0 = .3
dd = 5.
h0 = 3.

hstep = .1
tstep = .01
sstep = .01
hmin = 0.
hmax = dd
tmin = 0.
tmax = 1.-t0
#t == t1-t0

wc = 1.
LL1 = 1.08
gamma = 3.8

h = np.arange(hmin+hstep/2.,hmax-hstep/2.,hstep)
t = np.arange(tmin+tstep/2.,tmax+tstep*3/2.,tstep)

t1 = t+t0

T,H = np.meshgrid(t,h)
T1 = T+t0

FT = LL1*(T**2*H**2)/(wc*(dd-H))

psi = T*H**2.
ht = T**2.*H**2.

#Define plot
mpl.rcParams['xtick.labelsize']       = 12
mpl.rcParams['ytick.labelsize']       = 12
mpl.rcParams['lines.linewidth']       = 4.
mpl.rcParams['axes.labelsize']        = 16
mpl.rcParams['axes.titlesize']        = 16
mpl.rcParams['legend.fontsize']       = 8
mpl.rcParams['figure.subplot.hspace'] = .3
mpl.rcParams['figure.subplot.wspace'] = .4
mpl.rcParams['figure.subplot.top'] = .9
mpl.rcParams['figure.subplot.bottom'] = .15
mpl.rcParams['figure.subplot.left'] = .1
mpl.rcParams['figure.subplot.right'] = .95
mpl.rcParams['figure.figsize']        = 8,4
mpl.rcParams['font.family']           ='serif'
mpl.rcParams['font.serif']            = 'palatino'
mpl.rcParams['text.usetex']          = True
mpl.rcParams['patch.linewidth']       = 0

red = np.array([215,45,38])/255.     #Red color
blu = np.array([66,118,180])/255.    #Blue color
pur = np.array([119,43,133])/255. #Purple

cmap = plt.get_cmap('gist_stern_r')
cmap2 = plt.get_cmap('inferno_r')

fig,ax = plt.subplots(1,2)

dax = ax[0]
dax.contourf(t1,h,FT**.3,20,cmap=cmap)

x = np.arange(0,1,.01)
tA = 1-.5*(1-t0)*x**.3
tB = 1-.5*(1-t0)*x
tC = 1-.5*(1-t0)*x**3

hA = (h0*(1-t0)-1.5*x)/(tA-t0)
hB = (h0*(1-t0)-1.5*x)/(tB-t0)
hC = (h0*(1-t0)-1.5*x)/(tC-t0)

dax.plot(tA,hA,color=blu,label='A')
dax.plot(tB,hB,color=pur,label='B')
dax.plot(tC,hC,color=red,label='C')
#ax.legend()

dax.scatter(tB[0],hB[0],150,color='.4',marker='o',zorder=9,clip_on=False)
dax.scatter(tB[-1],hB[-1],150,color='.4',marker='s',zorder=9,clip_on=False)

dax.set_xlim([t0,1.])
dax.set_ylim([hmin,hmax])
dax.set_xticks([t0,1.])
dax.set_yticks([hmin,1,hmax])

dax.set_xticklabels(['$t_0$',1])
dax.set_yticklabels([0,1,'$\delta$'])

dax.set_xlabel('AW temperature')
dax.set_ylabel('thermocline depth')
dax.invert_yaxis()

dax=ax[1]
ftA = LL1*(tA-t0)**2*hA**2/(wc*(dd-hA))
ftB = LL1*(tB-t0)**2*hB**2/(wc*(dd-hB))
ftC = LL1*(tC-t0)**2*hC**2/(wc*(dd-hC))

dax.plot(x,ftA,color=blu)
dax.plot(x,ftB,color=pur)
dax.plot(x,ftC,color=red)

dax.scatter(0,ftB[0],150,color='.5',marker='o',zorder=9,clip_on=False)
dax.scatter(1,ftB[-1],150,color='.5',marker='s',zorder=9,clip_on=False)

t0A = np.mean(ftA)/gamma
t0B = np.mean(ftB)/gamma
t0C = np.mean(ftC)/gamma

dax.text(.5,2.2,'$t_0=$ '+str(.01*int(100*t0A)),color=blu,fontsize=16)
dax.text(.5,1.8,'$t_0=$ '+str(.01*int(100*t0B)),color=pur,fontsize=16)
dax.text(.5,1.4,'$t_0=$ '+str(.01*int(100*t0C)),color=red,fontsize=16)

#dax.fill_between(x,0.*x,(tA-t0)**2*hA**2/(dd-hA),color=blu,alpha=.0)
#dax.fill_between(x,0.*x,(tB-t0)**2*hB**2/(dd-hB),color=pur,alpha=.67)
#dax.fill_between(x,0.*x,(tC-t0)**2*hC**2/(dd-hC),color=red,alpha=.33)
dax.set_xlim([0,1])
dax.set_xticks([0,1])
dax.set_ylim(bottom=0)
dax.set_yticks([0,1,2,3])

dax.set_xlabel('Distance from inflow')
dax.set_ylabel('Heat loss to interior')

fname = '../figures/draftplot_FT.png'
plt.savefig(fname,dpi=300)
os.system('eog '+fname)

