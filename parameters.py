# Parameters

L    = 4.5e6   # m
Ht   = 1000.   # m
Hs   = 200.    # m
Wc   = 100.e3  # m
Ws   = 100.e3  # m
cth  = 0.006   #
csg  = 0.025   #
g    = 9.8     # m/s2
a    = 1.2e-4  # /degC
b    = 8.e-4   # /(g/kg)
Tin  = 10      # degC
Ta   = -2      # degC
nth  = .5      # 
nsg  = .5      #
f    = 1.4e-4  # /s
k    = 1.e-4   # m2/s
G    = 20.     # W/(m2 degC)
R    = 2.25e-2   # m2/s
A    = L**2/(4*3.14) # m2
c    = 4.2e6   # J/(m3 degC)
S1   = 35.     # g/kg

r    = 1.      #

#Implicit parameters
dT = Tin-Ta
dH = Ht-Hs

# Length scales
L1 = Ws*Hs/(cth*dH)
L2 = nth*g*a*dT*Hs**2*dH/(2*f*k*Ws)
L3 = Ws/csg
L4 = nsg*g*a*dT*Hs**3/(2*f*k*Ws)
L5 = c*g*a*dT*Hs**2/(2*f*Ws*G)

# Nondim parameters
epsilon = a*dT/(b*S1)
delta = Ht/Hs
gamma = 2*f*G*A/(c*g*a*dT*Hs**2)

wc = Wc/Ws

R =r* g*a**2*dT**2*Hs**2/(2*f*b*S1*L)

#Print values
print 'L1 = ',L1*1.e-3,' km'
print 'L2 = ',L2*1.e-3,' km'
print 'L3 = ',L3*1.e-3,' km'
print 'L4 = ',L4*1.e-3,' km'
print 'L5 = ',L5*1.e-3,' km'
print '-------------------'
print 'L/L1 = ',L/L1
print 'L/L2 = ',L/L2
print 'L/L3 = ',L/L3
print 'L/L4 = ',L/L4
print 'L/L5 = ',L/L5
print '-------------------'
print 'epsilon = ',epsilon
print 'delta = ',delta
print 'gamma = ',gamma
print 'wc = ',wc
print 'r = ',r
print '-------------------'
print 'Psi = ', g*a*dT*Hs**2/(2*f)*1.e-6,' Sv'
print 'R = ',R, 'm2/s'
print 'FW = ',R*L*1.e-3, 'mSv'
print '-------------------'
print 'phi1 = ', L3/L4
print 'phi2 = ', r*L3/L
print 'phi3 = ', L/L3
