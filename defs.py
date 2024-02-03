import numpy as np
from numpy import pi, sqrt, exp, log


EPS = 1e-7

#--------------------------------------------------    
# constants
me    = 9.1093826e-28   # electron mass g
mp    = 1.672621637e-24 # proton mass g
e0    = 4.803250e-10    # electron charge
c     = 2.99792458e10   # speed of light cm/s
r_0   = 2.817940325e-13 # Bohr radius; cm
h_pc  = 6.6260693e-27   # Plank constnat; erg*s
sigT  = 6.65245873e-25  # Thomson optical depth; cm**2
sigSB = 5.670400e-5     # Stefan-Boltzman constant; erg/(s*cm**2*K**4)
kB    = 1.3806505e-16   # erg/K
a_rad = 7.565767e-15    # radiation constnat; erg/(cm**3*K**4)
bQ    = 4.41e13         # Swinger field
lamC  = 2*pi*3.86e-11   # Compton wavelength
alphaF= 1/137           # fine structure constant



#--------------------------------------------------    
# simulation parameters

# Coppi/Stern magnetized cascades
#R = 1e15 # cm
#coulombLog = 16.5
#l_B    = 2.9 # magnetic compacntess
#l_inj  = 30 # electron injection
#l_heat = 0
#l_BB   = 6 #10.0 #2.5*l_inj # soft photon injection
#kTbb   = 1e-5 # Coppi 1992


# Coppi diffusive heating 
#R = 1e14 # cm
#coulombLog = 16.5
#l_B    = 1.0 # magnetic compacntess; just for testing
#l_inj  = 10 # electron injection
#l_heat = 100  # stochastic heating
#l_BB   = 10 # soft photon injection
#kTbb   = 15e-3/511.0 # black body temp


# Ghisellini 98 SSA 
#R = 1e13 # cm
#coulombLog = 20.0
#l_B    = 30 # magnetic compacntess
#l_inj  = 10 # electron injection
#l_heat = 0  # stochastic heating
#l_BB   = 1e-2    # soft photon injection
#kTbb   = 7e-8 # black body temp


# merger flare; archived version
#R = 1e7 # cm / seems to work better than 1e7
#H = 0.3*R # slab thickness
#coulombLog = 17.0
#l_B    = 1e5 # magnetic compacntess
#l_heat = 0.1*l_B #0.1*l_B #0.2*l_B # stochastic heating
#l_inj  = 0.01*l_B # electron injection
#l_BB   = 0 #0.1*l_inj   # soft photon injection
#kTbb   = 3e-4 # black body temp


# TODO merger flare; new version
#R = 1e7 # cm / seems to work better than 1e7
#H = 0.1*R # slab thickness
#coulombLog = 17.0
#l_B    = 1e8 # magnetic compacntess
#l_heat = 0.001*l_B # stochastic heating
#l_inj  = 0.001*l_B # electron injection
#l_BB   = 0 #0.1*l_inj   # soft photon injection
#kTbb   = 3e-4 # black body temp


# Black hole corona; Cyg X-1
R = 3e7 # cm 10 Rs 
H = 1.0*R # slab thickness
vol = R*R*H

coulombLog = 17.0
l_inj  = 30.0 # electron injection
l_heat = 0.0  # stochastic heating
l_B    = 0.0 #0.25*l_inj # magnetic compacntess
l_BB   = 0.05 #*l_inj # soft photon injection
kTbb   = 1e-5 #1e-3 # black body temp
tau_tar= 0.0 # enforced optical depth

#--------------------------------------------------
# very small wide grid
if False:
    Nz   = 32   # grid resolution
    zmin = 1e-4 # momentum grid minimum
    zmax = 1e4  # momentum grid maximum
    Nx   = 64   # grid resolution
    xmin = 1e-6 # photon grid minimum
    xmax = 1e4  # photon grid maximum
    grid_res = '64x-6_4_32z-4_4' 

# small wide grid
if False:
    Nz    = 100 # grid resolution
    zmin = 1e-4 # momentum grid minimum
    zmax = 1e4  # momentum grid maximum
    Nx    = 100 # grid resolution
    xmin = 1e-15 # photon grid minimum
    xmax = 1e4   # photon grid maximum
    grid_res = '100x-15_4_100z-4_4' 


# small narrow grid
if False:
    Nz    = 128 # grid resolution
    zmin = 1e-4 # momentum grid minimum
    zmax = 1e4  # momentum grid maximum
    Nx    = 128 # grid resolution
    xmin = 1e-6 # photon grid minimum
    xmax = 1e4   # photon grid maximum
    grid_res = '128x-6_4_128z-4_4' 


# NOTE: default small narrow grid; 
if False:
    Nz    = 128 # grid resolution
    zmin = 1e-3 # momentum grid minimum
    zmax = 1e4  # momentum grid maximum
    Nx    = 128 # grid resolution
    xmin = 1e-6 # photon grid minimum
    xmax = 1e4   # photon grid maximum
    grid_res = '128x-6_4_128z-3_4' 


# large narrow grid
if False:
    Nz    = 200 # grid resolution
    zmin = 1e-4 # momentum grid minimum
    zmax = 1e4  # momentum grid maximum
    Nx    = 200 # grid resolution
    xmin = 1e-6 # photon grid minimum
    xmax = 1e4   # photon grid maximum
    grid_res = '200x-6_4_200z-4_4'

# wide medium grid
if True:
    Nz   = 128 # grid resolution
    zmin = 1e-3 # momentum grid minimum
    zmax = 1e4  # momentum grid maximum
    Nx   = 256 # grid resolution
    xmin = 1e-9 # photon grid minimum
    xmax = 1e4   # photon grid maximum
    grid_res = '256x-9_4_128z-3_4'


#--------------------------------------------------    
# derived quantities

Ub = l_B*(me*c**2)/(sigT*R) # magnetic field ene dens
ub = Ub/(me*c**2) # mag. en. dens. in units of mc^2

Bfield = sqrt(8*pi*Ub) # magnetic field
xb = e0*Bfield/(2*pi*me*c)*h_pc/(me*c**2)

b = Bfield/bQ

print('R      ', R, np.log10(R))
print('l_B    ', l_B, np.log10(l_B))
print('Ub     ', Ub, np.log10(Ub))
print('ub     ', ub)
print('Bfield ', Bfield, np.log10(Bfield))
print('b      ', b, np.log10(b))
print('xb     ', xb, np.log10(xb))

