import numpy as np

setup="rp"
#setup="msp"
#setup="msp_hot"
#setup="wd"


if setup=="rp":
    tkev = 0.015 #keV
    Bfield = 1e12
    r_star = 1e6
    period = 1.0
elif setup=="msp":
    tkev = 0.1
    Bfield = 1e8
    r_star = 1e6
    period = 0.002
elif setup=="msp_hot":
    tkev = 0.5
    Bfield = 1e8
    r_star = 1e6
    period = 0.002
elif setup=="wd":
    tkev = 0.07
    Bfield = 1e8
    r_star = 1e9
    period = 1000.0


hplanck = 6.6261e-27 #erg s
c = 2.9979e10 #cm/s
kev2erg = 1.602176633e-9
m_e = 9.10939e-28  # electron mass; g
echarge = 4.80325e-10  # electron charge; esu
lamC = 2.4240e-10   # Compton wavelength; cm
lamCbar = lamC/(2.0*np.pi)
alphaf = 1/137 #7.2974e-3
sigT  = 6.65246e-25  # Thomson cross section; cm^2
re = echarge**2/(m_e*c**2)
Bschw = 2.0*np.pi*m_e**2*c**3/(hplanck*echarge) #Schwinger field: ~ 4.414e13 G
bratio = Bfield/Bschw

r_lc = c*period/(2.0*np.pi)
beta_pc = (2.0*np.pi*r_star/(period*c))**(3.0/2.0)
r_pcap = r_star*(r_star/r_lc)**0.5
n_gj = beta_pc*Bfield/(4.0*np.pi*echarge*r_pcap)
skdepth = c/np.sqrt(4.0*np.pi*echarge**2*n_gj/m_e)
#rad_curv = np.sqrt(r_lc*r_star)
rad_curv = r_star #This is what is said we use in the draft
r_g = m_e*c**2/(echarge*Bfield)


print("beta: ",beta_pc)
print("r_pcap [cm]: ", r_pcap)
print("n_gj [cm^-3]: ", n_gj)
print("skin depth [cm]: ", skdepth)
print("r_pcap/skin_depth: ", r_pcap/skdepth)
print("r_g [cm]: ",r_g)
print("rad_curv [cm]: ",rad_curv)
print("r_g/rad_curv: ",r_g/rad_curv)
print("bratio: ",bratio)
print()

gam_gap = echarge*beta_pc*Bfield*r_pcap/(2*m_e*c**2)
gam_rad_syn = gam_gap**0.25*(r_g/rad_curv)**-0.5*bratio**-0.5*( 1.5*lamCbar/r_pcap/alphaf)**0.25
gam_thr_syn = (4.0*rad_curv/(3.0*bratio*r_g))**(1.0/3.0)

kTperhc = kev2erg*tkev/(hplanck*c)
nBB = 16.0*np.pi*1.202*kTperhc**3
av_ene = 2.7*tkev*kev2erg
x2 = 2.7*tkev/511.0
U_x = nBB*av_ene
gam_rad_compton = np.sqrt((3.0*m_e*c**2*gam_gap)/(4.0*sigT*U_x*r_pcap))

print("gam_gap: ", gam_gap)
print("gam_rad_syn: ", gam_rad_syn)
print("gam_thr_syn:", gam_thr_syn)
print("gam_rad_compton: ", gam_rad_compton)
print()
print("gam_rad_syn/gam_gap: ",gam_rad_syn/gam_gap)
print("gam_thr_syn/gam_gap:", gam_thr_syn/gam_gap)
print("gam_rad_compton/gam_gap: ",gam_rad_compton/gam_gap)
print()

t_esc = r_pcap/c
t_acc = m_e*c/(echarge*beta_pc*Bfield)

t_x_curv = lamCbar*rad_curv/(alphaf*bratio*c*r_g*gam_rad_syn)
lmfp_curv = t_x_curv*c

if(gam_rad_compton < gam_rad_syn):
    gam_rad_smaller = gam_rad_compton
else:
    gam_rad_smaller = gam_rad_syn

lmfp_compton = gam_rad_smaller*3.0*m_e*c**2/(4.0*sigT*gam_rad_smaller**2*U_x)
t_x_comp = lmfp_compton/c

t_p_gj = 1.0/np.sqrt(4.0*np.pi*echarge**2*n_gj/m_e)


print("t_acc [s]: ",t_acc)
print("t_x_curv [s]: ",t_x_curv)
print("t_x_comp [s]: ",t_x_comp)
print("t_esc [s]: ",t_esc)
print("t_p,co [s]: ",t_p_gj)

# characteristic synchrotron photon energy
xsyn = 1.5*bratio*(r_g/rad_curv)*gam_rad_syn**3
print()
print("xsyn [m_e x c^2]: ",xsyn)


if gam_rad_compton > gam_gap:
    gam_rad_compton = gam_gap

# characteristic one-time scattered photon energy
xcomp = x2+x2*(4.0/3.0)*gam_rad_compton**2
print("xcomp [m_e x c^2]: ",xcomp)
print()

if setup=="rp":
    x1 = xsyn
else:
    x1 = xcomp

#Next 1-photon mean free path:
from estimation_functions import comp_lmfp_1phot
lmfp_1phot, lmfp_1phot_approx = comp_lmfp_1phot(x1,bratio,rad_curv,period)
print("lmfps_1phot/Rpc: ", lmfp_1phot/r_pcap)

#Next 2-photon mean free path (between curvature and thermal photons):
from estimation_functions import comp_lmfp_2phot
lmfp_2phot = comp_lmfp_2phot(x1,tkev)
print("lmfps_2phot/Rpc: ", lmfp_2phot/r_pcap)

print("lmfp_compton/R_pc: ",lmfp_compton/r_pcap)
print("lmfp_curv/R_pc: ",lmfp_curv/r_pcap)

t_pp2 = lmfp_2phot/c
t_pp1 = lmfp_1phot/c

print("t_pp1 [s]: ",t_pp1)
print("t_pp2 [s]: ",t_pp2)

#Estimate limiting temperatures:

Const1 = (0.029*np.pi**0.5*echarge*hplanck**3*c**1.5/(2**0.5*sigT))**0.5
Const2 = (3.0*echarge**3*lamCbar*Bschw**2/(4.0*alphaf*m_e**3))**0.25*((2.0*np.pi)**3/(c**15))**(1.0/8.0)
T_gam_rad_equal = (Const1*Const2**(-1)*Bfield**(1.0/4.0)*period**(-3.0/8.0)*r_star**(3.0/8.0)*rad_curv**(-0.5))**0.5/kev2erg

#Double checking that formulas match:
#print(gam_gap, 2.0*np.pi**2*echarge*Bfield*r_star**3/(m_e*c**4*period**2))
#print(gam_rad_compton, Const1*Bfield**0.5*period**(-0.75)*r_star**0.75*(tkev*kev2erg)**-2)
#print(gam_rad_syn, Const2*Bfield**0.25*period**(-3.0/8.0)*r_star**(3.0/8.0)*rad_curv**0.5)

print("Const1:",Const1)
print("(Const1/Const2)**0.5:",(Const1/Const2)**0.5)

print("T_gam_rad_equal [keV], T_gam_rad_equal [K]: ", T_gam_rad_equal, T_gam_rad_equal*11604525.0062)

prefac1 = (3.0*m_e*c**5*hplanck**3)**0.5/(4.0*51.9*np.pi*sigT*Const1)**0.5*(2.0*np.pi/c)**(-0.25)
T_comp_lmpf1 = prefac1*Bfield**(-0.25)*r_star**(-9.0/8.0)*period**(5.0/8.0)/kev2erg

print("Const3: ",prefac1)

print("T_comp_lmpf1 [keV], T_comp_lmpf1 [K]: ", T_comp_lmpf1, T_comp_lmpf1*11604525.0062)

