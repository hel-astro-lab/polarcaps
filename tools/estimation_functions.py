import numpy as np
from scipy.optimize import fsolve
import scipy.special

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


def comp_lmfp_1phot(x_in,b_in,r_curv,p_in):

    period = p_in 
    R_LC = c*period/(2.0*np.pi)
    rhoc = r_curv #curvature radius
    epsy = x_in
    bs = b_in
    A = 0.1425*alphaf*rhoc/(lamC*epsy**2*bs)
        
    def func(x):
        return A*x**(8.0/3.0)*np.exp(-8.0/(3.0*x)) - 1    
    root = fsolve(func, 1)
    chia = root[0]
    lmfp = rhoc*chia/(bs*epsy)
    Afac = lamC*epsy**2*bs/(alphaf*rhoc)

    LW = 0.333/(Afac)**(3.0/8.0)
    chiaW = 1.0 / scipy.special.lambertw(LW)
    LWa2 = (LW/(1.0+LW))*(1.0+2.0*LW)/(LW+np.exp(LW/(1.0+LW)))
    chiaWa = 1.0/(np.log(0.333/Afac**(3.0/8.0))-np.log(np.log(0.333/Afac**(3.0/8.0))))

    lmfp_approx = (rhoc/(bs*epsy))*chiaWa
    Tchi = 0.38*chia**(-1.0/3.0)*np.exp(-8.0/(3.0*chia))
    
    return lmfp, lmfp_approx


def comp_lmfp_2phot(x1,Tkev):
    sigT  = 6.65246e-25  # Thomson cross section; cm^2
    hplanck = 6.6261e-27 #erg s
    c = 2.9979e10 #cm/s
    kev2erg = 1.602176633e-9    
    mec2_inkev = 511.0
    
    kTperhc = kev2erg*Tkev/(hplanck*c)    
    nBB = 16.0*np.pi*1.202*kTperhc**3
    x2_av_ene = 2.7*Tkev/mec2_inkev
    
    xfac = x1*x2_av_ene/(np.log(x1*x2_av_ene))
    
    #This should be more accurate approx (working also outside ultra-rel limit):
    xfac_acc_inv = ((x1*x2_av_ene)**2-1.0)*np.log(x1*x2_av_ene)*np.heaviside(x1*x2_av_ene-1.0,1.0)/(x1*x2_av_ene)**3
    xfac_acc = 1.0/xfac_acc_inv
    
    lmfp_2phot = (1.0/(0.652*sigT))*xfac_acc/nBB
    
    return lmfp_2phot
