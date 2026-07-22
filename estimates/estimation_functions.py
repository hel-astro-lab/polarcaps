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

def comp_lmfp_1phot(x_in, b_in, r_curv, p_in):

    x_arr = np.atleast_1d(x_in)  # ensure array

    period = p_in
    R_LC = c * period / (2.0 * np.pi)
    rhoc = r_curv
    bs = b_in

    lmfp_out = np.zeros_like(x_arr, dtype=float)
    lmfp_approx_out = np.zeros_like(x_arr, dtype=float)

    for i, epsy in enumerate(x_arr):
        A = 0.1425 * alphaf * rhoc / (lamC * epsy**2 * bs)

        def func(x):
            return A * x**(8.0/3.0) * np.exp(-8.0/(3.0*x)) - 1

        root = fsolve(func, 1.0)[0]
        chia = root

        lmfp = rhoc * chia / (bs * epsy)

        Afac = lamC * epsy**2 * bs / (alphaf * rhoc)
        chiaWa = 1.0 / (
            np.log(0.333 / Afac**(3.0/8.0))
            - np.log(np.log(0.333 / Afac**(3.0/8.0)))
        )
        lmfp_approx = (rhoc / (bs * epsy)) * chiaWa

        lmfp_out[i] = lmfp
        lmfp_approx_out[i] = lmfp_approx

    # return scalar if scalar input
    if np.isscalar(x_in):
        return lmfp_out[0], lmfp_approx_out[0]

    return lmfp_out, lmfp_approx_out


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
    
    
import numpy as np
from scipy.optimize import brentq, newton
from scipy.special import spence

def g(b):
    term1 = (0.5 * b + 6.0 + 6.0 / b) * np.log(1.0 + b)
    term2 = (11.0 / 12.0 * b**3 + 6.0 * b**2 + 9.0 * b + 4.0) / (1.0 + b)**2
    term3 = -2.0 + 2.0 * spence(1.0+b)
    return term1 - term2 + term3

# -------------------------------------------------
# KN suppression function (Moderski approximation)
# -------------------------------------------------
def f_kn(gamma, a):
    return (1.0 + a * gamma)**(-1.5)


def f_kn2(gamma, a):
    b = a*gamma
    return (9.0/(2*b**2))*(np.log(b)-11.0/6.0)

def f_kn_exact(gamma, a):  
    b = a*gamma
    return 9.0 * g(b) / b**3

# -------------------------------------------------
# Function whose root we want
# f(gamma) = gamma^2 F_KN - C
# -------------------------------------------------
def balance_eq(gamma, a, C):
    return gamma**2 * f_kn_exact(gamma, a) - C
    #return gamma**2 * f_kn(gamma, a) - C
    #return gamma**2 - C      

# -------------------------------------------------
# Robust solver (recommended)
# -------------------------------------------------
#def solve_gamma_brent(a, C, gmin=1e-6, gmax=1e12):
#def solve_gamma_brent(a, C, gmin=1e3, gmax=1e18):
#    """
#    Solve gamma^2 (1 + a gamma)^(-3/2) = C
#    using Brent's method (very robust).
#    """   
#    return brentq(balance_eq, gmin, gmax, args=(a, C))

def solve_gamma_brent(a, C, gmin=1e-1, gmax=1e18):

    a = np.asarray(a)
    C = np.asarray(C)

    # broadcast shapes
    a, C = np.broadcast_arrays(a, C)

    gamma = np.empty_like(a, dtype=float)

    for i in np.ndindex(a.shape):

        ai = a[i]
        Ci = C[i]

        g0 = gmin
        g1 = gmax

        f0 = balance_eq(g0, ai, Ci)
        f1 = balance_eq(g1, ai, Ci)

        success = True
        # expand bracket if needed
        while f0 * f1 > 0:
            g1 *= 10
            f1 = balance_eq(g1, ai, Ci)

            if g1 > 1e50:
                success = False
                break
                #raise RuntimeError("Could not bracket root")

        if success:
            gamma[i] = brentq(balance_eq, g0, g1, args=(ai, Ci))
        else:
            gamma[i] = 1e50


    return gamma    
# -------------------------------------------------
# Fast Newton solver (optional)
# -------------------------------------------------
def balance_eq_prime(gamma, a, C):
    """Derivative with respect to gamma."""
    term1 = 2 * gamma * (1 + a*gamma)**(-1.5)
    term2 = gamma**2 * (-1.5) * a * (1 + a*gamma)**(-2.5)
    return term1 + term2


def solve_gamma_newton(a, C):
    """
    Faster but needs a good initial guess.
    Thomson solution is excellent.
    """
    gamma0 = np.sqrt(C)  # Thomson guess
    return newton(balance_eq, gamma0, fprime=balance_eq_prime, args=(a, C))    
       
