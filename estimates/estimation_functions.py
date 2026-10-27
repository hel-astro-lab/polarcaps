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

    rhoc = r_curv
    bs = b_in

    lmfp_out = np.zeros_like(x_arr, dtype=float)
    lmfp_approx_out = np.zeros_like(x_arr, dtype=float)

    for i, epsy in enumerate(x_arr):

        # Footnote 5:
        # tau_BW ≈ 0.086 * alpha_f R_curv/(lambda_C x^2 b)
        #          * chi_a^3 * exp[-8/(3 chi_a)]
        A = 0.086 * alphaf * rhoc / (lamCbar * epsy**2 * bs)

        def func(chi):
            return A * chi**3 * np.exp(-8.0/(3.0*chi)) - 1.0

        chia = fsolve(func, 1.0)[0]

        lmfp = rhoc * chia / (bs * epsy)

        # Eq. (21), using W(z) ≈ ln(z) - ln[ln(z)]
        Warg = (
            alphaf * rhoc /
            (16.0 * lamCbar * epsy**2 * bs)
        )**(1.0/3.0)

        chia_approx = (8.0/9.0) / (
            np.log(Warg) - np.log(np.log(Warg))
        )

        lmfp_approx = rhoc * chia_approx / (bs * epsy)

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
       


# -------------------------------------------------
# Angle-dependent Compton and two-photon estimates
# -------------------------------------------------
def comp_eps_ic(gamma, eps_bb, psi):
    """
    Maximum (back-scattered) IC photon energy [m_e c^2] for a soft photon of
    energy eps_bb [m_e c^2] colliding with a particle of Lorentz factor gamma
    at angle psi [rad] between the particle velocity and the photon direction
    (psi = pi is head-on). Valid also in the KN regime:
        eps_IC = 2 gamma^2 eps (1 - beta cos psi) / (1 + 2 gamma eps (1 - beta cos psi)).
    For psi = pi (and beta ~ 1) this reduces to 4 gamma^2 eps / (1 + 4 gamma eps).
    """
    beta = np.sqrt(1.0 - 1.0/gamma**2)
    a = gamma*eps_bb*(1.0 - beta*np.cos(psi))
    return 2.0*gamma*a/(1.0 + 2.0*a)


def sigma_bw(s):
    """
    Breit-Wheeler two-photon pair-production cross section in units of sigma_T
    as a function of s = eps1 eps2 (1 - cos theta)/2 (threshold s = 1),
    e.g. Gould & Schreder (1967).
    """
    s = np.asarray(s, dtype=float)
    out = np.zeros_like(s)
    m = s > 1.0
    b = np.sqrt(1.0 - 1.0/s[m])
    out[m] = (3.0/16.0)*(1.0 - b**2)*((3.0 - b**4)*np.log((1.0 + b)/(1.0 - b))
                                     - 2.0*b*(2.0 - b**2))
    return out


def comp_lmfp_2phot_angle(x1, Tkev, theta):
    """
    Two-photon pair-production mean free path [cm] of a photon with energy
    x1 [m_e c^2] moving through a blackbody photon field of temperature
    Tkev [keV] whose photons all collide with it at the fixed angle theta [rad]
    (theta = pi is head-on). The target photons are taken to be monoenergetic,
    eps2 = 2.7 kT (as in comp_lmfp_2phot), with the blackbody number density n_BB:
        1/l = (1 - cos theta) n_BB sigma_T sigma_BW(eps1 eps2 (1 - cos theta)/2),
    where sigma_BW is given by sigma_bw (in units of sigma_T).
    """
    sigT = 6.65246e-25
    hplanck = 6.6261e-27
    c = 2.9979e10
    kev2erg = 1.602176633e-9
    x1, Tkev = np.broadcast_arrays(np.atleast_1d(x1).astype(float),
                                   np.atleast_1d(Tkev).astype(float))
    Theta = Tkev/511.0
    kTperhc = kev2erg*Tkev/(hplanck*c)
    nBB = 16.0*np.pi*1.202*kTperhc**3
    omc = 1.0 - np.cos(theta)
    s = x1*2.7*Theta*omc/2.0
    rate = omc*nBB*sigT*sigma_bw(s)
    with np.errstate(divide="ignore"):
        return np.where(rate > 0.0, 1.0/np.where(rate > 0.0, rate, 1.0), np.inf)


def sigma_KN(x):
    """Total Klein-Nishina cross section in units of sigma_T; x = photon energy in the electron
    rest frame [m_e c^2]."""
    x = np.asarray(x, dtype=float)
    return 0.75*(((1+x)/x**3)*(2*x*(1+x)/(1+2*x) - np.log1p(2*x))
                 + np.log1p(2*x)/(2*x) - (1+3*x)/(1+2*x)**2)


def comp_lmfp_scatter_iso(gamma, eps, nBB, nmu=2001):
    """Compton scattering length [cm] (mean distance between scatterings) of a particle with Lorentz
    factor gamma in an isotropic, monoenergetic photon field (energy eps [m_e c^2], density nBB):
        1/l = nBB sigma_T < (1 - beta mu) sigma_KN(eps') >_mu ,  eps' = gamma eps (1 - beta mu).
    gamma, eps and nBB may be arrays of equal shape."""
    sigT = 6.65246e-25
    gamma, eps, nBB = np.broadcast_arrays(np.asarray(gamma, float), np.asarray(eps, float), np.asarray(nBB, float))
    mu = np.linspace(-1.0, 1.0, nmu)
    beta = np.sqrt(1.0 - 1.0/gamma**2)[..., None]
    fac = (1.0 - beta*mu)
    avg = np.trapezoid(fac*sigma_KN(gamma[..., None]*eps[..., None]*fac), mu, axis=-1)/2.0
    return 1.0/(nBB*sigT*avg)
