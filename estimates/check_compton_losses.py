"""
Estimates of the inverse-Compton losses of the most energetic particles in the gap (Appendix A.3 of the paper):

  (1) mean fractional energy loss <dE/E> of a particle with gamma ~ gamma_max in ONE inverse-Compton
      scattering (Klein-Nishina regime), from the scattered-photon spectrum of Blumenthal & Gould (1970),
      Eqs. (2.48)-(2.49), for an isotropic photon field (a head-on version is printed for comparison),
  (2) the Compton scattering length l_scat and the probability P = 1 - exp(-H_gap/l_scat)
      that such a particle scatters at least once while crossing the gap,
  (3) for comparison, the Compton cooling (energy-loss) length of Eq. (A11) in the paper,
      l_cool = gamma m_e c^3 / P_C, with P_C the Klein-Nishina corrected Compton power in an
      isotropic photon field (as in estimate_pulsar_quantities.py). Since a particle loses on average
      the fraction <dE/E> of its energy per scattering, l_cool ~ l_scat(iso) / <dE/E>.

Assumptions (same as in the other appendix estimates):
  * thermal photons are monoenergetic, eps = 2.7 kT/(m_e c^2), with the blackbody number density n_BB,
  * H_gap = R_pc = R_* (R_*/R_LC)^(1/2),
  * gamma_max = min(gamma_gap, gamma_rad,curv, gamma_rad,Compton) as in estimate_pulsar_quantities.py.

Scattering rate for a particle (velocity beta along +x) and photons moving at angle psi to +x:
    dN/dt = c n (1 - beta cos psi) sigma_KN(eps'),   eps' = gamma eps (1 - beta cos psi)
where sigma_KN(x) is the total Klein-Nishina cross section with x the photon energy in the
electron rest frame. It is averaged over (a) an isotropic photon field, (b) head-on photons, and
(c) an outward hemisphere (photons moving away from the star, for a particle moving outward).

Usage:  python check_compton_losses.py
"""
import numpy as np
from estimation_functions import solve_gamma_brent, f_kn_exact, sigma_KN

# cgs constants
hplanck = 6.6261e-27; c = 2.9979e10; kev2erg = 1.602176633e-9
m_e = 9.10939e-28; echarge = 4.80325e-10; sigT = 6.65246e-25
alphaf = 1/137; lamCbar = 2.4240e-10/(2*np.pi)
Bschw = 2*np.pi*m_e**2*c**3/(hplanck*echarge)

SETUPS = {  # B [G], R_* [cm], P [s], temperatures kT [keV] to check
    "MSP": dict(B=1e8, rs=1e6, P=0.002, kT=[0.1, 0.3, 1.0]),
    "WD":  dict(B=1e8, rs=1e9, P=1000., kT=[0.01, 0.03, 0.1]),
}

def gamma_max(B, rs, P, kT):
    """gamma_max = min(gamma_gap, gamma_rad (curvature), gamma_rad (Compton)); as in estimate_pulsar_quantities.py"""
    r_lc = c*P/(2*np.pi)
    beta_pc = (2*np.pi*rs/(P*c))**1.5
    r_pc = rs*(rs/r_lc)**0.5
    bratio = B/Bschw
    r_g = m_e*c**2/(echarge*B)
    rad_curv = rs
    g_gap = echarge*beta_pc*B*r_pc/(2*m_e*c**2)
    g_curv = g_gap**0.25*(r_g/rad_curv)**-0.5*bratio**-0.5*(1.5*lamCbar/r_pc/alphaf)**0.25
    n = 16*np.pi*1.202*(kev2erg*kT/(hplanck*c))**3
    U = n*2.7*kT*kev2erg
    g_comp_thom = np.sqrt(3*m_e*c**2*g_gap/(4*sigT*U*r_pc))
    g_comp = solve_gamma_brent(4*2.7*kT/511.0, g_comp_thom**2)
    return float(min(g_gap, g_curv, g_comp)), r_pc, n

def cooling_length(g, eps, n):
    """Eq. (A11): l_cool = (3/4) m_e c^2 / (gamma sigma_T U_BB f_KN(gamma eps)), with U_BB = n eps m_e c^2
    and f_KN the exact Klein-Nishina correction (Moderski et al. 2005) with b = 4 gamma eps."""
    U = n*eps*m_e*c**2
    return 0.75*m_e*c**2/(g*sigT*U*f_kn_exact(g, 4*eps))

def scattering_length(g, eps, n, mu):
    """1/l = n sigma_T < (1 - beta mu) sigma_KN(eps') >, averaged uniformly over the cosines mu."""
    b = np.sqrt(1 - 1/g**2)
    return 1.0/(n*sigT*np.mean((1 - b*mu)*sigma_KN(g*eps*(1 - b*mu))))

def mean_fractional_loss_BG(g, eps, n_E=2000001):
    """Mean fraction of the particle energy given to the photon in one scattering, for an isotropic,
    monoenergetic photon field, from Blumenthal & Gould (1970), Eqs. (2.48)-(2.49):
        dN/dt dE1 ~ 2 q ln q + (1 + 2q)(1 - q) + (1/2) (Gamma q)^2 (1 - q)/(1 + Gamma q),
        Gamma = 4 eps gamma,  q = E1/[Gamma (1 - E1)],  1/(4 gamma^2) <= q <= 1,
    where E1 = eps1/(gamma m_e c^2) is the scattered photon energy in units of the particle energy.
    Returns <E1> = int E1 dN / int dN."""
    Gam = 4.0*eps*g
    E1max = Gam/(1.0 + Gam)                          # q = 1
    E1 = np.linspace(1e-12, E1max*(1 - 1e-12), n_E)
    q = E1/(Gam*(1.0 - E1))
    F = 2*q*np.log(q) + (1 + 2*q)*(1 - q) + 0.5*(Gam*q)**2/(1 + Gam*q)*(1 - q)
    F = np.where(q >= 1.0/(4.0*g**2), F, 0.0)        # kinematic lower limit of q
    return np.trapezoid(F*E1, E1)/np.trapezoid(F, E1)

def mean_fractional_loss(g, eps):
    """Mean fraction of the particle energy given to the photon in one scattering, for head-on photons.
    In the rest frame the photon (energy eps' = gamma eps (1+beta)) arrives along -beta and scatters
    by Theta with the Klein-Nishina differential cross section; the lab energy of the scattered
    photon is eps1 = gamma eps1' (1 - beta cos Theta)."""
    b = np.sqrt(1 - 1/g**2)
    ep = g*eps*(1 + b)
    mu = np.linspace(-1, 1, 400001)                 # cos Theta in the rest frame
    e1p = ep/(1 + ep*(1 - mu))                      # Compton formula in the rest frame
    w = (e1p/ep)**2*(e1p/ep + ep/e1p - (1 - mu**2)) # KN differential cross section (unnormalized)
    e1 = g*e1p*(1 - b*mu)                           # lab energy of the scattered photon
    return np.trapezoid(w*e1, mu)/np.trapezoid(w, mu)/g

if __name__ == "__main__":
    mu_iso = np.linspace(-1, 1, 200001)
    mu_out = np.linspace(0, 1, 100001)
    for name, s in SETUPS.items():
        print(f"\n{name}:")
        print("  kT[keV]  gamma_max   <dE/E> (B&G iso, head-on)   l_scat/H (iso, head-on, outward)    P_scatter (iso, head-on, outward)"
              "   l_cool/H (Eq. A11)   l_scat(iso)/<dE/E>/H")
        for kT in s["kT"]:
            g, H, n = gamma_max(s["B"], s["rs"], s["P"], kT)
            eps = 2.7*kT/511.0
            ls = [scattering_length(g, eps, n, mu)/H for mu in (mu_iso, np.array([-1.0]), mu_out)]
            Ps = [1 - np.exp(-1/l) for l in ls]
            dE = mean_fractional_loss_BG(g, eps)      # isotropic photon field (Blumenthal & Gould 1970)
            dE_head = mean_fractional_loss(g, eps)    # head-on photons, for comparison
            lc = cooling_length(g, eps, n)/H
            print("  %6.3f  %10.3g      %5.2f  %5.2f            %6.3g %6.3g %6.3g             %5.2f %5.2f %5.2f          %8.3g             %8.3g"
                  % (kT, g, dE, dE_head, *ls, *Ps, lc, ls[0]/dE))
