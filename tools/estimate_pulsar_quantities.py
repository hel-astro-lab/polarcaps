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

if setup=="wd":
    tkelvin_grid = np.logspace(np.log10(2e4), np.log10(1e7), 200)
else:
    tkelvin_grid = np.logspace(np.log10(5e4), np.log10(1e7), 200)

tgrid_plot = True
keV2Kelvin = 11604525.0062

if tgrid_plot:
    tkev = tkelvin_grid/keV2Kelvin
else:
    tkev = np.array([tkev])

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


#a = 0.0001 #4.0*x2 #1e-6
a = 4.0*x2 #0.00002 #0.0001 #4.0*x2
C = gam_rad_compton**2 #1e10
#print("a, b: ",a, a*np.sqrt(C))
print("b: ",a*np.sqrt(C))

from estimation_functions import solve_gamma_brent as solve_gamma_brent
from estimation_functions import solve_gamma_newton as solve_gamma_newton

from estimation_functions import f_kn as f_kn
from estimation_functions import f_kn2 as f_kn2
from estimation_functions import f_kn_exact as f_kn_exact

gam_rad_compton_KN = solve_gamma_brent(a, C)

gam_rad_compton = gam_rad_compton_KN

#print("gam_rad_compton_KN:", gam_rad_compton_KN)
#print(gam_rad_compton)
#print(tkev)

#exit()

#gam_rad_compton_KN = np.array([0.9e13,1e13,2e13,3e13,4e13,5e13,6e13,7e13])
#gam_rad_compton_KN = np.array([1e12,1e13,1e14,1e15,1e20,1e40,1e60,1e80,1e100])
#gam_rad_compton_KN = np.array([1e5,1e6,1e7,1e8,1e9,1e10,1e11,1e12,1e13])

#print(tkev[130])
#print(gam_rad_compton[130])
#print(C[130])
#print(a[130]*gam_rad_compton_KN)

#print("Residual:", gam_rad_compton_KN**2*f_kn(gam_rad_compton_KN,a[130])-C[130])
#print("Residual:", gam_rad_compton_KN**2*f_kn2(gam_rad_compton_KN,a[130])-C[130])
#print("Residual:", gam_rad_compton_KN**2*f_kn_exact(gam_rad_compton_KN,a[130])-C[130])

#gam_rad_compton_KN_approx = a[120]**3*gam_rad_compton[120]**4
#print("gam_rad_compton_KN_approx:", gam_rad_compton_KN_approx)

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

if setup=="rp":
    gam_smallest = np.minimum(gam_rad_compton, gam_rad_syn)
    gam_smallest = np.minimum(gam_smallest,np.array([gam_gap]*len(gam_rad_compton)))
else:
    gam_smallest = np.minimum(gam_rad_compton, np.array([gam_gap]*len(gam_rad_compton)))

#gam_smallest = np.minimum(gam_rad_compton,np.array([gam_gap]*len(gam_rad_compton)))
#gam_smallest = gam_gap

#lmfp_compton = gam_smallest*3.0*m_e*c**2/(4.0*sigT*gam_smallest**2*U_x)

#exit()

P_C_KN = 4.0*sigT*c*gam_smallest**2*U_x*f_kn_exact(gam_smallest,a)/3.0
lmfp_compton_KN = gam_smallest*m_e*c**3/P_C_KN
lmfp_compton = lmfp_compton_KN

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


#if gam_rad_compton > gam_gap:
#    gam_rad_compton = gam_gap

gam_rad_compton = np.minimum(gam_rad_compton, gam_gap)

# characteristic one-time scattered photon energy
#xcomp = x2+x2*(4.0/3.0)*gam_rad_compton**2
#xcomp = x2*(4.0/3.0)*gam_rad_compton**(2*3)
#xcomp = x2*gam_rad_compton**2/(1+gam_rad_compton*x2)

#gam_scaling = [1.0]
gam_scaling = [1e-4,1e-3,1e-2,1e-1,1e0]
linestyles = ["solid","dashed","dashdot","dotted",(0, (1, 1))]

lmfp_1phot_arr = []*len(gam_scaling)
lmfp_2phot_arr = []*len(gam_scaling)
lmfp_compton_arr = []*len(gam_scaling)

for ig in range(0,len(gam_scaling)):

    xcomp = x2*(gam_scaling[ig]*gam_smallest)**2/(1+gam_scaling[ig]*gam_smallest*x2)

    #xcomp = x2+x2*(4.0/3.0)*(0.1*gam_rad_compton)**2
    #xcomp = x2*(0.05*gam_rad_compton)**2/(1+(0.05*gam_rad_compton)*x2)
    #xcomp = np.minimum(xcomp,gam_rad_compton)
    print("xcomp [m_e x c^2]: ",xcomp)
    print()

    P_C_KN = 4.0*sigT*c*(gam_scaling[ig]*gam_smallest)**2*U_x*f_kn_exact(gam_scaling[ig]*gam_smallest,a)/3.0
    lmfp_compton_KN = gam_scaling[ig]*gam_smallest*m_e*c**3/P_C_KN
    lmfp_compton_arr.append(lmfp_compton_KN)


    if setup=="rp":
        xsyn = 1.5*bratio*(r_g/rad_curv)*(gam_scaling[ig]*gam_smallest)**3    
        x1 = xsyn
    else:
        x1 = xcomp

    #Next 1-photon mean free path:
    from estimation_functions import comp_lmfp_1phot
    lmfp_1phot, lmfp_1phot_approx = comp_lmfp_1phot(x1,bratio,rad_curv,period)
    print("lmfps_1phot/Rpc: ", lmfp_1phot/r_pcap)
    lmfp_1phot_arr.append(lmfp_1phot)

    #Next 2-photon mean free path (between curvature and thermal photons):
    from estimation_functions import comp_lmfp_2phot
    lmfp_2phot = comp_lmfp_2phot(x1,tkev)
    print("lmfps_2phot/Rpc: ", lmfp_2phot/r_pcap)
    lmfp_2phot_arr.append(lmfp_2phot)

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
    #print(gam_gap, 2.0*np.pi**2*echarge*Bfield*r_stlinestylesar**3/(m_e*c**4*period**2))
    #print(gam_rad_compton, Const1*Bfield**0.5*period**(-0.75)*r_star**0.75*(tkev*kev2erg)**-2)
    #print(gam_rad_syn, Const2*Bfield**0.25*period**(-3.0/8.0)*r_star**(3.0/8.0)*rad_curv**0.5)
    #print("Const1:",Const1)
    #print("(Const1/Const2)**0.5:",(Const1/Const2)**0.5)

    #print("T_gam_rad_equal [keV], T_gam_rad_equal [K]: ", T_gam_rad_equal, T_gam_rad_equal*keV2Kelvin)

    prefac1 = (3.0*m_e*c**5*hplanck**3)**0.5/(4.0*51.9*np.pi*sigT*Const1)**0.5*(2.0*np.pi/c)**(-0.25)
    T_comp_lmpf1 = prefac1*Bfield**(-0.25)*r_star**(-9.0/8.0)*period**(5.0/8.0)/kev2erg

    #print("Const3: ",prefac1)

    #print("T_comp_lmpf1 [keV], T_comp_lmpf1 [K]: ", T_comp_lmpf1, T_comp_lmpf1*keV2Kelvin)

if tgrid_plot:

    import matplotlib
    matplotlib.use('agg')
    from pylab import *

    rc("text", usetex=True)
    fig = figure(figsize=(10,10), dpi=300)
    plt.figure(1)
    lbfontsz = 30
    lwidth= 2.0
    rc("xtick", labelsize=lbfontsz)
    rc("ytick", labelsize=lbfontsz)
    rc("axes", linewidth=lwidth)
    ax = plt.subplot(1,1,1)

    for ig in range(0,len(gam_scaling)):
        ax.loglog(tkev*keV2Kelvin, lmfp_2phot_arr[ig]/r_pcap, color="g", linestyle = linestyles[ig], linewidth=lwidth, label="2-photon")
        ax.loglog(tkev*keV2Kelvin, lmfp_1phot_arr[ig]/r_pcap, color="b", linestyle = linestyles[ig], linewidth=lwidth, label="1-photon")
        #ax.loglog(tkev*keV2Kelvin, lmfp_compton_arr[ig]/r_pcap, color="r", linestyle = linestyles[ig], linewidth=lwidth, label="Compton")

        if setup=="wd":
            x_target1 = 0.07*keV2Kelvin
            y_target1 = np.interp(x_target1, tkev*keV2Kelvin, lmfp_2phot_arr[ig]/r_pcap)
            y_target2 = np.interp(x_target1, tkev*keV2Kelvin, lmfp_1phot_arr[ig]/r_pcap)    
            ax.loglog(x_target1, y_target1, 'og', ms=8)
            ax.loglog(x_target1, y_target2, 'ob', ms=8)
        else:
            x_target1 = 0.1*keV2Kelvin
            x_target2 = 0.5*keV2Kelvin        
            y_target1 = np.interp(x_target1, tkev*keV2Kelvin, lmfp_2phot_arr[ig]/r_pcap)
            y_target2 = np.interp(x_target1, tkev*keV2Kelvin, lmfp_1phot_arr[ig]/r_pcap)    
            ax.loglog(x_target1, y_target1, 'og', ms=8)
            ax.loglog(x_target1, y_target2, 'ob', ms=8)    
            y_target1 = np.interp(x_target2, tkev*keV2Kelvin, lmfp_2phot_arr[ig]/r_pcap)
            y_target2 = np.interp(x_target2, tkev*keV2Kelvin, lmfp_1phot_arr[ig]/r_pcap)    
            ax.loglog(x_target2, y_target1, 'og', ms=8)
            ax.loglog(x_target2, y_target2, 'ob', ms=8)  

    ax.axhline(y=1.0, color='grey', linestyle='-',linewidth=lwidth)

    ax.set_xlabel("$T_{\mathrm{BB}} [K]$",fontsize=lbfontsz)
    ax.set_ylabel("$l_{\mathrm{mfp}}/H_{\mathrm{gap}}$",fontsize=lbfontsz)
    #ax.set_ylim(1e-3,1e3)
    ax.set_ylim(1e-3,1e6)
    
    if setup=="wd":
          ax.set_xlim(5e4,1e6)  
    else:
        ax.set_xlim(8e4,1e7)


    # --- add top axis ---
    secax = ax.secondary_xaxis(
        'top',
        functions=(lambda K: K / keV2Kelvin,
                   lambda keV: keV * keV2Kelvin)
    )
    secax.set_xlabel("$k_\mathrm{B}T_{\mathrm{BB}} [\mathrm{keV}]$",fontsize=lbfontsz,labelpad=12) 

    if setup=="wd":
        secax.set_xticks([0.005, 0.01, 0.02, 0.05])
        secax.set_xticklabels(["$0.005$","$0.01$","$0.02$","$0.05$"])    
    else:
        secax.set_xticks([0.01, 0.03, 0.1, 0.2, 0.5])
        secax.set_xticklabels(["$0.01$","$0.03$","$0.1$","$0.2$", "$0.5$"])

    ax.tick_params(axis='both', which='both', width=lwidth, length=5.0) 
    secax.tick_params(axis='both', which='both', width=lwidth,length=5.0) 
    #ax.set_title("$\gamma_{\mathrm{rad}}$",fontsize=lbfontsz)

    #legend_lines = [
    #    Line2D([0], [0], color="g", label="2-photon"),
    #    Line2D([0], [0], color="b", label="1-photon"),
    #]
    legend_lines = [
        Line2D([0], [0], color='black', linestyle=linestyles[0], label="$0.0001\gamma_{\mathrm{max}}$"),
        Line2D([0], [0], color='black', linestyle=linestyles[1], label="$0.001\gamma_{\mathrm{max}}$"),
        Line2D([0], [0], color='black', linestyle=linestyles[2], label="$0.01\gamma_{\mathrm{max}}$"),
        Line2D([0], [0], color='black', linestyle=linestyles[3], label="$0.1\gamma_{\mathrm{max}}$"),
        Line2D([0], [0], color='black', linestyle=linestyles[4], label="$\gamma_{\mathrm{max}}$"),                
    ]     
    plt.legend(handles=legend_lines, fontsize=22)

    #plt.legend(fontsize=22)
    
    boxlabel = "MSP"
    if setup=="wd":
        boxlabel = "WD"
    elif setup=="rp":
        boxlabel = "RP"
    
    ax.text(
        0.5, 0.95,                 # position
        boxlabel,
        transform=ax.transAxes,     # use axes coordinates
        fontsize=22,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8)
    )    
    
    plt.tight_layout()
    fig.savefig("lmfp_vs_x.png",dpi=300.0)

