import numpy as np
import sys, os

import matplotlib

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import get_cmap
from matplotlib import colorbar

from numpy import pi, log, log10, exp, sin, cos 
from scipy.signal import savgol_filter


# runko tools
import pytools
from init_problem import Configuration_Pulsar as Configuration
import h5py as h5

from qed_toolset import QEDToolset


import warnings
#warnings.filterwarnings("once") # suppress warnings after reporting them once
warnings.filterwarnings('ignore')

#-------------------------------------------------- 
if __name__ == "__main__":


    #--------------------------------------------------
    # command line driven version
    args = pytools.parse_args()
    conf = Configuration(args.conf_filename, do_print=False)#False)
    fdir = conf.outdir + "/"

    print("reading {}".format(fdir))

    toolset = QEDToolset(conf)



    #--------------------------------------------------
    fig = plt.figure(1, figsize=(2.0, 8.0)) # narrow figure

    #fig = plt.figure(1, figsize=(3.25, 8.0)) # single figure
    #fig = plt.figure(1, figsize=(7.0,  5.5)) # two-column figure

    plt.rc('xtick', top   = True)
    plt.rc('ytick', right = True)

    plt.rc('font',  family='serif',)
    plt.rc('text',  usetex=True)

    plt.rc('xtick', labelsize=6)
    plt.rc('ytick', labelsize=6)
    plt.rc('axes',  labelsize=8)
    plt.rc('legend',  handlelength=4.0)

    #--------------------------------------------------
    if True: # regular gridspec
        nrow_fig = 7
        ncol_fig = 1

        gs = plt.GridSpec(nrow_fig, ncol_fig)

        gs.update(wspace = 0.2)
        gs.update(hspace = 0.12)

        axs = np.empty( (nrow_fig,ncol_fig), dtype=object)

        for j in range(ncol_fig):
            for i in range(nrow_fig):
                axs[i,j] = plt.subplot(gs[i,j])
                axs[i,j].minorticks_on()

                #axs[i,j].set_yscale('log')

                if i <= nrow_fig-2:
                    axs[i,j].tick_params(labelbottom=False)

            #axs[0,0].tick_params(labeltop=True)


    axs[0,0].set_ylabel(r"$p_-$ ($m_e c$)")
    axs[1,0].set_ylabel(r"$p_+$ ($m_e c$)")
    axs[2,0].set_ylabel(r"$p_i$ ($m_e c$)")
    axs[3,0].set_ylabel(r"$x$ ($m_e c^2$)")

    axs[4,0].set_ylabel(r"$m_{\pm/x}$")
    #axs[4,0].set_ylabel(r"$\mathrm{d} m/\mathrm{d} x$")
    #axs[3,0].set_ylabel(r"$m(x)/ \Delta x$")

    axs[5,0].set_ylabel(r"$\langle \gamma \rangle$")
    axs[6,0].set_ylabel(r"$E/E_\mathrm{rot}$")

    #axs[4,0].set_ylabel(r"$j/j_m$")
    #axs[5,0].set_ylabel(r"$B/B_0 - 1$")
    #axs[5,0].set_ylabel(r"$\langle E \rangle_\mathrm{LF}/E_\mathrm{rot}$")


    hmin = -0.05  #0.0
    hmax = 1.25 #1.2 #0.05 #0.05 #7.0
    #hmin = -0.001
    #hmax = 0.005 #1.2 #0.05 #0.05 #7.0
    #hmin = 0.55
    #hmax = 0.65 
    #hmin = 1.15
    #hmax = 1.30

    #hmin = -0.01
    #hmax = +0.01
    #hmin = 1.25-0.01
    #hmax = 1.25+0.01

    for j in range(ncol_fig):
        for i in range(nrow_fig):
            axs[i,j].set_xlim((hmin, hmax))

    #for j in range(ncol_fig):
    #    for i in range(nrow_fig):
    #        #axs[i,j].set_xlim((0, 0.1)) # left start
    #        axs[i,j].set_xlim((0.9, 1.0)) #right end

    #axs[4,0].set_ylim((1e-1, 1e3))
    axs[4,0].set_ylim((1e-3, 1e8))
    axs[5,0].set_ylim((1.0, 1e8))

    axs[6,0].set_ylim((-0.1, 0.6))
    #axs[6,0].set_ylim((-1.1, 1.1))
    #axs[6,0].set_ylim((0.0, 0.2))
    #axs[6,0].set_ylim((-0.01, 0.01))


    #axs[5,0].set_ylim((-0.00005, 0.00005))

    #axs[5,0].set_ylim((-1, 1))
    #axs[5,0].set_yscale("symlog", linthresh=1e-7)

    axs[4,0].set_yscale("log")
    axs[5,0].set_yscale("log")
    axs[6,0].set_xlabel(r"$x/H_\mathrm{gap}$")



    # grid configuration
    if conf.oneD:
        cenx   = -conf.rad_star + conf.rad_curv_shift
        ceny   = 0 
        cenz   = 0 
        Lh = conf.Lx
    elif conf.twoD:
        cenx   = conf.Lx//2 + 0.5
        ceny   = -conf.rad_star + conf.rad_curv_shift
        cenz   = 0 
        Lh = conf.Ly
    elif conf.threeD:
        cenx   = conf.Lx//2 + 0.5
        ceny   = -conf.rad_star + conf.rad_curv_shift
        Lh = conf.Lz


    #--------------------------------------------------
    # particle spectra
    if True:

        # read from output file
        fname = conf.outdir + '/qed_{}.h5'.format(str(args.lap))
        f5 = h5.File(fname,'r')

        hem = f5['h2_ene_e-'][()]
        hep = f5['h2_ene_e+'][()]
        hip = f5['h2_ene_p'] [()]
        hph = f5['h2_ene_ph'][()]

        f5.close()

        #--------------------------------------------------

        hhlims = [1,1]
        hhlims[0] = 0 #-(conf.rad_curv_shift + conf.height_atms)/conf.rad_pcap #(0 - ceny)/conf.rad_pcap
        hhlims[1] = Lh/conf.rad_pcap

        hh = np.linspace(hhlims[0], hhlims[1], toolset.Nhist) # height grid for the histogram data

        print('hhlims', hhlims)
        print(toolset.pxlims[1],  toolset.pxlims[0] )
        px_log_extent = toolset.pxlims[1] - toolset.pxlims[0] 
        print('px_log_extent:', px_log_extent)

        for i in range(0,4): 
            axs[i,0].set_ylim((-px_log_extent, px_log_extent))

        print(toolset.xxlims[1],  toolset.xxlims[0] )
        xx_log_extent = toolset.xxlims[1] - toolset.xxlims[0] 
        print('xx_log_extent:', xx_log_extent)

        #electrons
        imem = axs[0,0].imshow( np.zeros((toolset.Nhist, 2*toolset.Nhist)),
                         extent=[hhlims[0], hhlims[1], 
                                 -px_log_extent, px_log_extent], 
                         origin='lower',
                         cmap='turbo',
                         aspect='auto',
                         interpolation='nearest',
                         vmin=np.log10(toolset.pylims[0]),
                         vmax=np.log10(toolset.pylims[1]),
                         zorder=1,
                         rasterized=True,
                         )

        #positrons
        imep = axs[1,0].imshow( np.zeros((toolset.Nhist, 2*toolset.Nhist)),
                         extent=[hhlims[0], hhlims[1], 
                                 -px_log_extent, px_log_extent], 
                         origin='lower',
                         cmap='turbo',
                         aspect='auto',
                         interpolation='nearest',
                         vmin=np.log10(toolset.pylims[0]),
                         vmax=np.log10(toolset.pylims[1]),
                         zorder=1,
                         rasterized=True,
                         )

        #protons
        imip = axs[2,0].imshow( np.zeros((toolset.Nhist, 2*toolset.Nhist)),
                         extent=[hhlims[0], hhlims[1], 
                                 -px_log_extent, px_log_extent], 
                         origin='lower',
                         cmap='turbo',
                         aspect='auto',
                         interpolation='nearest',
                         vmin=np.log10(toolset.pylims[0]),
                         vmax=np.log10(toolset.pylims[1]),
                         zorder=1,
                         rasterized=True,
                         )


        # height vs ene; photons
        imph = axs[3,0].imshow( np.zeros((toolset.Nhist, 2*toolset.Nhist)),
                         extent=[hhlims[0], hhlims[1], 
                                 -xx_log_extent, xx_log_extent], 
                          origin='lower',
                          cmap='turbo',
                          aspect='auto',
                          interpolation='nearest',
                          vmin=np.log10( toolset.xylims[0]),
                          vmax=np.log10( toolset.xylims[1]),
                          zorder=1,
                          rasterized=True,
                          )

        imem.set_data(np.log10( hem.T)) 
        imep.set_data(np.log10( hep.T)) 
        imip.set_data(np.log10( hip.T)) 
        imph.set_data(np.log10( hph.T)) 


        #--------------------------------------------------

        # momentum axis
        lnzs  = toolset.lnzs
        lnzs2 = np.array( [-np.flip(lnzs), lnzs ] ).flatten()
        zs    = toolset.zs
        gs    = np.sqrt( 1 + toolset.zs**2 )

        # photon energy axes; turn into [-val,+val]
        lnxs = toolset.lnxs
        xs = toolset.xs

        lnxs2 = np.array( [-np.flip(lnxs), lnxs ] ).flatten()
        xs2   = np.array( [-np.flip(xs), xs ] ).flatten()

        #debugging mask for low-energy photon cut
        xlow_mask2 = np.ones_like(xs2)

        if conf.qed_mode_rp:
            for i in range(len(xs2)):
                if np.abs(xs2[i]) < 2.0:
                    xlow_mask2[i] = 0

        #--------------------------------------------------
        def integrate(xs, ys): # trapedzoidal
            dxs = xs[1] - xs[0] # assume uniform grid
            return np.sum(ys*dxs)


        # histogram into units of n_GJ
        n_units = toolset.N_box/toolset.N_wgt # de-unitize what we have in the qed_toolset
        n_units *= 1.0/conf.ppc # normalize to n_GJ
        n_units *= toolset.Nhist/conf.Lx # normalize to per cell 

        # conversion factor into units of n_GJ
        nx_units = toolset.N_box*toolset.N_time/toolset.N_wgt # de-unitize what we have in qed_toolset
        nx_units *= 1/conf.ppc      # normalize to n_GJ
        nx_units *= conf.Lx/toolset.Nhist # normalize to per cell 

        #NOTE: after these units, integral over the x axis gives total number of particles in units of multiplicity


        #--------------------------------------------------
        # projections as a function of height
        N = toolset.Nhist

        he_p =          hem[:, N:]  # positive velocities
        he_m = np.flip( hem[:, :N], axis=1) # negative velocities

        hp_p =          hep[:, N:]  # positive velocities
        hp_m = np.flip( hep[:, :N], axis=1) # negative velocities

        hi_p =          hip[:, N:]  # positive velocities
        hi_m = np.flip( hip[:, :N], axis=1) # negative velocities

        # NOTE: we also mask x < 2 values out from the count
        #hx = np.sum(hph[hmin2:hmax2,:], axis=0)

        hph = xlow_mask2*hph/xs2**2 # mask low energy photons and normalize hist to dN/dx units
        hx_p =          hph[:, N:]          # positive vels
        hx_m = np.flip( hph[:, :N], axis=1) # neg vels


        # multiplicities for positive velocities (_p) and negative (_m) velocities
        me_p = np.zeros(N)
        me_m = np.zeros(N)
        mp_p = np.zeros(N)
        mp_m = np.zeros(N)
        mi_p = np.zeros(N)
        mi_m = np.zeros(N)
        mx_p = np.zeros(N)
        mx_m = np.zeros(N)

        # mean gamma factor <\gamma> for positive/negative vels

        # he/hm/hx are particle distribution functions represented as h(x, ln z) 
        gam_me_p = np.zeros(N)
        gam_me_m = np.zeros(N)
        gam_mp_p = np.zeros(N)
        gam_mp_m = np.zeros(N)
        gam_mi_p = np.zeros(N)
        gam_mi_m = np.zeros(N)
        gam_mx_p = np.zeros(N)
        gam_mx_m = np.zeros(N)


        # <\gamma^-3>
        g3_me_p = np.zeros(N)
        g3_me_m = np.zeros(N)
        g3_mp_p = np.zeros(N)
        g3_mp_m = np.zeros(N)
        g3_mi_p = np.zeros(N)
        g3_mi_m = np.zeros(N)


        # loop over spatial coordinate
        for i in range(N):
            me_p[i] = integrate(lnzs, he_p[i,:])*n_units
            me_m[i] = integrate(lnzs, he_m[i,:])*n_units

            mp_p[i] = integrate(lnzs, hp_p[i,:])*n_units
            mp_m[i] = integrate(lnzs, hp_m[i,:])*n_units

            mi_p[i] = integrate(lnzs, hi_p[i,:])*n_units
            mi_m[i] = integrate(lnzs, hi_m[i,:])*n_units

            mx_p[i] = integrate(lnxs, hx_p[i,:])*nx_units
            mx_m[i] = integrate(lnxs, hx_m[i,:])*nx_units

            gam_me_p[i] = integrate(lnzs, gs*he_p[i,:])/integrate(lnzs, he_p[i,:])
            gam_me_m[i] = integrate(lnzs, gs*he_m[i,:])/integrate(lnzs, he_m[i,:])

            gam_mp_p[i] = integrate(lnzs, gs*hp_p[i,:])/integrate(lnzs, hp_p[i,:])
            gam_mp_m[i] = integrate(lnzs, gs*hp_m[i,:])/integrate(lnzs, hp_m[i,:])

            gam_mi_p[i] = integrate(lnzs, gs*hi_p[i,:])/integrate(lnzs, hi_p[i,:])
            gam_mi_m[i] = integrate(lnzs, gs*hi_m[i,:])/integrate(lnzs, hi_m[i,:])

            gam_mx_p[i] = integrate(lnzs, xs*hx_p[i,:])/integrate(lnzs, hx_p[i,:])
            gam_mx_m[i] = integrate(lnzs, xs*hx_m[i,:])/integrate(lnzs, hx_m[i,:])
 
            #g3 = integrate(lnzs, (zs**(-3))*he_m[i,:])*n_units
            #g  = integrate(lnzs,         zs*he_m[i,:])*n_units
            #g3_me_p[i] = g3/g**-3
            #print(i, g3, g**-3, g3**(-1/3), g)

            #g3_me_p[i] = ( integrate(lnzs, (zs**(-3))*he_p[i,:])*n_units ) / ( integrate(lnzs, zs*he_p[i,:])*n_units )
            g3_me_p[i] = integrate(lnzs, (gs**(-3))*he_p[i,:]) / integrate(lnzs, he_p[i,:])
            g3_me_m[i] = integrate(lnzs, (gs**(-3))*he_m[i,:]) / integrate(lnzs, he_m[i,:])
            g3_mp_p[i] = integrate(lnzs, (gs**(-3))*hp_p[i,:]) / integrate(lnzs, hp_p[i,:])
            g3_mp_m[i] = integrate(lnzs, (gs**(-3))*hp_m[i,:]) / integrate(lnzs, hp_m[i,:])
            g3_mi_p[i] = integrate(lnzs, (gs**(-3))*hi_p[i,:]) / integrate(lnzs, hi_p[i,:])
            g3_mi_m[i] = integrate(lnzs, (gs**(-3))*hi_m[i,:]) / integrate(lnzs, hi_m[i,:])



        if True: # smooth values 
            def smooth(v): # smooth and mask nan's away
                m = np.argwhere(np.isnan(v))
                v1 = np.array(v, copy=True)
                v1[m] = 0.0
                v1 = savgol_filter(v1, int(0.05*len(v)), 1)
                v1[m] = np.nan
                return v1

            me_p = smooth(me_p)
            me_m = smooth(me_m)

            mp_m = smooth(mp_m)
            mp_p = smooth(mp_p)

            mi_m = smooth(mi_m)
            mi_p = smooth(mi_p)

            mx_m = smooth(mx_m)
            mx_p = smooth(mx_p)

            gam_me_m = smooth(gam_me_m)
            gam_me_p = smooth(gam_me_p)
            gam_mp_m = smooth(gam_mp_m)
            gam_mp_p = smooth(gam_mp_p)
            gam_mi_m = smooth(gam_mi_m)
            gam_mi_p = smooth(gam_mi_p)

        #print(g3_me_p)

        lw = 0.8
        axs[4,0].plot(hh, me_p, color="C0", alpha=0.8, lw=lw, linestyle="solid")
        axs[4,0].plot(hh, me_m, color="C0", alpha=0.8, lw=lw, linestyle="dotted")

        axs[4,0].plot(hh, mp_p, color="C1", alpha=0.8, lw=lw, linestyle="solid")
        axs[4,0].plot(hh, mp_m, color="C1", alpha=0.8, lw=lw, linestyle="dotted")

        axs[4,0].plot(hh, mx_p, color="C2", alpha=0.8, lw=lw, linestyle="solid")
        axs[4,0].plot(hh, mx_m, color="C2", alpha=0.8, lw=lw, linestyle="dotted")

        axs[4,0].plot(hh, mi_p, color="C3", alpha=0.8, lw=lw, linestyle="solid")
        axs[4,0].plot(hh, mi_m, color="C3", alpha=0.8, lw=lw, linestyle="dotted")



        axs[5,0].plot(hh, gam_me_p, color="C0", lw=lw, linestyle="solid")
        axs[5,0].plot(hh, gam_me_m, color="C0", lw=lw, linestyle="dashed")

        axs[5,0].plot(hh, gam_mp_p, color="C1", lw=lw, linestyle="solid")
        axs[5,0].plot(hh, gam_mp_m, color="C1", lw=lw, linestyle="dashed")

        axs[5,0].plot(hh, gam_mi_p, color="C3", lw=lw, linestyle="solid")
        axs[5,0].plot(hh, gam_mi_m, color="C3", lw=lw, linestyle="dashed")

        #print("gam+", gam_mi_p)
        #print("gam-", gam_mi_m)

        #axs[3,0].plot(hh, gam_mx_p, color="C2", lw=lw, linestyle="solid")
        #axs[3,0].plot(hh, gam_mx_m, color="C2", lw=lw, linestyle="dashed")


        #axs[4,0].plot(hh, g3_me_p/gam_me_p**(-3), color="C0", lw=lw, linestyle="solid")
        #axs[4,0].plot(hh, g3_me_m/gam_me_m**(-3), color="C0", lw=lw, linestyle="dashed")
        #axs[4,0].plot(hh, g3_mp_p/gam_mp_p**(-3), color="C1", lw=lw, linestyle="solid")
        #axs[4,0].plot(hh, g3_mp_m/gam_mp_m**(-3), color="C1", lw=lw, linestyle="dashed")



        #-----------------------------------------------------------
        # end of prtcl analysis




    # custom tick formatter to manipulate the mangled values back to real values
    class CustomScalarFormatter(matplotlib.ticker.ScalarFormatter):
        def __init__(self, useOffset=None, useMathText=None, useLocale=None, replace_values=([],[])):
            super().__init__(useOffset=None, useMathText=None, useLocale=None)
            self.replace_values = replace_values
    
        def __call__(self, x, pos=None):
            """
            Return the format for tick value *x* at position *pos*.
            """
            if len(self.locs) == 0:
                return ''
            #elif x in self.replace_values[0]:
            #    idx = self.replace_values[0].index(x)
            #    return str(self.replace_values[1][idx])
            else:
                xp = (x - self.offset) / (10. ** self.orderOfMagnitude)
                if abs(xp) < 1e-8:
                    xp = 0

                if xp > 0:
                    xps = xp - 2 # NOTE: this assumes pxlims left value is -2
                    s = int(xps)
                    return f"$10^{s}$"
                if xp == 0:
                    return "$10^{-2}$" # NOTE: same here
                else:
                    xps = -xp - 2
                    s = int(xps)
                    return f"-$10^{s}$"

    axs[0,0].set_ylim((-px_log_extent, px_log_extent))
    axs[1,0].set_ylim((-px_log_extent, px_log_extent))
    axs[2,0].set_ylim((-px_log_extent, px_log_extent))
    axs[3,0].set_ylim((-xx_log_extent, xx_log_extent))

    #yticks = [-10, -6, -2, 2, 6, 10] # narrower grid for previous pxlims (-2, 6)
    #yticks = [-10, -6, -2, 2, 6, 10] # expanded grid for pxlims (-2, -7)
    yticks = [-10, -6, -2, 2, 6, 10] # expanded grid for pxlims (-2, -12)

    axs[0,0].set_yticks(yticks)
    axs[1,0].set_yticks(yticks)
    axs[2,0].set_yticks(yticks)

    majorformatter = CustomScalarFormatter()
    axs[0,0].yaxis.set_major_formatter(majorformatter)
    axs[1,0].yaxis.set_major_formatter(majorformatter)
    axs[2,0].yaxis.set_major_formatter(majorformatter)


    yticks = [-8, -6, -4, -2, 2, 4, 6, 8] # expanded grid for pxlims (-2, -12)
    axs[3,0].set_yticks(yticks)
    axs[3,0].yaxis.set_major_formatter(majorformatter)


    axs[0,0].axhline(y=0, linestyle='dashed', color='k', lw=0.4)
    axs[1,0].axhline(y=0, linestyle='dashed', color='k', lw=0.4)
    axs[2,0].axhline(y=0, linestyle='dashed', color='k', lw=0.4)
    axs[3,0].axhline(y=0, linestyle='dashed', color='k', lw=0.4)


    #--------------------------------------------------
    # debug additions
    if True:

        # find indixes of elements nearest to x
        def find_arg_nearest(arr, x):
            d = np.abs(arr - x)
            i = np.argmin(d)
            return i

        # real momenta array for pairs
        pvals = np.linspace(-px_log_extent, px_log_extent, 2*toolset.Nhist)
        garr = np.zeros_like(pvals)

        n = len(garr)
        print('mid-1', pvals[n//2-1]) # neg values
        print('mid  ', pvals[n//2])   # pos values

        # create a real energy array from the mangled xarr
        garr[0:n//2] = -10**(-pvals[0:n//2 ] - 2)
        garr[n//2:]  = +10**(+pvals[n//2:] - 2  )

        #print('garr:', garr)
        print('mid-1', garr[n//2-1]) # neg values
        print('mid  ', garr[n//2])   # pos values

        # same for photons
        xvals = np.linspace(-xx_log_extent, xx_log_extent, 2*toolset.Nhist)
        xarr = np.zeros_like(xvals)
        xarr[0:n//2] = -10**(-xvals[0:n//2 ] - 2)
        xarr[n//2:]  = +10**(+xvals[n//2:] - 2  )


        # now plot varios thresholds
        for gref in [conf.gam_gap, conf.gam_rad]:
            i1 = find_arg_nearest(garr, -gref)
            i2 = find_arg_nearest(garr, +gref)

            print('i1:', i1, 'g', garr[i1], 'ggap', -gref, np.log10(gref))
            print('i2:', i2, 'g', garr[i2], 'ggap', +gref, np.log10(gref))

            axs[0,0].axhline(y=pvals[i1], linestyle='solid', color='k', lw=0.4)
            axs[0,0].axhline(y=pvals[i2], linestyle='solid', color='k', lw=0.4)

            axs[1,0].axhline(y=pvals[i1], linestyle='solid', color='k', lw=0.4)
            axs[1,0].axhline(y=pvals[i2], linestyle='solid', color='k', lw=0.4)

            axs[2,0].axhline(y=pvals[i1], linestyle='solid', color='k', lw=0.4)
            axs[2,0].axhline(y=pvals[i2], linestyle='solid', color='k', lw=0.4)
        for xref in [conf.xsyn]:
            i1 = find_arg_nearest(xarr, -xref)
            i2 = find_arg_nearest(xarr, +xref)

            print('i1:', i1, 'x', xarr[i1], 'xref', -xref, np.log10(xref))
            print('i2:', i2, 'x', xarr[i2], 'xref', +xref, np.log10(xref))

            axs[3,0].axhline(y=xvals[i1], linestyle='solid', color='k', lw=0.4)
            axs[3,0].axhline(y=xvals[i2], linestyle='solid', color='k', lw=0.4)


    axs[5,0].axhline(y=conf.gam_gap, linestyle='solid', color='k', lw=0.4)
    axs[5,0].axhline(y=conf.gam_rad, linestyle='solid', color='k', lw=0.4)


    #--------------------------------------------------
    # field values
    if True:

        # read from output file
        fname = conf.outdir + '/flds_{}.h5'.format(str(args.lap))
        f5 = h5.File(fname,'r')

        ex = pytools.read_h5_array(f5, 'ex')/conf.e_norm
        ey = pytools.read_h5_array(f5, 'ey')/conf.e_norm
        ez = pytools.read_h5_array(f5, 'ez')/conf.e_norm

        jx = pytools.read_h5_array(f5, 'jx')/conf.j_norm
        jy = pytools.read_h5_array(f5, 'jy')/conf.j_norm
        jz = pytools.read_h5_array(f5, 'jz')/conf.j_norm

        bx = pytools.read_h5_array(f5, 'bx')/conf.b_norm
        by = pytools.read_h5_array(f5, 'by')/conf.b_norm
        bz = pytools.read_h5_array(f5, 'bz')/conf.b_norm

        print('shape ex', np.shape(ex))

        # reduce dimensions
        if conf.oneD:
            ex = np.mean(ex, axis=(1,2))
            ey = np.mean(ey, axis=(1,2))
            ez = np.mean(ez, axis=(1,2))

            jx = np.mean(jx, axis=(1,2))
            jy = np.mean(jy, axis=(1,2))
            jz = np.mean(jz, axis=(1,2))

            bx = np.mean(bx, axis=(1,2))
            by = np.mean(by, axis=(1,2))
            bz = np.mean(bz, axis=(1,2))
        else:
            print('TODO')

        e = ex + ey + ez
        b = bx + by + bz
        j = jx + jy + jz
        hh = np.linspace(hhlims[0], hhlims[1], Lh)

        print('e', e)
        print('j', j)
        print('b', b)

        # smooth current 
        if True:
            #jx = savgol_filter(jx, 300, 2)
            #jx = savgol_filter(jx, 300, 2)
            ex_flt = savgol_filter(ex, 2000, 1)
            ex[15000:18000] = ex_flt[15000:18000]
            #ex[16000:17000] = 1.0
        
        axs[6,0].plot(hh, ex, lw=0.8, linestyle='solid',  color='C0', alpha=0.8)
        #axs[6,0].plot(hh, ey, lw=0.8, linestyle='solid',  color='C1', alpha=0.8)
        #axs[6,0].plot(hh, ez, lw=0.8, linestyle='solid',  color='C2', alpha=0.8)
        #axs[6,0].plot(hh, bx, lw=0.8, linestyle='dotted', color='C3', alpha=0.8)

        #axs[6,0].plot(hh, ey, lw=0.8, linestyle='solid', color='C1', alpha=0.8)
        #axs[6,0].plot(hh, ez, lw=0.8, linestyle='solid', color='C2', alpha=0.8)


        #v = jx - jx[0]
        #print("dj:", v)
        #axs[6,0].plot(hh, v, lw=0.3, linestyle='solid', color='C1', alpha=0.8)

        #axs[6,0].plot(hh, jx, lw=0.3, linestyle='solid', color='C1', alpha=0.8)
        #axs[6,0].plot(hh, jy, lw=0.3, linestyle='solid', color='C1', alpha=0.8)
        #axs[6,0].plot(hh, jz, lw=0.3, linestyle='solid', color='C2', alpha=0.8)

        #ex_flt = savgol_filter(ex, 300, 1)
        #axs[6,0].plot(hh, ex_flt, lw=0.8, linestyle='solid', color='C0', alpha=0.8)
        #axs[6,0].plot(hh, ey, lw=0.8, linestyle='solid', color='C1', alpha=0.8)
        #axs[6,0].plot(hh, ez, lw=0.8, linestyle='solid', color='C2', alpha=0.8)
        #axs[6,0].plot(hh, bx-1.0, lw=0.8, linestyle='dashed', color='C0', alpha=0.8)

        print("jx L:", jx[0:15])
        print("jx R:", jx[-15:-1])

    #--------------------------------------------------
    for j in range(ncol_fig):
        for i in range(nrow_fig):
            axs[i,j].fill_between([-0.5, 0.0], -100, 1e20, color='gray', alpha=0.4, edgecolor=None) 

    #--------------------------------------------------
    # time stamp

    #print('r_pc:', conf.rad_pcap)
    print('t:', args.lap/conf.t_norm)
    stitle = r"$t/t_\mathrm{esc}$ = " + "{:3.1f}".format(args.lap/conf.t_norm)
    axs[0,0].set_title(stitle, fontsize=10)


    #--------------------------------------------------
    axleft    = 0.25
    axbottom  = 0.05
    axright   = 0.97
    axtop     = 0.95

    if False:
        pos1 = axs[0,0].get_position()

        axwidth  = axright - axleft
        axheight = (axtop - axbottom)*0.03
        axpad = 0.02
        cax = fig.add_axes([axleft, axtop + axpad, axwidth, axheight])

        cb1 = matplotlib.colorbar.ColorbarBase(
                cax,
                cmap=cmap,
                norm=norm,
                orientation='horizontal',
                ticklocation='top')
        cb1.set_label(r'power-law index $a$')

    pos = axs[0,0].get_position()
    print('ax pos:', pos)
    fig.subplots_adjust(left=axleft, bottom=axbottom, right=axright, top=axtop)

    print('-----lap:', args.lap, " tc/H:", args.lap/conf.t_norm)


    slap = str(args.lap).rjust(8, '0')

    fname = fdir + 'fig_gap2_' + slap + '_v2.pdf' 
    plt.savefig(fname)

    #fname = fdir + 'fig_zoom_gap2_' + slap + '.png' 
    #fname = fdir + 'fig_gap2_' + slap + '.png' 
    plt.savefig(fname, dpi=300)



