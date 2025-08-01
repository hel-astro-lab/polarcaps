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
    fig = plt.figure(1, figsize=(2.2, 6.0)) # single figure

    #fig = plt.figure(1, figsize=(3.25, 8.0)) # single figure
    #fig = plt.figure(1, figsize=(7.0,  5.5)) # two-column figure

    plt.rc('xtick', top   = True)
    plt.rc('ytick', right = True)

    plt.rc('font',  family='serif',)
    plt.rc('text',  usetex=False)

    plt.rc('xtick', labelsize=6)
    plt.rc('ytick', labelsize=6)
    plt.rc('axes',  labelsize=8)
    plt.rc('legend',  handlelength=4.0)

    #--------------------------------------------------
    if True: # regular gridspec
        nrow_fig = 5
        ncol_fig = 1

        gs = plt.GridSpec(nrow_fig, ncol_fig)

        gs.update(wspace = 0.2)
        gs.update(hspace = 0.15)

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
    #axs[1,0].set_ylabel(r"$p_+$ ($m_e c$)")
    axs[1,0].set_ylabel(r"$x$ ($m_e c^2$)")

    axs[2,0].set_ylabel(r"$\mathrm{d} m/\mathrm{d} h$")
    axs[3,0].set_ylabel(r"$\langle p \rangle$, $\langle x \rangle$")

    axs[4,0].set_ylabel(r"$\langle \gamma^{-3} \rangle / \langle \gamma \rangle^{-3}$")
    #axs[4,0].set_ylabel(r"$\langle \gamma^{-3} \rangle$")

    hmin = -0.1
    hmax = 1.0 #7.0
    for j in range(ncol_fig):
        for i in range(nrow_fig):
            axs[i,j].set_xlim((hmin, hmax))

    axs[2,0].set_yscale("log")
    axs[3,0].set_yscale("log")
    axs[4,0].set_yscale("log")

    if conf.qed_mode_msp:
        axs[2,0].set_ylim((1e-1, 1e5))
    else:
        axs[2,0].set_ylim((1e-1, 1e3))

    axs[3,0].set_ylim((1e0, 1e8))

    axs[4,0].set_ylim((1e0, 1e12))


    #axs[4,0].set_ylim((-1.2, 3.2))
    #axs[5,0].set_ylim((-0.00005, 0.00005))

    axs[4,0].set_xlabel(r"$h/h_\mathrm{pc}$")

    for j in range(ncol_fig):
        for i in range(nrow_fig):
            axs[i,j].fill_between([-0.5, 0.0], -100, 1e20, color='gray', alpha=0.4, edgecolor=None) 


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
        hph = f5['h2_ene_ph'][()]

        f5.close()

        #--------------------------------------------------

        hhlims = [1,1]
        hhlims[0] = -(conf.rad_curv_shift + conf.height_atms)/conf.rad_pcap #(0 - ceny)/conf.rad_pcap
        hhlims[1] = Lh/conf.rad_pcap

        hh = np.linspace(hhlims[0], hhlims[1], toolset.Nhist) # height grid for the histogram data


        print('hhlims', hhlims)
        print(toolset.pxlims[1],  toolset.pxlims[0] )
        px_log_extent = toolset.pxlims[1] - toolset.pxlims[0] 
        print('px_log_extent:', px_log_extent)

        for i in range(0,2): 
            axs[i,0].set_ylim((-px_log_extent, px_log_extent))

        print(toolset.xxlims[1],  toolset.xxlims[0] )
        xx_log_extent = toolset.xxlims[1] - toolset.xxlims[0] 
        print('xx_log_extent:', xx_log_extent)


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
        def integrate(xs, ys):
            # trapedzoidal
            #dxs = np.diff(xs)
            dxs = xs[1] - xs[0] # assume uniform grid
            return np.sum(ys*dxs)



        # histogram into units of n_GJ
        n_units = toolset.N_box/toolset.N_wgt # de-unitize what we have in the qed_toolset
        n_units *= 1.0/conf.ppc # normalize to n_GJ
        n_units *= 1.0/toolset.Nhist # normalize to per histogram cell 

        # conversion factor into units of n_GJ
        nx_units = toolset.N_box*toolset.N_time/toolset.N_wgt # de-unitize what we have in qed_toolset
        nx_units *= 1/conf.ppc      # normalize to n_GJ
        nx_units *= 1/toolset.Nhist # normalize by area into xx per cell

        #NOTE: after these units, integral over the x axis gives total number of particles in units of multiplicity


        #--------------------------------------------------
        # projections as a function of height
        N = toolset.Nhist

        he_p =          hem[:, N:]  # positive velocities
        he_m = np.flip( hem[:, :N], axis=1) # negative velocities

        hp_p =          hep[:, N:]  # positive velocities
        hp_m = np.flip( hep[:, :N], axis=1) # negative velocities


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
        mx_p = np.zeros(N)
        mx_m = np.zeros(N)

        # mean gamma factor <\gamma> for positive/negative vels

        # he/hm/hx are particle distribution functions represented as h(x, ln z) 
        gam_me_p = np.zeros(N)
        gam_me_m = np.zeros(N)
        gam_mp_p = np.zeros(N)
        gam_mp_m = np.zeros(N)
        gam_mx_p = np.zeros(N)
        gam_mx_m = np.zeros(N)


        # <\gamma^-3>/<\gamma>^-3
        g3_me_p = np.zeros(N)
        g3_me_m = np.zeros(N)
        g3_mp_p = np.zeros(N)
        g3_mp_m = np.zeros(N)


        # loop over spatial coordinate
        for i in range(N):
            me_p[i] = integrate(lnzs, he_p[i,:])*n_units
            me_m[i] = integrate(lnzs, he_m[i,:])*n_units

            mp_p[i] = integrate(lnzs, hp_p[i,:])*n_units
            mp_m[i] = integrate(lnzs, hp_m[i,:])*n_units

            mx_p[i] = integrate(lnxs, hx_p[i,:])*nx_units
            mx_m[i] = integrate(lnxs, hx_m[i,:])*nx_units

            gam_me_p[i] = integrate(lnzs, zs*he_p[i,:])/integrate(lnzs, he_p[i,:])
            gam_me_m[i] = integrate(lnzs, zs*he_m[i,:])/integrate(lnzs, he_m[i,:])

            gam_mp_p[i] = integrate(lnzs, zs*hp_p[i,:])/integrate(lnzs, hp_p[i,:])
            gam_mp_m[i] = integrate(lnzs, zs*hp_m[i,:])/integrate(lnzs, hp_m[i,:])

            gam_mx_p[i] = integrate(lnzs, xs*hx_p[i,:])/integrate(lnzs, hx_p[i,:])
            gam_mx_m[i] = integrate(lnzs, xs*hx_m[i,:])/integrate(lnzs, hx_m[i,:])
 
            #g3 = integrate(lnzs, (zs**(-3))*he_m[i,:])*n_units
            #g  = integrate(lnzs,         zs*he_m[i,:])*n_units
            #g3_me_p[i] = g3/g**-3
            #print(i, g3, g**-3, g3**(-1/3), g)

            #g3_me_p[i] = ( integrate(lnzs, (zs**(-3))*he_p[i,:])*n_units ) / ( integrate(lnzs, zs*he_p[i,:])*n_units )

            g3_me_p[i] = integrate(lnzs, (gs**(-3))*he_p[i,:])*n_units / integrate(lnzs, he_p[i,:])
            g3_me_m[i] = integrate(lnzs, (gs**(-3))*he_m[i,:])*n_units / integrate(lnzs, he_m[i,:])
            g3_mp_p[i] = integrate(lnzs, (gs**(-3))*hp_p[i,:])*n_units / integrate(lnzs, hp_p[i,:])
            g3_mp_m[i] = integrate(lnzs, (gs**(-3))*hp_m[i,:])*n_units / integrate(lnzs, hp_m[i,:])



        #print(g3_me_p)

        lw = 0.5
        axs[2,0].plot(hh, me_p, color="C0", lw=lw, linestyle="solid")
        axs[2,0].plot(hh, me_m, color="C0", lw=lw, linestyle="dashed")

        axs[2,0].plot(hh, mp_p, color="C1", lw=lw, linestyle="solid")
        axs[2,0].plot(hh, mp_m, color="C1", lw=lw, linestyle="dashed")

        axs[2,0].plot(hh, mx_p, color="C2", lw=lw, linestyle="solid")
        axs[2,0].plot(hh, mx_m, color="C2", lw=lw, linestyle="dashed")


        axs[3,0].plot(hh, gam_me_p, color="C0", lw=lw, linestyle="solid")
        axs[3,0].plot(hh, gam_me_m, color="C0", lw=lw, linestyle="dashed")

        axs[3,0].plot(hh, gam_mp_p, color="C1", lw=lw, linestyle="solid")
        axs[3,0].plot(hh, gam_mp_m, color="C1", lw=lw, linestyle="dashed")

        axs[3,0].plot(hh, gam_mx_p, color="C2", lw=lw, linestyle="solid")
        axs[3,0].plot(hh, gam_mx_m, color="C2", lw=lw, linestyle="dashed")


        axs[4,0].plot(hh, g3_me_p/gam_me_p**(-3), color="C0", lw=lw, linestyle="solid")
        axs[4,0].plot(hh, g3_me_m/gam_me_m**(-3), color="C0", lw=lw, linestyle="dashed")
        axs[4,0].plot(hh, g3_mp_p/gam_mp_p**(-3), color="C1", lw=lw, linestyle="solid")
        axs[4,0].plot(hh, g3_mp_m/gam_mp_m**(-3), color="C1", lw=lw, linestyle="dashed")


        #--------------------------------------------------
        # histograms

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
                         )

        #imep = axs[1,0].imshow( np.zeros((toolset.Nhist, 2*toolset.Nhist)),
        #                 extent=[hhlims[0], hhlims[1], 
        #                         -px_log_extent, px_log_extent], 
        #                 origin='lower',
        #                 cmap='turbo',
        #                 aspect='auto',
        #                 interpolation='nearest',
        #                 vmin=np.log10(toolset.pylims[0]),
        #                 vmax=np.log10(toolset.pylims[1]),
        #                 zorder=1,
        #                 )

        # height vs ene; photons
        imph = axs[1,0].imshow( np.zeros((toolset.Nhist, 2*toolset.Nhist)),
                         extent=[hhlims[0], hhlims[1], 
                                 -xx_log_extent, xx_log_extent], 
                          origin='lower',
                          cmap='turbo',
                          aspect='auto',
                          interpolation='nearest',
                          vmin=np.log10( toolset.xylims[0]),
                          vmax=np.log10( toolset.xylims[1]),
                          zorder=1,
                          )

        imem.set_data(np.log10( hem.T)) 
        #imep.set_data(np.log10( hep.T)) 
        imph.set_data(np.log10( hph.T)) 


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

    #yticks = [-10, -6, -2, 2, 6, 10] # narrower grid for previous pxlims (-2, 6)
    yticks = [-10, -6, -2, 2, 6, 10] # expanded grid for pxlims (-2, -7)
    axs[0,0].set_yticks(yticks)
    axs[1,0].set_yticks(yticks)

    majorformatter = CustomScalarFormatter()
    axs[0,0].yaxis.set_major_formatter(majorformatter)
    axs[1,0].yaxis.set_major_formatter(majorformatter)


    #--------------------------------------------------
    # debug additions
    if True:

        # find indixes of elements nearest to x
        def find_arg_nearest(arr, x):
            d = np.abs(arr - x)
            i = np.argmin(d)
            return i

        # real momenta array for pairs
        xvals = np.linspace(-px_log_extent, px_log_extent, 2*toolset.Nhist)
        garr = np.zeros_like(xvals)

        n = len(garr)
        print('mid-1', xvals[n//2-1]) # neg values
        print('mid  ', xvals[n//2])   # pos values

        # create a real energy array from the mangled xarr
        garr[0:n//2] = -10**(-xvals[0:n//2 ] - 2)
        garr[n//2:]  = +10**(+xvals[n//2:] - 2  )

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

            print('i1:', i1, 'g', garr[i1], 'ggap', -gref)
            print('i2:', i2, 'g', garr[i2], 'ggap', +gref)

            axs[0,0].axhline(y=xvals[i1], linestyle='solid', color='C0', lw=0.4)
            axs[0,0].axhline(y=xvals[i2], linestyle='solid', color='C0', lw=0.4)

        for xref in [conf.xsyn]:
            i1 = find_arg_nearest(xarr, -gref)
            i2 = find_arg_nearest(xarr, +gref)

            print('i1:', i1, 'x', xarr[i1], 'xref', -xref)
            print('i2:', i2, 'x', xarr[i2], 'xref', +xref)

            axs[1,0].axhline(y=xvals[i1], linestyle='solid', color='C0', lw=0.4)
            axs[1,0].axhline(y=xvals[i2], linestyle='solid', color='C0', lw=0.4)




    #--------------------------------------------------
    # time stamp

    #print('r_pc:', conf.rad_pcap)
    print('t:', args.lap/conf.t_norm)
    stitle = r"$t c/H_\mathrm{pc}$ = " + "{:3.1f}".format(args.lap/conf.t_norm)
    axs[0,0].set_title(stitle, fontsize=10)


    #--------------------------------------------------
    axleft    = 0.25
    axbottom  = 0.09
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


    slap = str(args.lap).rjust(8, '0')

    #fname = fdir + 'fig_casc_' + slap + '.pdf' 
    #plt.savefig(fname)

    fname = fdir + 'fig_prtcl_profiles_' + slap + '.png' 
    plt.savefig(fname, dpi=300)



