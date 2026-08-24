import numpy as np
import sys, os

import matplotlib

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import get_cmap
from matplotlib import colorbar

from numpy import pi, log, log10, exp, sin, cos 

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
    conf1 = Configuration(args.conf_filename, do_print=True)#False)
    fdir = conf1.outdir + "/"

    print("reading {}".format(fdir))

    #conf2 = Configuration("msp-jmfix_jm15_LT_1pp_rcsmall_iphot.ini", do_print=False)
    conf2 = Configuration("msp_1ph.ini", do_print=False)
    fdir2 = conf2.outdir + "/"
    print("reading {}".format(fdir2))


    confs = [conf1,conf2]

    #--------------------------------------------------
    # Increased vertical size to accommodate 2 subplots comfortably
    fig = plt.figure(1, figsize=(3.25, 4.5)) 

    plt.rc('xtick', top    = True)
    plt.rc('ytick', right = True)

    plt.rc('font',  family='serif',)
    plt.rc('text',  usetex=True)

    plt.rc('xtick', labelsize=7)
    plt.rc('ytick', labelsize=7)
    plt.rc('axes',  labelsize=8)
    plt.rc('legend',  handlelength=4.0)

    #--------------------------------------------------
    if True: # regular gridspec
        nrow_fig = 2 # Changed to 2 vertical subplots
        ncol_fig = 1

        gs = plt.GridSpec(nrow_fig, ncol_fig)

        gs.update(wspace = 0.2)
        gs.update(hspace = 0.0) # Reduced hspace to pack them closer

        axs = np.empty( (nrow_fig,ncol_fig), dtype=object)

        for j in range(ncol_fig):
            for i in range(nrow_fig):
                axs[i,j] = plt.subplot(gs[i,j])
                axs[i,j].minorticks_on()

                #axs[i,j].set_yscale('log')
                
                # Hide x-tick labels for the top subplot to prevent overlap
                if i <= nrow_fig-2:
                    axs[i,j].tick_params(labelbottom=False)

                axs[i,j].set_xlim((-2, 6.0))
            #axs[0,0].tick_params(labeltop=True)


    for iconf, conf in enumerate(confs):

        toolset = QEDToolset(conf)

        #--------------------------------------------------
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

        hhlims = [1,1]
        hhlims[0] = -(conf.rad_curv_shift + conf.height_atms)/conf.rad_pcap #(0 - ceny)/conf.rad_pcap
        hhlims[1] = Lh/conf.rad_pcap


        #--------------------------------------------------
        # time loop

        laps = np.array( list( range(0, conf.Nt, conf.interval) ) )

        hmin2 = 1
        hmax2 = 10
        
        # limits for the height-energy histogram
        hs = np.linspace(0, Lh, toolset.Nhist)
        print("hmin2/hmax2", hmin2, hmax2, hs[hmin2], hs[hmax2])

        # momentum axis
        lnzs  = toolset.lnzs
        zs    = toolset.zs
        lnzs2 = np.array( [-np.flip(lnzs), lnzs ] ).flatten()
        zs2   = np.array( [-np.flip(zs),     zs ] ).flatten()

        # photon energy axes; turn into [-val,+val]
        lnxs = toolset.lnxs
        xs = toolset.xs

        lnxs2 = np.array( [-np.flip(lnxs), lnxs ] ).flatten()
        xs2   = np.array( [-np.flip(xs), xs ] ).flatten()

        #debugging mask for low-energy photon cut
        xlow_mask2 = np.ones_like(xs2)
        for i in range(len(xs2)):
            if np.abs(xs2[i]) < 2.0:
                xlow_mask2[i] = 0

        #--------------------------------------------------
        def integrate(xs, ys):
            dxs = xs[1] - xs[0] # assume uniform grid
            return np.sum(ys*dxs)

        tmin = 0.0 #6.0 #0.0
        tmax = 4.0 #10.0 #2.0
        tspan = tmax - tmin

        norm = matplotlib.colors.Normalize(vmin=0, vmax=tspan)
        cmap = matplotlib.colormaps['turbo_r']

        merge = 12 # use this to tune how many curves are merged together
        i = 0

        # total time-integrated spec
        N = toolset.Nhist
        he_int = np.zeros(N)
        hp_int = np.zeros(N)
        hx_int = np.zeros(N)
        he_thermal_int = np.zeros(N)
        
        # Initialize array for total fake spec
        fake_int = np.zeros_like(lnzs) 
        
        Nlaps = 0
        for lap in laps:

            if not(tmin < lap/conf.t_norm < tmax):
                continue

            col = cmap(norm( lap/conf.t_norm - tmin ))

            #--------------------------------------------------
            # particle spectra
            fname = conf.outdir + '/qed_{}.h5'.format(str(lap))

            if not(os.path.isfile(fname)): continue

            print("lap", lap, "t", lap/conf.t_norm)

            f5 = h5.File(fname,'r')

            hem = f5['h2_ene_e-'][()]
            hep = f5['h2_ene_e+'][()]
            hph = f5['h2_ene_ph'][()]

            f5.close()

            #--------------------------------------------------
            #particles
            n_units = toolset.N_box/toolset.N_wgt 
            n_units *= 1/(hmax2-hmin2) 
            n_units *= toolset.Nhist/conf.Lx 
            n_units *= 1/conf.ppc 

            n_units_thermal = n_units*(hmax2-hmin2)

            # in units of dN/dlnp 
            he = np.sum(hem[hmin2:hmax2,:], axis=0)*n_units
            he = np.flip(he[:N]) # down

            hp = np.sum(hep[hmin2:hmax2,:], axis=0)*n_units
            hp = np.flip(hp[:N]) # down

            he_thermal = np.sum(hem[0:1,:], axis=0)*n_units_thermal
            he_thermal = he_thermal[N:] # up

            #--------------------------------------------------
            # photons
            nx_units = toolset.N_box*toolset.N_time/toolset.N_wgt
            nx_units *= 1/(hmax2-hmin2) 
            nx_units *= toolset.Nhist/conf.Lx 
            nx_units *= 1/conf.ppc 

            hx = np.sum(hph[hmin2:hmax2,:], axis=0)*nx_units
            hx = np.flip(hx[:N]) #"down

            # sum for total spec
            he_int[:] += he
            hp_int[:] += hp
            hx_int[:] += hx

            he_thermal_int[:] += he_thermal

            Nlaps = Nlaps + 1

            # coarse grain spec with every merge step
            if i == 0:
                heM = he
                hpM = hp
                hxM = hx
            else:
                heM += he
                hpM += hp
                hxM += hx
            i += 1

            if i == merge:
                axs[iconf,0].plot(lnzs, (heM + hpM)/merge, color=col, lw=0.8)
                #axs[1,0].plot(lnzs, (heM + hpM)/merge, color=col, lw=0.8)
                i = 0

        # plot total spec as well
        axs[iconf,0].plot(lnzs, (he_int + hp_int)/Nlaps, color=col, lw=1.5)
        axs[iconf,0].plot(lnzs, (he_thermal_int)/Nlaps, color=col, linestyle="dashdot", lw=1.5)
        

        #axs[1,0].plot(lnzs, (he_int + hp_int)/Nlaps, color=col, lw=1.5)
        #axs[1,0].plot(lnzs, (he_thermal_int)/Nlaps, color=col, linestyle="dashdot", lw=1.5)


        print("Energy integrated multiplicity: ",integrate(lnzs, (he_int + hp_int)/Nlaps))
        print("integrated he:", he_int)
        print("integrated hp:", hp_int)

        gamma = 10.0**lnzs
        mpl_avg = (he_int + hp_int)/Nlaps
        gamma_avg = integrate(lnzs,gamma * mpl_avg) / integrate(lnzs,mpl_avg)

        print("Average Lorentz factor: ",gamma_avg)


    # manual slope
    gams = np.array([1e2, 1e5])
    dmdg = gams**-0.8
    dmdg[:] *= 1e2/dmdg[0]
    axs[0,0].plot(np.log10(gams), dmdg, color="k", linestyle="dashed", lw=1.5)

    #--------------------------------------------------
    #manual peak
    gams = 50.0*np.ones(2)
    dmdg = [5e1, 5e2]
    axs[0,0].plot(np.log10(gams), dmdg, color="k", linestyle="dashed", lw=1.5)


    #-------------------------------------------------- 
    # Moved the X-axis label from axs[0,0] to axs[1,0]
    axs[1,0].set_xlabel(r"$\log_{10} \gamma\beta$")

    axs[0,0].set_ylabel(r"$  \mathrm{d} m_\pm/ \mathrm{d} \, \log_{10}(\gamma\beta) $")
    axs[1,0].set_ylabel(r"$  \mathrm{d} m_\pm/ \mathrm{d} \, \log_{10}(\gamma\beta) $")

    axs[0,0].set_yscale("log")
    axs[1,0].set_yscale("log")

    axs[0,0].set_ylim((1e-2, 1e3))
    axs[1,0].set_ylim((1e-2, 5e1))

    #--------------------------------------------------


    #--------------------------------------------------
    axleft    = 0.17
    axbottom  = 0.15
    axright   = 0.97
    axtop     = 0.80

    if True:
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
        cb1.set_label(r'$t/t_\mathrm{esc}$')

    pos = axs[0,0].get_position()
    print('ax pos:', pos)
    fig.subplots_adjust(left=axleft, bottom=axbottom, right=axright, top=axtop)

    fname = fdir + 'fig_return_prtcl_spec.pdf'
    plt.savefig(fname, dpi=300)
