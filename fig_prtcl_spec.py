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
    conf = Configuration(args.conf_filename, do_print=True)#False)
    fdir = conf.outdir + "/"

    print("reading {}".format(fdir))

    toolset = QEDToolset(conf)

    #--------------------------------------------------
    fig = plt.figure(1, figsize=(3.25, 2.5)) # single figure

    #fig = plt.figure(1, figsize=(3.25, 8.0)) # single figure
    #fig = plt.figure(1, figsize=(7.0,  5.5)) # two-column figure


    plt.rc('xtick', top   = True)
    plt.rc('ytick', right = True)

    plt.rc('font',  family='serif',)
    plt.rc('text',  usetex=False)

    plt.rc('xtick', labelsize=7)
    plt.rc('ytick', labelsize=7)
    plt.rc('axes',  labelsize=8)
    plt.rc('legend',  handlelength=4.0)

    #--------------------------------------------------
    if True: # regular gridspec
        nrow_fig = 1
        ncol_fig = 1

        gs = plt.GridSpec(nrow_fig, ncol_fig)

        gs.update(wspace = 0.2)
        gs.update(hspace = 0.3)

        axs = np.empty( (nrow_fig,ncol_fig), dtype=object)

        for j in range(ncol_fig):
            for i in range(nrow_fig):
                axs[i,j] = plt.subplot(gs[i,j])
                axs[i,j].minorticks_on()
                #axs[i,j].set_xlim((-2, 6.0))
                axs[i,j].set_xlim((10**-2, 10**6.0))                
            #axs[0,0].tick_params(labeltop=True)



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
    #tt = laps/conf.t_norm

    # narrow region above the star
    #hmin = int( 2 )
    #hmax = int( 20 )

    ## limits for the height-energy histogram
    hs = np.linspace(0, Lh, toolset.Nhist)

    hmin2 = 1
    hmax2 = 10
    print("hmin2/hmax2", hmin2, hmax2, hs[hmin2], hs[hmax2])

    hmin3i = int( conf.rad_pcap*0.6 )
    hmax3i = int( conf.rad_pcap*1.0 )

    #hmin3i = int( conf.rad_pcap*0.6 )
    #hmax3i = int( conf.rad_pcap*0.8 )

    # limits for the height-energy histogram
    hs = np.linspace(0, Lh, toolset.Nhist)
    for i in range(len(hs)):
        if hs[i] > hmin3i: break
        hmin3 = i

    for i in range(len(hs)):
        if hs[i] > hmax3i: break
        hmax3 = i

    print("hmin3/hmax3", hmin3, hmax3, hs[hmin3], hs[hmax3])

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
    #print(xlow_mask2)

    #--------------------------------------------------
    def integrate(xs, ys):
        # trapedzoidal
        #dxs = np.diff(xs)
        dxs = xs[1] - xs[0] # assume uniform grid
        return np.sum(ys*dxs)


    tmin = 5.5
    tmax = 7.5
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
    
    Nlaps = 0
    for lap in laps:

        if not(tmin < lap/conf.t_norm < tmax):
            continue

        col = cmap(norm( lap/conf.t_norm - tmin ))

        #--------------------------------------------------
        # particle spectra

        # read from output file
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
        # he and hp store particles as histograms of dn/dlnp

        # histogram into units of n_GJ
        n_units_return = toolset.N_box/toolset.N_wgt # de-unitize what we have in the qed_toolset
        n_units_return *= 1/(hmax2-hmin2) # normalize to per hmax2-hmin2 histrogram bins
        n_units_return *= toolset.Nhist/conf.Lx # normalize to per cell
        n_units_return *= 1/conf.ppc # normalize to n_GJ

        n_units_thermal = n_units_return*(hmax2-hmin2)
        n_units_gap = n_units_thermal/(hmax3-hmin3)

        # in units of dN/dlnp 
        he_return = np.sum(hem[hmin2:hmax2,:], axis=0)*n_units_return #*zs2
        he_return = np.flip(he_return[:N]) # down

        #--------------------------------------------------
        hp_return = np.sum(hep[hmin2:hmax2,:], axis=0)*n_units_return #*zs2
        hp_return = np.flip(hp_return[:N]) # down

        he_thermal = np.sum(hem[0:1,:], axis=0)*n_units_thermal
        he_thermal = he_thermal[N:] # up

        he_gap = np.sum(hem[hmin3:hmax3,:], axis=0)*n_units_gap
        he_gap = he_gap[N:] + np.flip(he_gap[:N])

        hp_gap = np.sum(hep[hmin3:hmax3,:], axis=0)*n_units_gap
        hp_gap = hp_gap[N:] + np.flip(hp_gap[:N])

        # sum for total return spec
        he_int[:] += he_return
        hp_int[:] += hp_return

        he_thermal_int[:] += he_thermal

        Nlaps = Nlaps + 1

        # coarse grain spec with every merge step
        if i == 0:
            heM = he_gap
            hpM = hp_gap
        else:
            heM += he_gap
            hpM += hp_gap
        i += 1

        if i == merge:
            #axs[0,0].plot(10**lnzs, (heM + hpM)/merge, color=col, lw=0.8)
            axs[0,0].plot(10**lnzs, (heM + hpM)/merge/np.log(10), color=col, lw=0.8)            
            i = 0


    # plot total spec as well
    #axs[0,0].plot(10**lnzs, (he_int + hp_int)/Nlaps, color=col, lw=1.5)
    #axs[0,0].plot(10**lnzs, (he_thermal_int)/Nlaps, color=col, linestyle="dashdot", lw=1.5)

    axs[0,0].plot(10**lnzs, (he_int + hp_int)/Nlaps/np.log(10), color=col, lw=1.5)
    axs[0,0].plot(10**lnzs, (he_thermal_int)/Nlaps/np.log(10), color=col, linestyle="dashdot", lw=1.5)

    print("Energy integrated multiplicity: ",integrate(lnzs, (he_int + hp_int)/Nlaps))

    #print("integrated he:", he_int)
    #print("integrated hp:", hp_int)

    #-------------------------------------------------- 
    #axs[0,0].set_xlabel(r"$\log_{10} \gamma\beta$")
    axs[0,0].set_xlabel(r"$\gamma\beta$")    

    #axs[0,0].set_ylabel(r"$  \mathrm{d} m_\pm/ \mathrm{d} \, \log_{10}(\gamma\beta) $")
    axs[0,0].set_ylabel(r"$  \mathrm{d} m_\pm/ \mathrm{d} \, \ln(\gamma\beta) $")

    axs[0,0].set_yscale("log")
    axs[0,0].set_xscale("log")

    axs[0,0].set_ylim((1e-2, 1e3))

    #--------------------------------------------------
    # manual slope
    #gams = np.array([1e4, 8e6])
    #dmdg = gams**-0.3
    gams = np.array([1e2, 1e5])
    dmdg = gams**-0.9
    dmdg[:] *= 1e2/dmdg[0]
    #axs[0,0].plot(np.log10(gams), dmdg, color="k", linestyle="dashed", lw=1.5)
    #axs[0,0].plot(gams, dmdg, color="k", linestyle="dashed", lw=1.5)
    axs[0,0].plot(gams, dmdg/np.log(10), color="k", linestyle="dashed", lw=1.5)

    #--------------------------------------------------
    #manual peak
    #gams = 400.0*np.ones(2)
    #dmdg = [3e5, 5e6]
    gams = 30.0*np.ones(2)
    dmdg = [9.5e1, 9e2]
    #axs[0,0].plot(np.log10(gams), dmdg, color="k", linestyle="dashed", lw=1.5)
    #axs[0,0].plot(gams, dmdg, color="k", linestyle="dashed", lw=1.5)
    axs[0,0].plot(gams, dmdg/np.log(10), color="k", linestyle="dashed", lw=1.5)

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
        cb1.set_label(r'$(t-t_{0})/t_\mathrm{esc}$')

    pos = axs[0,0].get_position()
    print('ax pos:', pos)
    fig.subplots_adjust(left=axleft, bottom=axbottom, right=axright, top=axtop)

    fname = fdir + 'fig_prtcl_spec.pdf'
    plt.savefig(fname, dpi=300)



