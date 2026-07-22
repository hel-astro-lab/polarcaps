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
    conf = Configuration(args.conf_filename, do_print=False)
    fdir = conf.outdir + "/"

    print("reading {}".format(fdir))

    toolset = QEDToolset(conf)



    #--------------------------------------------------
    #fig = plt.figure(1, figsize=(3.25, 3.5)) # single figure
    fig = plt.figure(1, figsize=(3.25, 2.0)) # single figure

    #fig = plt.figure(1, figsize=(3.25, 2.5)) # single figure
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
        nrow_fig = 1 #2
        ncol_fig = 1

        gs = plt.GridSpec(nrow_fig, ncol_fig)

        gs.update(wspace = 0.2)
        gs.update(hspace = 0.33)

        axs = np.empty( (nrow_fig,ncol_fig), dtype=object)

        for j in range(ncol_fig):
            for i in range(nrow_fig):
                axs[i,j] = plt.subplot(gs[i,j])
                axs[i,j].minorticks_on()

                #axs[i,j].set_yscale('log')

                #if i <= nrow_fig-2:
                #    axs[i,j].tick_params(labelbottom=False)

            #axs[0,0].tick_params(labeltop=True)

    #axs[0,0].set_ylabel(r"$E/E_\mathrm{rot}$")
    #axs[1,0].set_ylabel(r"$k \hat{S}_{k} (k)$")
    axs[0,0].set_ylabel(r"$|\hat{E}_{k}^2|$")

    #axs[0,0].set_xlabel(r"$h/H_\mathrm{pc}$")
    axs[0,0].set_xlabel(r"$k H_\mathrm{gap}$")

    hmin = 0.0
    hmax = 1.2 #7.0
    #axs[0,0].set_xlim((hmin, hmax))

    #axs[i,j].set_xlim((hmin, hmax))

    #axs[0,0].set_ylim((-1, 1))
    #axs[0,0].set_yscale("symlog", linthresh=1e-7)


    axs[0,0].set_xscale("log")
    axs[0,0].set_yscale("log")
    axs[0,0].set_xlim((1.0, 1e4))
    #axs[0,0].set_ylim((1e-8, 1e8))
    axs[0,0].set_ylim((1e-8, 1e2))


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

        #print('e', e)
        #print('j', j)
        #print('b', b)

        lw = 0.5

        ex_flt = savgol_filter(ex, 1000, 1)

        #axs[0,0].plot(hh, ey, lw=lw, linestyle='solid', color='C1', alpha=0.8)
        #axs[0,0].plot(hh, ez, lw=lw, linestyle='solid', color='C2', alpha=0.8)

        #axs[0,0].plot(hh, bx-bx[0], lw=lw, linestyle='dashed', color='C0', alpha=0.8)
        #axs[0,0].plot(hh, by, lw=lw, linestyle='dashed', color='C1', alpha=0.8)
        #axs[0,0].plot(hh, bz, lw=lw, linestyle='dashed', color='C2', alpha=0.8)

        #axs[0,0].plot(hh, ex_flt, lw=lw, linestyle='dashed', color='C0', alpha=0.8)
        #axs[0,0].plot(hh, ex, lw=lw, linestyle='solid', color='C0', alpha=0.8)


        #--------------------------------------------------
        # spec

        nx = len(ex)
        dx = 1.0/Lh
        ks = np.fft.rfftfreq(nx, d=dx)

        # fourier spec of EM fields
        #ex_s  = ks*np.abs( np.fft.rfft(ex) )
        #ey_s  = ks*np.abs( np.fft.rfft(ey) )
        #ez_s  = ks*np.abs( np.fft.rfft(ez) )
        #ex_s2 = ks*np.abs( np.fft.rfft(ex_flt) )

        # energy in fields
        #ex_s  = np.abs( np.fft.rfft(ex**2) )
        #ey_s  = np.abs( np.fft.rfft(ey**2) )*1e6
        #ez_s  = np.abs( np.fft.rfft(ez**2) )
        #ex_s2 = np.abs( np.fft.rfft(ex_flt**2) )

        # poynting flux
        #ex_s  = np.abs( np.fft.rfft(ex*bx) )
        #ey_s  = np.abs( np.fft.rfft(ey*bx) ) #*1e6
        #ez_s  = np.abs( np.fft.rfft(ez*bx) )
        #ex_s2 = np.abs( np.fft.rfft(ex_flt*bx) )


        hwin = np.hamming(len(ex)) # add Hamming window to the data since this is not a periodic function in x

        # E field energy
        ex_s  = np.fft.rfft(hwin*ex)**2
        ey_s  = np.fft.rfft(hwin*ey)**2
        ez_s  = np.fft.rfft(hwin*ez)**2
        ex_s2 = np.fft.rfft(hwin*ex_flt)**2

        #axs[0,0].plot(ks, np.abs(ey_s), lw=lw, linestyle='solid', color='C1', alpha=0.8)
        #axs[0,0].plot(ks, np.abs(ez_s), lw=lw, linestyle='solid', color='C2', alpha=0.8)

        #axs[0,0].plot(ks, np.abs(ey_s)*1e6, lw=lw, linestyle='solid', color='C1', alpha=0.8)
        #axs[0,0].plot(ks, np.abs(ez_s)*1e6, lw=lw, linestyle='solid', color='C2', alpha=0.8)
        ##axs[1,0].plot(ks, np.abs(ex_s2), lw=lw,linestyle='solid', color='C3', alpha=0.8)

        axs[0,0].plot(ks, np.abs(ex_s), lw=lw, linestyle='solid', color='C0', alpha=0.8)

    #--------------------------------------------------
    # time stamp

    #print('r_pc:', conf.rad_pcap)
    print('t:', args.lap/conf.t_norm)
    stitle = r"$t/t_\mathrm{esc}$ = " + "{:3.1f}".format(args.lap/conf.t_norm)
    axs[0,0].set_title(stitle, fontsize=10)

    #--------------------------------------------------
    # manual slopes

    ks = np.array([1e2, 1e3])
    ek = ks**-1.0
    ek[:] *= 5e0/ek[0]
    axs[0,0].plot(ks, ek, color="k", linestyle="dashed", lw=1.5)

    ks = np.array([2e3, 6e3])
    ek = ks**-5.0
    ek[:] *= 2e-1/ek[0]
    axs[0,0].plot(ks, ek, color="k", linestyle="dashed", lw=1.5)

    #--------------------------------------------------
    #manual peak

    tacc = conf.cfl/conf.vrot/conf.bstar
    print("tacc", tacc)
    tesc = conf.t_norm
    print("tesc", tesc)
    print("tesc/tacc", tesc/tacc)
    print("tacc/tesc", tacc/tesc)

    M_low = 1e1
    M_high = 1e2
    g3_low = 1e-2
    g3_high = 1e-1

    wosc_low = ( M_low * (tesc/tacc)*g3_low )**0.5
    wosc_high = ( M_high * (tesc/tacc)*g3_high )**0.5   

    #Note: this is actually kosc*h_pcap = (wosc/cfl)*h_pcap = wosc*tesc    
    print("wosc_low", wosc_low)
    print("wosc_high", wosc_high)

    ks_low = wosc_low*np.ones(2)
    ks_high = wosc_high*np.ones(2)

    print(( M_low * 1.0e9*g3_low )**0.5 )
    ks_low_msp_phys = ( M_low * 1.0e9*g3_low )**0.5 * np.ones(2)
    ks_high_msp_phys = ( M_high * 1.0e9*g3_high )**0.5 * np.ones(2)

    ek = [5e-1, 5e1]
    axs[0,0].plot(ks_low, ek, color="r", linestyle="dashed", lw=1.5)
    axs[0,0].plot(ks_high, ek, color="r", linestyle="dashed", lw=1.5)

    #axs[0,0].plot(ks_low_msp_phys, ek, color="b", linestyle="dashed", lw=1.5)
    #axs[0,0].plot(ks_high_msp_phys, ek, color="b", linestyle="dashed", lw=1.5)

    #--------------------------------------------------
    axleft    = 0.18
    axbottom  = 0.12
    axright   = 0.97
    axtop     = 0.93

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

    #fname = fdir + 'fig_kspec_PRL_' + slap + '.png' 
    fname = fdir + 'fig_kspec_PRL_' + slap + '.pdf'     
    plt.savefig(fname, dpi=300,bbox_inches='tight')



