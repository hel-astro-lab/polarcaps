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
    fig = plt.figure(1, figsize=(2.2, 8.0)) # single figure

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
        nrow_fig = 6
        ncol_fig = 1

        gs = plt.GridSpec(nrow_fig, ncol_fig)

        gs.update(wspace = 0.2)
        gs.update(hspace = 0.2)

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
    axs[2,0].set_ylabel(r"$x$ ($m_e c^2$)")
    axs[3,0].set_ylabel(r"$E/E_\mathrm{rot}$")
    axs[4,0].set_ylabel(r"$j/j_m$")
    axs[5,0].set_ylabel(r"$B/B_0$")


    hmin = -0.1
    hmax = 1.5 #7.0
    for j in range(ncol_fig):
        for i in range(nrow_fig):
            axs[i,j].set_xlim((hmin, hmax))

    axs[3,0].set_ylim((-1, 1))
    axs[3,0].set_ylim((-1.2, 1.2))
    axs[4,0].set_ylim((-1.2, 3.2))
    axs[5,0].set_ylim((0, 1.2))

    axs[5,0].set_xlabel(r"$h/h_\mathrm{pc}$")

    for j in range(ncol_fig):
        for i in range(nrow_fig):
            axs[i,j].fill_between([-0.5, 0.0], -100, 100, color='gray', alpha=0.4, edgecolor=None) 


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

        print('hhlims', hhlims)
        print(toolset.pxlims[1],  toolset.pxlims[0] )
        px_log_extent = toolset.pxlims[1] - toolset.pxlims[0] 
        print('px_log_extent:', px_log_extent)

        for i in range(0,3): 
            axs[i,0].set_ylim((-px_log_extent, px_log_extent))


        print(toolset.xxlims[1],  toolset.xxlims[0] )
        xx_log_extent = toolset.xxlims[1] - toolset.xxlims[0] 
        print('xx_log_extent:', xx_log_extent)

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
                         )

        # height vs ene; photons
        imph = axs[2,0].imshow( np.zeros((toolset.Nhist, 2*toolset.Nhist)),
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
        imep.set_data(np.log10( hep.T)) 
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
                    xps = xp - 2
                    return f"$10^{int(xps)}$"
                if xp == 0:
                    return "$10^{-2}$"
                else:
                    xps = -xp - 2
                    return f"-$10^{int(xps)}$"

    yticks = [-6, -2, 2, 6]
    axs[0,0].set_yticks(yticks)
    axs[1,0].set_yticks(yticks)
    axs[2,0].set_yticks(yticks)

    majorformatter = CustomScalarFormatter()
    axs[0,0].yaxis.set_major_formatter(majorformatter)
    axs[1,0].yaxis.set_major_formatter(majorformatter)
    axs[2,0].yaxis.set_major_formatter(majorformatter)


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
        
            axs[1,0].axhline(y=xvals[i1], linestyle='solid', color='C0', lw=0.4)
            axs[1,0].axhline(y=xvals[i2], linestyle='solid', color='C0', lw=0.4)

        for xref in [conf.xsyn]:
            i1 = find_arg_nearest(xarr, -gref)
            i2 = find_arg_nearest(xarr, +gref)

            print('i1:', i1, 'x', xarr[i1], 'xref', -xref)
            print('i2:', i2, 'x', xarr[i2], 'xref', +xref)

            axs[2,0].axhline(y=xvals[i1], linestyle='solid', color='C0', lw=0.4)
            axs[2,0].axhline(y=xvals[i2], linestyle='solid', color='C0', lw=0.4)



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
        
        axs[3,0].plot(hh, ex, lw=0.8, linestyle='solid', color='C0', alpha=0.8)
        axs[3,0].plot(hh, ey, lw=0.8, linestyle='solid', color='C1', alpha=0.8)
        axs[3,0].plot(hh, ez, lw=0.8, linestyle='solid', color='C2', alpha=0.8)

        axs[4,0].plot(hh, jx, lw=0.8, linestyle='solid', color='C0', alpha=0.8)
        axs[4,0].plot(hh, jy, lw=0.8, linestyle='solid', color='C1', alpha=0.8)
        axs[4,0].plot(hh, jz, lw=0.8, linestyle='solid', color='C2', alpha=0.8)

        axs[5,0].plot(hh, bx, lw=0.8, linestyle='solid', color='C0', alpha=0.8)


    #--------------------------------------------------
    # time stamp

    #print('r_pc:', conf.rad_pcap)
    print('t:', args.lap/conf.t_norm)
    stitle = r"$t c/R_\mathrm{pc}$ = " + "{:3.1f}".format(args.lap/conf.t_norm)
    axs[0,0].set_title(stitle, fontsize=10)


    #--------------------------------------------------
    axleft    = 0.22
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


    slap = str(args.lap).rjust(5, '0')

    #fname = fdir + 'fig_casc_' + slap + '.pdf' 
    #plt.savefig(fname)

    fname = fdir + 'fig_gap_' + slap + '.png' 
    plt.savefig(fname, dpi=300)



