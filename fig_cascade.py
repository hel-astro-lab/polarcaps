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
from init_problem import Configuration_Turbulence as Configuration
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
    fig = plt.figure(1, figsize=(3.25, 3.5)) # single figure
    #fig = plt.figure(1, figsize=(7.0,  5.5)) # two-column figure

    plt.rc('xtick', top   = True)
    plt.rc('ytick', right = True)

    plt.rc('font',  family='serif',)
    plt.rc('text',  usetex=False)

    plt.rc('xtick', labelsize=7)
    plt.rc('ytick', labelsize=7)
    plt.rc('axes',  labelsize=8)
    plt.rc('legend',  handlelength=4.0)

    nrow_fig = 3
    ncol_fig = 1

    gs = plt.GridSpec(nrow_fig, ncol_fig)

    gs.update(wspace = 0.25)
    gs.update(hspace = 0.0)

    axs = np.empty( (nrow_fig,ncol_fig), dtype=object)

    for j in range(ncol_fig):
        for i in range(nrow_fig):
            axs[i,j] = plt.subplot(gs[i,j])
            axs[i,j].minorticks_on()


    axs[0,0].set_xlabel(r"$h/R_{\mathrm{pc}}$")
    axs[0,0].set_ylabel(r"$p$")

    hmin =  0.0
    hmax =  2.5

    for j in range(ncol_fig):
        for i in range(nrow_fig):
            axs[i,j].set_xlim((hmin, hmax))

    #axs[i,j].set_ylim((1, 1e3))

    #axs[0,0].set_xscale('log')
    #axs[0,0].set_yscale('log')

    #tmin = 0.0
    #tmax = 10.0
    #norm = matplotlib.colors.Normalize(vmin=tmin, vmax=tmax)
    #cmap = matplotlib.colormaps['turbo_r']

    #--------------------------------------------------

    # read from output file
    fname = conf.outdir + '/qed_{}.h5'.format(str(args.lap))
    f5 = h5.File(fname,'r')

    hem = f5['h2_ene_e-'][()]
    hep = f5['h2_ene_e+'][()]
    hph = f5['h2_ene_ph'][()]

    f5.close()

    #--------------------------------------------------

    if conf.twoD:
        cenx   = conf.Lx//2 + 0.5
        ceny   = -conf.rad_star + conf.rad_curv_shift
        cenz   = 0 

        Lh = conf.Ly
    elif conf.threeD:
        cenx   = conf.Lx//2 + 0.5
        #ceny   = conf.Ly//2 + 0.5
        ceny   = -conf.rad_star + conf.rad_curv_shift
        Lh = conf.Lz

    hhlims = [1,1]
    hhlims[0] = 0 #(0 - ceny)/conf.rad_pcap
    hhlims[1] = Lh/conf.rad_pcap

    print('hhlims', hhlims)

    imem = axs[0,0].imshow( np.zeros((toolset.Nhist, 2*toolset.Nhist)),
                     extent=[hhlims[0], hhlims[1], 
                             -toolset.Nhist, toolset.Nhist], # TODO
                           #toolset.pxlims[0], toolset.pxlims[1]],
                     origin='lower',
                     cmap='turbo',
                     aspect='auto',
                     interpolation='nearest',
                     vmin=np.log10(toolset.pylims[0]),
                     vmax=np.log10(toolset.pylims[1]),)

    imep = axs[1,0].imshow( np.zeros((toolset.Nhist, 2*toolset.Nhist)),
                     extent=[hhlims[0], hhlims[1], 
                             -toolset.Nhist, toolset.Nhist], # TODO
                             #toolset.pxlims[0], toolset.pxlims[1]],
                     origin='lower',
                     cmap='turbo',
                     aspect='auto',
                     interpolation='nearest',
                     vmin=np.log10(toolset.pylims[0]),
                     vmax=np.log10(toolset.pylims[1]),)

    # height vs ene; photons
    imph = axs[2,0].imshow( np.zeros((toolset.Nhist, 2*toolset.Nhist)),
                     extent=[hhlims[0], hhlims[1], 
                             -toolset.Nhist, toolset.Nhist], # TODO
                           #toolset.xxlims[0], toolset.xxlims[1]],
                         origin='lower',
                         cmap='turbo',
                         aspect='auto',
                         interpolation='nearest',
                         vmin=np.log10( toolset.xylims[0]),
                         vmax=np.log10( toolset.xylims[1]),)



    imem.set_data(np.log10( hem.T)) 
    imep.set_data(np.log10( hep.T)) 
    imph.set_data(np.log10( hph.T)) 


    #--------------------------------------------------
    axleft    = 0.15
    axbottom  = 0.15
    axright   = 0.97
    axtop     = 0.82

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

    fig.subplots_adjust(left=axleft, bottom=axbottom, right=axright, top=axtop)


    slap = str(args.lap).rjust(5, '0')

    #fname = fdir + 'fig_casc_' + slap + '.pdf' 
    #plt.savefig(fname)

    fname = fdir + 'fig_casc_' + slap + '.png' 
    plt.savefig(fname, dpi=300)



