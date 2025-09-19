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
#warnings.filterwarnings('ignore')

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
    #fig = plt.figure(1, figsize=(2.0, 5.0)) # narrow figure

    fig = plt.figure(1, figsize=(3.25, 5.0)) # single figure
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
        nrow_fig = 1
        ncol_fig = 1

        gs = plt.GridSpec(nrow_fig, ncol_fig)

        #gs.update(wspace = 0.2)
        #gs.update(hspace = 0.12)

        axs = np.empty( (nrow_fig,ncol_fig), dtype=object)

        for j in range(ncol_fig):
            for i in range(nrow_fig):
                axs[i,j] = plt.subplot(gs[i,j])
                axs[i,j].minorticks_on()

                #axs[i,j].set_yscale('log')
                #if i <= nrow_fig-2:
                #    axs[i,j].tick_params(labelbottom=False)
            #axs[0,0].tick_params(labeltop=True)

    axs[0,0].set_ylabel(r"$t/t_\mathrm{esc}$")
    axs[0,0].set_xlabel(r"$x/H_\mathrm{gap}$")



    # box limits
    Lh = conf.Lx
    xmin = 0
    xmax = Lh/conf.rad_pcap


    laps = np.array( list( range(0, conf.Nt, conf.interval) ) )
    tt = laps/conf.t_norm

    nx = int( conf.Nx*conf.NxMesh/conf.stride )
    ny = int( len(laps) )
    data = np.zeros((ny,nx))

    tmax = 0.0 # monitor maxium time 

    for i, lap in enumerate(laps):

        # read from output file
        fname = conf.outdir + '/flds_{}.h5'.format(str(lap))
        if not(os.path.isfile(fname)): continue

        tmax = tt[i] # update max time

        print(lap)

        f5 = h5.File(fname,'r')

        nn = pytools.read_h5_array(f5, 'rho')/conf.p_norm

        ex = pytools.read_h5_array(f5, 'ex')/conf.e_norm
        ey = pytools.read_h5_array(f5, 'ey')/conf.e_norm
        ez = pytools.read_h5_array(f5, 'ez')/conf.e_norm

        jx = pytools.read_h5_array(f5, 'jx')/conf.j_norm
        jy = pytools.read_h5_array(f5, 'jy')/conf.j_norm
        jz = pytools.read_h5_array(f5, 'jz')/conf.j_norm

        bx = pytools.read_h5_array(f5, 'bx')/conf.b_norm
        by = pytools.read_h5_array(f5, 'by')/conf.b_norm
        bz = pytools.read_h5_array(f5, 'bz')/conf.b_norm

        f5.close()

        #print('shape ex', np.shape(ex))

        # reduce dimensions
        if conf.oneD:
            nn = np.mean(nn, axis=(1,2))

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

        #e = ex + ey + ez
        #b = bx + by + bz
        #j = jx + jy + jz
        #hh = np.linspace(xmin, xmax, Lh)

        if args.var == 'ex':
            data[i, :] = ex
        if args.var == 'ex2':
            data[i, :] = ex**2
        if args.var == 'dens':
            data[i, :] = nn

    #--------------------------------------------------
    # done reading files
    print('-------------------------------------------------- ')
    print('data min/max:', np.min(data), np.max(data))
    print('data mean:', np.mean(data) )

    # scale plot to visible area
    axs[0,0].set_xlim((0.0, 1.25))
    axs[0,0].set_ylim((0.0, tmax))

    #--------------------------------------------------
    # now plot

    ymin = tt[0]
    ymax = tt[-1]

    #zmin = 0.0
    #zmax = 1.0
    #norm = matplotlib.colors.Normalize(vmin=-1.0, vmax=1.0)
    #norm = matplotlib.colors.SymLogNorm(1e-4, linscale=1, vmin=-1.0, vmax=1.0)
    #cmap = matplotlib.colormaps['RdBu']

    #norm = matplotlib.colors.LogNorm(vmin=1e-7, vmax=1)


    if args.var == 'ex':
        norm = matplotlib.colors.SymLogNorm(1e-4, linscale=1, vmin=-1.0, vmax=1.0)
        cmap = matplotlib.colormaps['RdBu']
        clabel = r'$E_x$'
    if args.var == 'ex2':
        data[i, :] = ex**2
        norm = matplotlib.colors.LogNorm(vmin=1e-7, vmax=1e0)
        cmap = matplotlib.colormaps['magma']
        clabel = r'$E_x^2$'
    if args.var == 'dens':
        norm = matplotlib.colors.LogNorm(vmin=1e0, vmax=1e4)
        cmap = matplotlib.colormaps['viridis']
        clabel = r'$m_\pm$'

    extent = [ xmin, xmax, ymin, ymax ]

    im = axs[0,0].imshow(
            data,
            extent=extent,
            origin='lower',
            interpolation='nearest',
            cmap = cmap,
            aspect="auto",
            norm = norm,
            )

    #--------------------------------------------------
    # save

    axleft    = 0.12
    axbottom  = 0.07
    axright   = 0.99
    axtop     = 0.88

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
        cb1.set_label(clabel)

    pos = axs[0,0].get_position()
    print('ax pos:', pos)
    fig.subplots_adjust(left=axleft, bottom=axbottom, right=axright, top=axtop)

    #slap = str(args.lap).rjust(5, '0')
    #fname = fdir + 'fig_casc_' + slap + '.pdf' 
    #plt.savefig(fname)

    fname = fdir + 'fig_' + args.var + '_waterfall.pdf'
    plt.savefig(fname, dpi=300)

