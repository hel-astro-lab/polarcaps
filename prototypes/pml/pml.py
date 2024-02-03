import numpy as np
import sys, os

import matplotlib

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import get_cmap
from matplotlib import colorbar

from numpy import pi, log, log10, exp, sin, cos 


#-------------------------------------------------- 
if __name__ == "__main__":

    #--------------------------------------------------
    fig = plt.figure(1, figsize=(3.25, 2.5)) # single figure
    #fig = plt.figure(1, figsize=(7.0,  5.5)) # two-column figure

    plt.rc('xtick', top   = True)
    plt.rc('ytick', right = True)

    plt.rc('font',  family='serif',)
    plt.rc('text',  usetex=False)

    plt.rc('xtick', labelsize=7)
    plt.rc('ytick', labelsize=7)
    plt.rc('axes',  labelsize=8)
    plt.rc('legend',  handlelength=4.0)

    nrow_fig = 1
    ncol_fig = 1

    gs = plt.GridSpec(nrow_fig, ncol_fig)

    gs.update(wspace = 0.25)
    gs.update(hspace = 0.35)

    axs = np.empty( (nrow_fig,ncol_fig), dtype=object)

    for j in range(ncol_fig):
        for i in range(nrow_fig):
            axs[i,j] = plt.subplot(gs[i,j])
            axs[i,j].minorticks_on()

    axs[0,0].set_xlabel(r"$x$")
    axs[0,0].set_ylabel(r"$y$")

    #axs[0,0].set_xlabel(r"$x~(\mathrm{cm})$")
    #axs[0,0].set_ylabel(r"$y_\mathrm{c}~(\mathrm{cm}\,\mathrm{s}^{-1})$")

    for j in range(ncol_fig):
        for i in range(nrow_fig):
            axs[i,j].set_xlim((0, 1.2))
            #axs[i,j].set_ylim((1, 1e3))

    #axs[0,0].set_xscale('log')
    #axs[0,0].set_yscale('log')

    #--------------------------------------------------

    rad_lim = 0.6

    z = np.linspace(0, 1.2, 100)

    lam = ( (z - rad_lim)/(1 - rad_lim) )**3
    for i in range(100):
        lam[i] = -10*np.minimum(lam[i], 1.0)

        if z[i] < rad_lim:
            lam[i] = 0.0

    #axs[0,0].plot(z, lam)

    bx = 1.0*np.ones(100)
    dbx = 0.5

    new_bx = np.zeros(100)
    for i in range(100):
        new_bx[i] = ( bx[i]*(1+lam[i]) + dbx )/(1 - lam[i])

        #new_bx[i] = bx[i] + (dbx  + bx[i]*lam[i])
        #new_bx[i] = bx[i] + (1+lam[i])*dbx


    axs[0,0].plot(z, bx, color='C0')
    axs[0,0].plot(z, new_bx, color='C1')






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


    #slap = str(args_cli.lap).rjust(4, '0')

    fname = 'pml.pdf' 
    plt.savefig(fname)

    fname = 'pml.png' 
    plt.savefig(fname, dpi=300)



