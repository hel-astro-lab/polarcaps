import numpy as np
import sys, os
import h5py as h5py

import matplotlib

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm
from matplotlib.cm import get_cmap
from matplotlib import colorbar

from numpy import pi, log, log10, exp, sin, cos, tanh

import pytools


#-------------------------------------------------- 
if __name__ == "__main__":

    #--------------------------------------------------
    fig = plt.figure(1, figsize=(3.25, 10.0)) # single figure
    #fig = plt.figure(1, figsize=(7.0,  6.0)) # two-column figure
    #fig = plt.figure(1, figsize=(9.0,  9.0)) # two-column figure

    plt.rc('xtick', top   = True)
    plt.rc('ytick', right = True)

    plt.rc('font',  family='serif',)
    plt.rc('text',  usetex=False)

    plt.rc('xtick', labelsize=7)
    plt.rc('ytick', labelsize=7)
    plt.rc('axes',  labelsize=8)
    plt.rc('legend',  handlelength=4.0)

    nrow_fig = 5
    ncol_fig = 1

    gs = plt.GridSpec(nrow_fig, ncol_fig)

    gs.update(wspace = 0.25)
    gs.update(hspace = 0.30)

    axs = np.empty( (nrow_fig,ncol_fig), dtype=object)

    for j in range(ncol_fig):
        for i in range(nrow_fig):
            axs[i,j] = plt.subplot(gs[i,j])
            axs[i,j].minorticks_on()

    axs[0,0].set_xlabel(r"$\chi_\pm$")
    axs[1,0].set_xlabel(r"$\chi_\pm$")

    axs[0,0].set_ylabel(r"$a$")
    axs[1,0].set_ylabel(r"$n$")

    axs[2,0].set_xlabel(r"$\chi_\gamma/\chi_e$")
    axs[2,0].set_ylabel(r"$\zeta$")

    axs[3,0].set_xlabel(r"$x/\gamma_e$")
    axs[3,0].set_ylabel(r"$dP/dx$")

    axs[4,0].set_xlabel(r"$x/\gamma_e$")
    axs[4,0].set_ylabel(r"$dP/dx$")

    for j in range(ncol_fig):
        for i in range(nrow_fig):
            axs[i,j].set_xscale('log')
            axs[i,j].set_yscale('log')
            axs[i,j].set_xlim((1e-4, 1e3))

    #axs[0,0].set_ylim((-2, 6))
    #axs[1,0].set_ylim((-0.03, 0.03))
    #axs[1,0].set_ylim((-1.0, 1.0))

    axs[0,0].set_ylim((1e-4, 1e4))
    axs[1,0].set_ylim((1e-1, 1e1))
    axs[2,0].set_ylim((1e-3, 1.0))

    axs[2,0].set_xlim((1e-8, 1e3))

    #--------------------------------------------------

    def a_fit(x):
        r  = 7.635/x**0.996
        r *= 0.5 - 0.4995*np.tanh(0.1271*(x - 21.93))

        #r  = 1.5/x**0.5
        return r

    def n_fit(x):
        return 0.313*(1 + 3.2*np.exp(-x**0.4332))

        #return 0.3333*np.ones_like(x)
        #return 0.313*(1 + 3.2)*np.ones_like(x) 


    x = np.logspace(-4, 3, 100)
    a = a_fit(x)
    n = n_fit(x)

    axs[0,0].plot(x, a)

    axs[1,0].plot(x, n)

    def xg_fit(x, z):
        a = a_fit(x)
        n = n_fit(x)
        return  x*(z/(1 + a*(1-z)))**(1/n)

    z = np.linspace(0, 1, 100)

    for xe in np.logspace(-3, 3, 20):
        xg = xg_fit(xe, z)

        axs[2,0].plot(xg/xe, z)


    # spectrum
    axs[3,0].set_xlim((1e-6, 1e1))
    axs[3,0].set_yscale('linear')

    axs[4,0].set_xlim((1e-6, 1e1))
    axs[4,0].set_ylim((1e-6, 1e-1))


    N = int(1e6)
    gam = 1000
    zs = np.random.rand(N) # random numbers


    # Maxwell Juttner initial beam (instead of cold beam with gamma = gam)
    #theta = 1 # 2.5e-3
    #gs = np.zeros(N) # random numbers
    #for i in range(N):
    #    ux,uy,uz, u = pytools.sample_boosted_maxwellian(theta, gam, dims=3)
    #    gs[i] = np.sqrt(1 + u*u)

    def get_col(x):
        if x == 0.01: return 'C0'
        if x == 0.1:  return 'C1'
        if x == 1.0:  return 'C2'
        if x ==10.0:  return 'C3'


    #for Xe in [0.1, 1.0]:
    for Xe in [0.01, 0.1, 1.0, 10]:
        #Xe = 0.1 # electron chi

        col = get_col(Xe)

        Xg = xg_fit(Xe, zs) # photon chi
        xs = gam*Xg/Xe # photon energy; e = mc^2 gam_e \chi_x/\chi_e

        # spectrum of x/\gamma 
        xbin = np.logspace(-6, 3, 200)
        hist, bins = np.histogram(xs/gam, bins=xbin)

        hist = hist/N # normalize by number of particles

        #print(Xg)
        #print(xs/gam)
        #print(hist)

        xxs = xbin[1:] # photon energy array (in correct length for the histogram)

        # note that we bin using log x so histogram is d(N)/d(log(x)) = x dN/dx = d(N*x)/dx \propto dP/dx
        axs[3,0].plot(xxs, hist, drawstyle='steps-pre', alpha=1.0, lw=1.5, linestyle='solid', color=col)
        axs[4,0].plot(xxs, hist, drawstyle='steps-pre', alpha=1.0, lw=1.5, linestyle='solid', color=col)


    #--------------------------------------------------
    # classical synchrotron spectrum
    if False:
        for Xe in [0.01, 0.1, 1.0, 10]:

            col = get_col(Xe)

            b = Xe/gam
            xcrit = 1.5*b*gam**2
            y = xbin/xcrit #*gam #/0.01/3
            dPdx_cl = 0.05*np.exp(-y)*y**(1/3)
        
            axs[3,0].plot(xbin/gam, dPdx_cl, linestyle='dotted', lw=1.5, color=col)
            axs[4,0].plot(xbin/gam, dPdx_cl, linestyle='dotted', lw=1.5, color=col)
        
    

    #--------------------------------------------------
    # plot Smilei reference data 
    # from https://smileipic.github.io/Smilei/Understand/radiation_loss.html#equation-radiatedpowerspectrum (fig28)
    if False:
        for (fname,Xe) in [
                ("smilei_syn_dPdx_x0.1.csv", 0.1),
                ("smilei_syn_dPdx_x1.csv", 1.0),
                ]:
            d = np.loadtxt(fname, delimiter=",")

            col = get_col(Xe)

            xg = d[:, 0]  
            dP = d[:, 1]*0.5 # normalize with arbitrary units

            axs[3,0].plot(xg, dP, linestyle='dashed', lw=0.8, color=col)
            axs[4,0].plot(xg, dP, linestyle='dashed', lw=0.8, color=col)





    #--------------------------------------------------
    axleft    = 0.17
    axbottom  = 0.10
    axright   = 0.97
    axtop     = 0.98

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

    fname = 'syn_revised_fitting_formula.pdf' 
    plt.savefig(fname)

    #fname = 'pml.png' 
    #plt.savefig(fname, dpi=300)
















