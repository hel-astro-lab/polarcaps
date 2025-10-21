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
    fig = plt.figure(1, figsize=(3.25, 4.5)) # single figure

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
        nrow_fig = 2
        ncol_fig = 1

        gs = plt.GridSpec(nrow_fig, ncol_fig)

        gs.update(wspace = 0.2)
        gs.update(hspace = 0.3)

        axs = np.empty( (nrow_fig,ncol_fig), dtype=object)

        for j in range(ncol_fig):
            for i in range(nrow_fig):
                axs[i,j] = plt.subplot(gs[i,j])
                axs[i,j].minorticks_on()

                #axs[i,j].set_yscale('log')
                #if i <= nrow_fig-2:
                #    axs[i,j].tick_params(labelbottom=False)

                axs[i,j].set_xlim((0.0, 5.0))
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
    hmin = int( 2 )
    hmax = int( 20 )

    # limits for the height-energy histogram
    hs = np.linspace(0, Lh, toolset.Nhist)

    for i in range(len(hs)):
        if hs[i] > hmin: break
        hmin2 = i

    for i in range(len(hs)):
        if hs[i] > hmax: break
        hmax2 = i

    print("hs", hs)
    print("hmin/hmax", hmin, hmax)
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
    #print(xlow_mask2)

    #--------------------------------------------------
    def integrate(xs, ys):
        # trapedzoidal
        #dxs = np.diff(xs)
        dxs = xs[1] - xs[0] # assume uniform grid
        return np.sum(ys*dxs)


    tmin = 0.0
    tmax = 2.0
    tspan = tmax - tmin

    norm = matplotlib.colors.Normalize(vmin=0, vmax=tspan)
    cmap = matplotlib.colormaps['turbo_r']

    merge = 6
    i = 0



    # total time-integrated spec
    N = toolset.Nhist
    he_tot = np.zeros(N)
    hp_tot = np.zeros(N)
    hx_tot = np.zeros(N)
    

    for lap in laps:

        if not(tmin < lap/conf.t_norm < tmax):
            continue

        print("lap", lap, "t", lap/conf.t_norm)

        col = cmap(norm( lap/conf.t_norm - tmin ))

        #--------------------------------------------------
        # particle spectra

        # read from output file
        fname = conf.outdir + '/qed_{}.h5'.format(str(lap))
        f5 = h5.File(fname,'r')

        hem = f5['h2_ene_e-'][()]
        hep = f5['h2_ene_e+'][()]
        hph = f5['h2_ene_ph'][()]

        f5.close()

        #--------------------------------------------------

        #--------------------------------------------------
        # field values

        ## read from output file
        #fname = conf.outdir + '/flds_{}.h5'.format(str(lap))
        #f5 = h5.File(fname,'r')

        #ex = pytools.read_h5_array(f5, 'ex')/conf.e_norm
        #ey = pytools.read_h5_array(f5, 'ey')/conf.e_norm
        #ez = pytools.read_h5_array(f5, 'ez')/conf.e_norm

        #jx = pytools.read_h5_array(f5, 'jx')/conf.j_norm
        #jy = pytools.read_h5_array(f5, 'jy')/conf.j_norm
        #jz = pytools.read_h5_array(f5, 'jz')/conf.j_norm

        #bx = pytools.read_h5_array(f5, 'bx')/conf.b_norm
        #by = pytools.read_h5_array(f5, 'by')/conf.b_norm
        #bz = pytools.read_h5_array(f5, 'bz')/conf.b_norm

        ##print('shape ex', np.shape(ex))

        ## reduce dimensions
        #if conf.oneD:
        #    ex = np.mean(ex, axis=(1,2))
        #    ey = np.mean(ey, axis=(1,2))
        #    ez = np.mean(ez, axis=(1,2))

        #    jx = np.mean(jx, axis=(1,2))
        #    jy = np.mean(jy, axis=(1,2))
        #    jz = np.mean(jz, axis=(1,2))

        #    bx = np.mean(bx, axis=(1,2))
        #    by = np.mean(by, axis=(1,2))
        #    bz = np.mean(bz, axis=(1,2))
        #else:
        #    print('TODO')

        #e = ex + ey + ez
        #b = bx + by + bz
        #j = jx + jy + jz

        ## compensate for external antenna
        #if lap > conf.rad_pcap/conf.cfl:
        #    jx[:] += 1.0

        #print('ex', ex[hmin])
        #print('j', j)
        #print('b', b)


        #--------------------------------------------------
        #particles
        # he and hp store particles as histograms of dn/dlnp

        # histogram into units of n_GJ
        n_units = toolset.N_box/toolset.N_wgt # de-unitize what we have in the qed_toolset
        n_units *= 1/(hmax2-hmin2) # normalize by area into xx per cell
        n_units *= 1/conf.ppc # normalize to n_GJ

        # in units of dN/dlnp 
        he = np.sum(hem[hmin2:hmax2,:], axis=0)*n_units #*zs2
        #he = he[N:] + np.flip(he[:N]) # up and down
        he = np.flip(he[:N]) # down

        #axs[0,0].plot(lnzs2, n_units*he, color=col)
        #axs[0,0].plot(lnzs, n_units*he, color=col, lw=0.8)

        #--------------------------------------------------
        hp = np.sum(hep[hmin2:hmax2,:], axis=0)*n_units #*zs2
        #hp = hp[N:] + np.flip(hp[:N]) # up and down
        hp = np.flip(hp[:N]) # down

        #numpp = integrate(lnzs2, hp)*n_units # integrate and de-unitize
        #axs[1,0].plot(lnzs2, n_units*he, color=col)
        #axs[1,0].plot(lnzs, n_units*hp, color=col, lw=0.8)

        #--------------------------------------------------
        # photons
        # hph stores photons as histograms of x^2 d(n_x x w)/dx = x^2 d(n_x w)/dlnx

        # conversion factor into units of n_GJ
        nx_units = toolset.N_box*toolset.N_time/toolset.N_wgt
        nx_units *= 1/(hmax2-hmin2) # normalize by area into xx per cell
        nx_units *= 1/conf.ppc # normalize to n_GJ

        hx = np.sum(hph[hmin2:hmax2,:], axis=0)*nx_units #/xs2**2
        #hx = hx[N:] + np.flip(hx[:N]) # up and down
        hx = np.flip(hx[:N]) #"down

        #numpx = integrate(lnxs2, xlow_mask2*hx/xs2**2)*nx_units # integrate and de-unitize
        #axs[2,0].plot(lnxs2, nx_units*hx, color=col)
        #axs[2,0].plot(lnxs, nx_units*hx, color=col, lw=0.8)


        # sum for total spec
        he_int[:] += he
        hp_int[:] += hp
        hx_int[:] += hx

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
            #axs[0,0].plot(lnzs, heM, color=col, lw=0.8)
            #axs[1,0].plot(lnzs, hpM, color=col, lw=0.8)

            axs[0,0].plot(lnzs, heM + hpM, color=col, lw=0.8)
            axs[1,0].plot(lnxs, hxM,       color=col, lw=0.8)
            i = 0


    # plot total spec as well
    axs[0,0].plot(lnzs, he_int + hp_int, color=col, lw=1.5)
    axs[1,0].plot(lnxs, hx_int,          color=col, lw=1.5)


    #-------------------------------------------------- 
    axs[0,0].set_xlabel(r"$\log_{10} p$")
    axs[1,0].set_xlabel(r"$\log_{10} x$")

    axs[0,0].set_ylabel(r"$  \mathrm{d} m_\pm/ \mathrm{d} \, \log p $")
    axs[1,0].set_ylabel(r"$x \mathrm{d}( x m_x )/ \mathrm{d} \, \log x $")

    axs[0,0].set_yscale("log")
    axs[1,0].set_yscale("log")

    axs[0,0].set_ylim((1e-1, 1e3))
    axs[1,0].set_ylim((1e1,  1e6))


    #--------------------------------------------------
    axleft    = 0.17
    axbottom  = 0.08
    axright   = 0.97
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
        cb1.set_label(r'Time $(t-t_0)/t_\mathrm{esc}$')

    pos = axs[0,0].get_position()
    print('ax pos:', pos)
    fig.subplots_adjust(left=axleft, bottom=axbottom, right=axright, top=axtop)

    #slap = str(args.lap).rjust(5, '0')
    #fname = fdir + 'fig_casc_' + slap + '.pdf' 
    #plt.savefig(fname)

    fname = fdir + 'fig_return_prtcl_spec.pdf'
    plt.savefig(fname, dpi=300)



