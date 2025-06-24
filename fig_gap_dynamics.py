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
        nrow_fig = 3
        ncol_fig = 1

        gs = plt.GridSpec(nrow_fig, ncol_fig)

        gs.update(wspace = 0.2)
        gs.update(hspace = 0.05)

        axs = np.empty( (nrow_fig,ncol_fig), dtype=object)

        for j in range(ncol_fig):
            for i in range(nrow_fig):
                axs[i,j] = plt.subplot(gs[i,j])
                axs[i,j].minorticks_on()

                #axs[i,j].set_yscale('log')

                if i <= nrow_fig-2:
                    axs[i,j].tick_params(labelbottom=False)

                axs[i,j].set_xlim((0.0, 20.0))

            #axs[0,0].tick_params(labeltop=True)

    axs[0,0].set_ylabel(r"$p_-$ ($m_e c$)")
    axs[1,0].set_ylabel(r"$p_+$ ($m_e c$)")


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
    tt = laps/conf.t_norm

    hmin = int( 15 )
    hmax = int( hmin + conf.rad_pcap*0.5 )


    epars = []
    epeak = []

    jps   = []
    mpps   = []
    mpes   = []
    mpxs   = []

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
    lnzs2 = np.array( [-np.flip(lnzs), lnzs ] ).flatten()

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



    for lap in laps:
        #print("lap", lap, "t =", lap/conf.t_norm)

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

        # read from output file
        fname = conf.outdir + '/flds_{}.h5'.format(str(lap))
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

        #print('shape ex', np.shape(ex))

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

        # compensate for external antenna
        if lap > conf.rad_pcap/conf.cfl:
            jx[:] += 1.0

        #print('e', e)
        #print('j', j)
        #print('b', b)


        if False:
        #if lap == 12000:
            #hh = np.linspace(0, Lh, Lh)
            #axs[0,0].plot(hh, ex, lw=0.8, linestyle='solid', color='C0', alpha=0.8)
            #axs[0,0].plot( [hmin, hmax], [1.0, 1.0])

            print("shape hem", np.shape(hem))

            
            h = np.sum(hph[hmin2:hmax2,:], axis=0)
            axs[0,0].plot(lnxs2, xlow_mask2*h/xs2**2)

            axs[0,0].plot(lnxs, np.flip(h[:256])/xs**2)
            axs[0,0].plot(lnxs, h[256:]/xs**2)

            break

        #--------------------------------------------------
        #particles
        # he and hp store particles as histograms of dn/dlnp

        N = toolset.Nhist

        # histogram into units of n_GJ
        n_units = toolset.N_box/toolset.N_wgt # de-unitize what we have in the qed_toolset
        n_units *= 1/(hmax2-hmin2) # normalize by area into xx per cell
        n_units *= 1/conf.ppc # normalize to n_GJ

        # in units of dN/dlnp 
        he = np.sum(hem[hmin2:hmax2,:], axis=0)

        #DONE: verified for same result from both versions

        # v1 integrating -vx and +vx parts separately
        #numpe = ( integrate(lnzs, he[N:]) + integrate(lnzs, np.flip(he[:N])) )*n_units

        # integrating on one go
        numpe = integrate(lnzs2, he)*n_units # integrate and de-unitize
        #print(numpe, numpe2)

        hp = np.sum(hep[hmin2:hmax2,:], axis=0)
        numpp = integrate(lnzs2, hp)*n_units # integrate and de-unitize

        #--------------------------------------------------
        # photons
        # hph stores photons as histograms of x^2 d(n_x x w)/dx = x^2 d(n_x w)/dlnx

        # conversion factor into units of n_GJ
        nx_units = toolset.N_box*toolset.N_time/toolset.N_wgt
        nx_units *= 1/(hmax2-hmin2) # normalize by area into xx per cell
        nx_units *= 1/conf.ppc # normalize to n_GJ

        # NOTE: we also mask x < 2 values out from the count
        hx = np.sum(hph[hmin2:hmax2,:], axis=0)
        numpx = integrate(lnxs2, xlow_mask2*hx/xs2**2)*nx_units # integrate and de-unitize

        mpes.append(numpe)
        mpps.append(numpp)
        mpxs.append(numpx)

        epeak.append( ex[hmin] )
        epars.append( np.mean(ex[hmin:hmax]) )
        jps.append( np.mean(jx[hmin:hmax]) )


    #--------------------------------------------------

    #v1 with electrons and positrons
    #axs[0,0].plot(tt, mpes, color="C0", lw=0.9, linestyle="solid")
    #axs[0,0].plot(tt, mpps, color="C0", lw=0.9, linestyle="dashed")

    # v2 with pairs and photons
    mpairs = np.array(mpes) + np.array(mpps)
    axs[0,0].plot(tt, mpairs, color="C0", lw=0.9, linestyle="solid")
    axs[0,0].plot(tt, mpxs, color="C2", lw=0.9, linestyle="dashed")

    # v1 with field quantities
    axs[1,0].plot(tt, epeak, color="C1", lw=0.9, linestyle="solid")
    #axs[2,0].plot(tt, jps/mpairs,   color="C3", lw=0.9)

    #axs[2,0].plot(tt, jps,   color="C3", lw=0.9)
    axs[2,0].plot(tt, jps/mpairs,     color="C3", lw=0.9)

    #axs[2,0].plot(tt, np.array(epeak),   color="C4", lw=0.9)
    axs[2,0].plot(tt, 0.2*(np.array(epeak)   ),   color="C4", lw=0.9, linestyle="dotted")
    #axs[2,0].plot(tt, 0.2*(np.array(epeak)**2),   color="C5", lw=0.9, linestyle="dotted")

    #de_dt = -0.1*np.gradient(epeak, tt)
    #axs[2,0].plot(tt, de_dt, color="C5", lw=0.9)

    #dn_dt = 1e-3*np.gradient(mpairs, tt)
    #axs[2,0].plot(tt, dn_dt, color="C6", lw=0.9)


    axs[0,0].set_ylabel(r"$m_\pm$, $m_x$")
    axs[1,0].set_ylabel(r"$\varepsilon = E_\parallel/\beta_\mathrm{rot} B_\star$")
    #axs[2,0].set_ylabel(r"$j_\pm/j_m$")
    axs[2,0].set_ylabel(r"$\langle v \rangle = j_\pm/m_\pm$")

    axs[2,0].set_xlabel(r"$t/t_\mathrm{esc}$")

    axs[0,0].set_yscale("log")

    axs[0,0].set_ylim((0.1, 1e2))
    axs[1,0].set_ylim((-0.1, 1.1))
    axs[2,0].set_ylim((-0.3, 0.3))
    #axs[2,0].set_ylim((-1.0, 1.0))


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

    pos = axs[0,0].get_position()
    print('ax pos:', pos)
    fig.subplots_adjust(left=axleft, bottom=axbottom, right=axright, top=axtop)

    #slap = str(args.lap).rjust(5, '0')
    #fname = fdir + 'fig_casc_' + slap + '.pdf' 
    #plt.savefig(fname)

    fname = fdir + 'fig_gap_dynamics.pdf'
    plt.savefig(fname, dpi=300)



