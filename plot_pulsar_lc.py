import numpy as np
import sys, os

import matplotlib

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import get_cmap
from matplotlib import colorbar

from numpy import pi, log, log10, exp, sin, cos, tan

# runko tools
import pytools
from init_problem import Configuration_Turbulence as Configuration
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
    conf = Configuration(args.conf_filename, do_print=False)
    fdir = conf.outdir + "/"

    print("reading {}".format(fdir))

    toolset = QEDToolset(conf)


    #--------------------------------------------------
    #fig = plt.figure(1, figsize=(3.25, 4.0)) # single figure
    fig = plt.figure(1, figsize=(7.0,  4.5)) # two-column figure

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
        nrow_fig = 2
        ncol_fig = 3

        gs = plt.GridSpec(nrow_fig, ncol_fig)

        gs.update(wspace = 0.3)
        gs.update(hspace = 0.3)

        axs = np.empty( (nrow_fig,ncol_fig), dtype=object)

        for j in range(ncol_fig):
            for i in range(nrow_fig):
                axs[i,j] = plt.subplot(gs[i,j])
                axs[i,j].minorticks_on()

    axs[0,0].set_title(r"center at $\mu$")
    axs[0,1].set_title(r"center at $\Omega$")

    #--------------------------------------------------
    # read slice file

    #--------------------------------------------------

    norm = {}
    norm['rho'] = conf.nGJ
    norm['j']   = abs(conf.qe)*conf.nGJ*conf.cfl**2
    norm['e']   = abs(conf.qe)*conf.nGJ*conf.rad_pcap

    xmin = -conf.Lx//2
    xmax = +conf.Lx//2

    ymin = -conf.Ly//2
    ymax = +conf.Ly//2

    xmin /= conf.rad_pcap
    ymin /= conf.rad_pcap
    xmax /= conf.rad_pcap
    ymax /= conf.rad_pcap

    xobs_tj = []
    yobs_tj = []
    axs[0,2].set_xlim(( 0.0, 1.0))
    axs[0,2].set_ylim((-0.1, 0.1))


    incls = np.linspace(0, 180, 180)

    nt = int(conf.period_star/conf.interval)
    vmap = np.zeros((nt, len(incls)))


    it = -1
    #for lap in [args.lap]:
    for lap in range(0, conf.Nt, conf.interval): # loop over time snapshots
        it += 1 # time step for light curve
        #it = int(lap)
        
        fname = fdir + 'slices-xy_' + str(lap) + '.h5'

        if not(os.path.isfile(fname)):
            continue

        print('lap:', lap, fname)

        f5F = h5.File(fname,'r')

        rho = pytools.read_h5_array(f5F, 'rho', stride=1 ) 

        jx = pytools.read_h5_array(f5F, 'jx', stride=1 ) 
        jy = pytools.read_h5_array(f5F, 'jy', stride=1 ) 
        jz = pytools.read_h5_array(f5F, 'jz', stride=1 ) 

        ex = pytools.read_h5_array(f5F, 'ex', stride=1 ) 
        ey = pytools.read_h5_array(f5F, 'ey', stride=1 ) 
        ez = pytools.read_h5_array(f5F, 'ez', stride=1 ) 

        bx = pytools.read_h5_array(f5F, 'bx', stride=1 ) 
        by = pytools.read_h5_array(f5F, 'by', stride=1 ) 
        bz = pytools.read_h5_array(f5F, 'bz', stride=1 ) 

        #--------------------------------------------------
        # image extent
        nx, ny, nz = np.shape(rho)

        #val = rho[:,:,0]/norm['rho']
        #val =  jz[:,:,0]/norm['j']
        val =  ez[:,:,0]/norm['e']

        #if True:
        if lap == args.lap:
            im = pytools.visualize.imshow(
               axs[0,0], 
               val,
               xmin, xmax, ymin, ymax,
               cmap = 'RdBu', #args['cmap'],
               vmin = -1,        #args['vmin'],
               vmax =  1,        #args['vmax'],
               clip = None,      #args['clip'],
               aspect=1,         #args['aspect'],
               #norm = matplotlib.colors.LogNorm(vmin=1e-4, vmax=10.0),
               norm = matplotlib.colors.SymLogNorm(1e-2, vmin=-10.0, vmax=10.0, base=10),
               )


        #--------------------------------------------------
        # pulsar's observer trajectory

        cenx = conf.Lx//2
        ceny = conf.Lx//2
        cenz = -conf.rad_star + conf.rad_curv_shift

        rad = conf.rad_star # R_*
        h = int(0.8*conf.Lz) - cenz # height of the slice

        chi = conf.chi # obliquity angle
        phase = lap*conf.Om_star # phase of the pulsar

        xom = sin(chi)*cos(phase) # x location of the \Omega vector
        yom = sin(chi)*sin(phase) # y loc

        xom *= h/rad
        yom *= h/rad

        #print('om x y', xom, yom)

        axs[0,0].plot( [0.0], [0.0], 'kx', alpha=0.3 ) # magnetic axis
        axs[0,0].plot( [xom], [yom], 'k.', alpha=0.5 ) # rotation vector axis

        #--------------------------------------------------

        xmin2 = xmin - xom
        xmax2 = xmax - xom
        ymin2 = ymin - yom
        ymax2 = ymax - yom

        if lap == args.lap:
        #if True:
            im = pytools.visualize.imshow(
               axs[0,1], 
               val,
               xmin2, xmax2, ymin2, ymax2,
               cmap = 'RdBu', #args['cmap'],
               vmin = -1,        #args['vmin'],
               vmax =  1,        #args['vmax'],
               clip = None,      #args['clip'],
               aspect=1,         #args['aspect'],
               #norm = matplotlib.colors.LogNorm(vmin=1e-4, vmax=10.0),
               norm = matplotlib.colors.SymLogNorm(1e-2, vmin=-10.0, vmax=10.0, base=10),
               )


            axs[0,1].set_xlim((-3,3))
            axs[0,1].set_ylim((-3,3))


        axs[0,1].plot( [0.0-xom], [0.0-yom], 'kx', alpha=0.3, zorder=1 ) # magnetic axis
        axs[0,1].plot( [xom-xom], [yom-yom], 'k.', alpha=0.5, zorder=1 ) # rotation vector axis

        #--------------------------------------------------
        # observer
        phase_obs = 0.0 # observer phase
        #incl = np.deg2rad(120.0) # observer inclination (colatitude from the \Omega vector)

        #ovec = np.array([0,0,1]) # rotation axis \Omega
        #kvec = np.array([sin(incl), 0, cos(incl)]) # vector to the direction of the observer
        #nvec = np.array([ sin(chi)*cos(phase_obs), sin(chi)*sin(phase_obs), cos(chi)]) # normal vector of the slice plane (spot)

        # observer location to data values
        def obs2data(x, y, conf):
            i = x*conf.rad_pcap + conf.Lx//2
            j = y*conf.rad_pcap + conf.Ly//2

            return int(i), int(j)


        icol = 0
        for iz, incl in enumerate(np.deg2rad(incls)):

            xobs = 0.0
            yobs = -h*tan( incl - chi )/rad

            #xobs = -h*sin(incl)/rad
            #yobs = 0.0

            ## observer is defined w.r.t. \Omega vector
            xobs += xom
            yobs += yom

            #xobs_tj.append(xobs)
            #yobs_tj.append(yobs)
            #axs[0,0].plot( xobs_tj, yobs_tj, color='C0', lw=0.3, alpha=0.5 ) # observer axis

            i,j = obs2data(xobs,yobs, conf)
            #print('ij', i,j)

            vobs = np.mean( val[i-2:i+2, j-2:j+2] )

            vmap[it, iz] = vobs 

            if iz % 20 == 0:
                axs[0,0].plot( [xobs], [yobs],               color='C'+str(icol), marker='.', markersize=0.9, alpha=0.8 ) # observer axis
                axs[0,1].plot( [xobs-xom], [yobs-yom],       color='C'+str(icol), marker='.', markersize=0.9, alpha=0.8 ) # observer axis
                axs[0,2].plot( [phase/(2*pi)], vmap[it, iz], color='C'+str(icol), marker='.', markersize=0.5, alpha=1.0)
                icol += 1

        #-------------------------------------------------- 
        # plot
        #for iz in range(0, len(incls), 10):
        #    axs[0,2].plot( [phase/(2*pi)], vmap[it, iz], 'k.')



    # skymap
    im = pytools.visualize.imshow(
           axs[1,2], 
           vmap,
           0.0,   1.0,   0.0,   180.0,
           #xmin2, xmax2, ymin2, ymax2,
           cmap = 'RdBu', #args['cmap'],
           vmin = -0.1,        #args['vmin'],
           vmax =  0.1,        #args['vmax'],
           clip = None,      #args['clip'],
           #aspect=1,         #args['aspect'],
           #norm = matplotlib.colors.LogNorm(vmin=1e-4, vmax=10.0),
           #norm = matplotlib.colors.SymLogNorm(1e-2, vmin=-10.0, vmax=10.0, base=10),
           )

    axs[1,2].set_title(r"skymap")
    axs[1,2].set_xlim((0.0, 1.0))
    axs[1,2].set_ylim((0.0, 180.0))

    axs[1,2].set_xlabel(r"phase $\phi/2\pi$")
    axs[1,2].set_ylabel(r"inclination $i$")


    #--------------------------------------------------
    axleft    = 0.14
    axbottom  = 0.08
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
    
    print('\n\n')
    print('ax pos:', pos)
    fig.subplots_adjust(left=axleft, bottom=axbottom, right=axright, top=axtop)


    #slap = str(args.lap).rjust(5, '0')
    slap = str(lap).rjust(5, '0')

    #fname = fdir + 'fig_casc_' + slap + '.pdf' 
    #plt.savefig(fname)

    fname = fdir + 'pulsar_lc_' + slap + '.png' 
    plt.savefig(fname, dpi=300)



