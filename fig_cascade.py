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
    fig = plt.figure(1, figsize=(3.25, 4.0)) # single figure
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
    if False: # regular gridspec
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

    #--------------------------------------------------
    else:

        nrow_fig = 7
        ncol_fig = 2
        axs = np.empty( (nrow_fig,ncol_fig), dtype=object)

        gs = plt.GridSpec(100, 100)
        gs.update(wspace = 0.0) # in between space
        gs.update(hspace = 0.0)

        cw = 22
        cr = 40
        axs[0,0] = plt.subplot(gs[0:cw, 0:cr]) # row, col
        axs[0,1] = plt.subplot(gs[0:cw, -cr:]) # row, col

        axs[0,0].minorticks_on()
        axs[0,1].minorticks_on()

        axs[0,0].tick_params(labelbottom=False)
        axs[0,1].tick_params(labelbottom=False)

        axs[0,0].tick_params(labeltop=True)
        axs[0,1].tick_params(labeltop=True)

        axs[0,0].set_yscale('log')
        axs[0,1].set_yscale('log')

        axs[0,0].set_ylim((1e-8, 1e-3))
        axs[0,1].set_ylim((1e-8, 1e-3))

        axs[0,0].set_ylabel(r"$p\, \mathrm{d} \tau/\mathrm{d}p$")
        axs[0,1].set_ylabel(r"$p\, \mathrm{d} \tau/\mathrm{d}p$")

        ph = 24 # panel height
        pad = 0 # padding between panels

        row1 = cw + 4
        row2 = row1 + ph 
        for i in range(1,4):
            print('row', row1, row2)
            axs[i,0] = plt.subplot(gs[row1:row2, 0:100]) # row, col

            row1 += pad + ph
            row2 += pad + ph

            axs[i,0].minorticks_on()

        #print(axs[2,0].spines['top'])
        #axs[1,0].set_xticklabels('')
        #axs[2,0].set_xticklabels('')

        axs[2,0].tick_params(top=False, which='both')
        axs[3,0].tick_params(top=False, which='both')

        axs[1,0].tick_params(labelbottom=False)
        axs[2,0].tick_params(labelbottom=False)

        axs[3,0].set_xlabel(r"$h/R_{\mathrm{pc}}$")

        axs[1,0].set_ylabel(r"$p_- ~ (m_e c)$")
        axs[2,0].set_ylabel(r"$p_+ ~ (m_e c)$")
        axs[3,0].set_ylabel(r"$x ~ (m_e c^2)$")



    if True:
        hmin = -0.4
        hmax =  2.4

        #for j in range(ncol_fig):
        #for i in range(nrow_fig):
        for i in [1,2,3]:
            axs[i,0].set_xlim((hmin, hmax))

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
        #hhlims[0] = 0 #(0 - ceny)/conf.rad_pcap
        hhlims[0] = -(conf.rad_curv_shift + conf.height_atms)/conf.rad_pcap #(0 - ceny)/conf.rad_pcap
        hhlims[1] = Lh/conf.rad_pcap

        print('hhlims', hhlims)

        print(toolset.pxlims[1],  toolset.pxlims[0] )
        px_log_extent = toolset.pxlims[1] - toolset.pxlims[0] 
        print('px_log_extent:', px_log_extent)


        for i in range(1,4): 
            axs[i,0].set_ylim((-px_log_extent, px_log_extent))
        axs[0,0].set_xlim((-px_log_extent, px_log_extent))
        axs[0,1].set_xlim((-px_log_extent, px_log_extent))


        print(toolset.xxlims[1],  toolset.xxlims[0] )
        xx_log_extent = toolset.xxlims[1] - toolset.xxlims[0] 
        print('xx_log_extent:', xx_log_extent)

        imem = axs[1,0].imshow( np.zeros((toolset.Nhist, 2*toolset.Nhist)),
                         extent=[hhlims[0], hhlims[1], 
                                 -px_log_extent, px_log_extent], 
                               #  -toolset.Nhist, toolset.Nhist], # TODO
                               #toolset.pxlims[0], toolset.pxlims[1]],
                         origin='lower',
                         cmap='turbo',
                         aspect='auto',
                         interpolation='nearest',
                         vmin=np.log10(toolset.pylims[0]),
                         vmax=np.log10(toolset.pylims[1]),
                         zorder=1,
                         )

        imep = axs[2,0].imshow( np.zeros((toolset.Nhist, 2*toolset.Nhist)),
                         extent=[hhlims[0], hhlims[1], 
                                 -px_log_extent, px_log_extent], 
                                 #-toolset.Nhist, toolset.Nhist], # TODO
                                 #toolset.pxlims[0], toolset.pxlims[1]],
                         origin='lower',
                         cmap='turbo',
                         aspect='auto',
                         interpolation='nearest',
                         vmin=np.log10(toolset.pylims[0]),
                         vmax=np.log10(toolset.pylims[1]),
                         zorder=1,
                         )

        # height vs ene; photons
        imph = axs[3,0].imshow( np.zeros((toolset.Nhist, 2*toolset.Nhist)),
                         extent=[hhlims[0], hhlims[1], 
                                 -xx_log_extent, xx_log_extent], 
                               #-toolset.Nhist, toolset.Nhist], # TODO
                               #toolset.xxlims[0], toolset.xxlims[1]],
                          origin='lower',
                          cmap='turbo',
                          aspect='auto',
                          interpolation='nearest',
                          vmin=np.log10( toolset.xylims[0]),
                          vmax=np.log10( toolset.xylims[1]),
                          zorder=1,
                          )


        xs1 = int(0.25*toolset.Nhist)
        xs2 = int(0.48*toolset.Nhist)

        #hem[xs1, :] = 10
        #hem[xs2, :] = 10

        imem.set_data(np.log10( hem.T)) 
        imep.set_data(np.log10( hep.T)) 
        imph.set_data(np.log10( hph.T)) 

        #--------------------------------------------------
        # get slices

        w = 3
        h1_em_1 = np.sum(hem[xs1-w:xs1+w,:], axis=0)
        h1_ep_1 = np.sum(hep[xs1-w:xs1+w,:], axis=0)
        h1_em_2 = np.sum(hem[xs2-w:xs2+w,:], axis=0)
        h1_ep_2 = np.sum(hep[xs2-w:xs2+w,:], axis=0)

        xvals = np.linspace(-px_log_extent, px_log_extent, 2*toolset.Nhist)
        hvals = np.linspace(hhlims[0], hhlims[1], toolset.Nhist)

        axs[1,0].axvline(hvals[xs1], lw=0.7, color='C0', linestyle='solid', zorder=0)
        axs[2,0].axvline(hvals[xs1], lw=0.7, color='C1', linestyle='solid', zorder=0)

        axs[1,0].axvline(hvals[xs2], lw=0.7, color='C0', linestyle='dashed', zorder=0)
        axs[2,0].axvline(hvals[xs2], lw=0.7, color='C1', linestyle='dashed', zorder=0)

        axs[0,0].plot(xvals, h1_em_1,
                      drawstyle='steps-pre',
                      color='C0',
                      alpha=1.0,
                      lw = 0.5,
                      linestyle='solid',
                      ) 
        axs[0,0].plot(xvals, h1_ep_1,
                      drawstyle='steps-pre',
                      color='C1',
                      alpha=1.0,
                      lw = 0.5,
                      linestyle='solid',
                      ) 

        axs[0,1].plot(xvals, h1_em_2,
                      drawstyle='steps-pre',
                      color='C0',
                      alpha=1.0,
                      lw = 0.5,
                      linestyle='solid',
                      ) 
        axs[0,1].plot(xvals, h1_ep_2,
                      drawstyle='steps-pre',
                      color='C1',
                      alpha=1.0,
                      lw = 0.5,
                      linestyle='solid',
                      ) 





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

    xticks = [-6, -2, 2, 6]
    axs[0,0].set_xticks(xticks)
    axs[0,1].set_xticks(xticks)

    xticks2 = [-6, -2, 2, 6]
    axs[1,0].set_yticks(xticks)
    axs[2,0].set_yticks(xticks)
    axs[3,0].set_yticks(xticks)


    majorformatter = CustomScalarFormatter()

    axs[0,0].xaxis.set_major_formatter(majorformatter)
    axs[0,1].xaxis.set_major_formatter(majorformatter)

    axs[1,0].yaxis.set_major_formatter(majorformatter)
    axs[2,0].yaxis.set_major_formatter(majorformatter)
    axs[3,0].yaxis.set_major_formatter(majorformatter)



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
    print('ax pos:', pos)
    fig.subplots_adjust(left=axleft, bottom=axbottom, right=axright, top=axtop)


    slap = str(args.lap).rjust(5, '0')

    #fname = fdir + 'fig_casc_' + slap + '.pdf' 
    #plt.savefig(fname)

    fname = fdir + 'fig_casc_' + slap + '.png' 
    plt.savefig(fname, dpi=300)



