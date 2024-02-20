import numpy as np
import sys, os
import h5py as h5py

import matplotlib

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm
from matplotlib.cm import get_cmap
from matplotlib import colorbar

from numpy import pi, log, log10, exp, sin, cos, tanh


#-------------------------------------------------- 
if __name__ == "__main__":

    #--------------------------------------------------
    #fig = plt.figure(1, figsize=(3.25, 6.0)) # single figure
    #fig = plt.figure(1, figsize=(7.0,  6.0)) # two-column figure
    fig = plt.figure(1, figsize=(9.0,  3.0)) # two-column figure

    plt.rc('xtick', top   = True)
    plt.rc('ytick', right = True)

    plt.rc('font',  family='serif',)
    plt.rc('text',  usetex=False)

    plt.rc('xtick', labelsize=7)
    plt.rc('ytick', labelsize=7)
    plt.rc('axes',  labelsize=8)
    plt.rc('legend',  handlelength=4.0)

    nrow_fig = 1
    ncol_fig = 3

    gs = plt.GridSpec(nrow_fig, ncol_fig)

    gs.update(wspace = 0.25)
    gs.update(hspace = 0.30)

    axs = np.empty( (nrow_fig,ncol_fig), dtype=object)

    for j in range(ncol_fig):
        for i in range(nrow_fig):
            axs[i,j] = plt.subplot(gs[i,j])
            axs[i,j].minorticks_on()

    axs[0,0].set_xlabel(r"$\chi_x$")
    #axs[1,0].set_xlabel(r"$\chi_\pm$")

    axs[0,0].set_ylabel(r"$\chi_\pm$")
    #axs[1,0].set_ylabel(r"$\chi_x$")

    axs[0,1].set_xlabel(r"$\chi_\pm$")
    #axs[1,1].set_xlabel(r"$\chi_x$")

    axs[0,1].set_ylabel(r"$F |_{\chi_x = \mathrm{const}}$")
    #axs[1,1].set_ylabel(r"$F |_{\chi_\pm = \mathrm{const}}$")

    axs[0,2].set_xlabel(r"$\chi_x$")
    #axs[1,2].set_xlabel(r"$\chi_\pm$")
    axs[0,2].set_ylabel(r"$F |_{\chi_\pm = \mathrm{const}}$")
    #axs[1,2].set_ylabel(r"$F |_{\chi_x = \mathrm{const}}$")


    #axs[0,0].set_xlabel(r"$x~(\mathrm{cm})$")
    #axs[0,0].set_ylabel(r"$y_\mathrm{c}~(\mathrm{cm}\,\mathrm{s}^{-1})$")

    for j in range(ncol_fig):
        for i in range(nrow_fig):
            axs[i,j].set_xscale('log')
            axs[i,j].set_yscale('log')
            axs[i,j].set_xlim((1e-2, 1e3))
            axs[i,j].set_ylim((1e-2, 1e3))

    #axs[0,0].set_ylim((-2, 6))
    #axs[1,0].set_ylim((-0.03, 0.03))
    #axs[1,0].set_ylim((-1.0, 1.0))

    axs[0,1].set_ylim((1e-3, 1.1))
    axs[0,2].set_ylim((1e-3, 1.1))
    #axs[1,1].set_ylim((1e-3, 1.1))
    #axs[2,1].set_ylim((1e-3, 1.1))

    #--------------------------------------------------

    # approximate function behavior
    def norm_fun(x0, x1):
        #return 10**(-3*3/7)*x1**(3/7) # at x = 1e3
        return 10**(-3/5)*x1**(3/5) # at x = 1e-4


    if True:
        path = "multiphoton_Breit_Wheeler_tables_64.h5"
        f = h5py.File(path, "r")

        if 'xi' in f:
        
            # TODO why does this only go up to 1/2 ?
            xi = 2.0*np.array(f['xi']).T

            size_photon_chi = f['xi'].attrs["size_photon_chi"]
            chiph_min = f['xi'].attrs["min_photon_chi"]
            chiph_max = f['xi'].attrs["max_photon_chi"]
            chipa_dim = f['xi'].attrs["size_particle_chi"]

            print("")
            print("Table xi")
            print("Dimension: {}".format(size_photon_chi))
            print("Min particle chi: {}".format(chiph_min))
            print("Max particle chi: {}".format(chiph_max))

            #chipa = np.linspace(1,chipa_dim+1,chipa_dim+1)
            #chipa = np.linspace(1,chipa_dim+1,chipa_dim)
            chipa = np.logspace(np.log10(chiph_min),np.log10(chiph_max),size_photon_chi)
            chiph = np.logspace(np.log10(chiph_min),np.log10(chiph_max),size_photon_chi)


            #im = axs[0,0].pcolormesh(chipa,chiph,xi[:,1:],
            #im = axs[0,0].pcolormesh(chipa,chiph,xi[:,:],
            #                cmap=get_cmap('plasma'),
            #                #shading='none',
            #                linewidth=0,
            #                rasterized=True,
            #                    )

            #im.set_norm(LogNorm())
            #cb = colorbar(im,ax=ax0)

    print('chi pa', chipa)
    print('chi ph', chiph)

    print('---------------------------------------------')
    print(' reduced xi table for storing ')


    N = 64
    sk = 2

    print('range')
    print(list(range(0,N,sk)))

    print('chi pa')
    for i in range(0, N, sk):
        #print(i, chipa[i])
        print('{:1.7e}'.format(chiph[i]))

    print()
    print('xi')
    print()

    for j in range(0, N, sk):
        print('{ ', end='')
        for i in range(0, N, sk):
            #print(xi[i,j], end=',')
            print('{:1.7e}'.format(xi[i,j]), end=',')
        print('1.0000000e+00},')

    print('---------------------------------------------')

    #chipa = chipa[::8]
    #chiph = chiph[::8]
    #xi = xi[::8,::8]
    #size_particle_chi = 32 
    #size_photon_chi = 32 


    #chiph = np.logspace(np.log10(chipa_min),np.log10(chipa_max),size_particle_chi)
    #xi_app = np.zeros_like(xi)
    #for i in range(len(chipa)):
    #    for j in range(len(chiph)):
    #        xi_app[i,j] = fun_xi_approx(chipa[i], chiph[j])

    #im2 = axs[1,0].pcolormesh(chipa,chiph,xi_app[:,:],
    #                        cmap=get_cmap('plasma'),
    #                        shading='none',
    #                        linewidth=0,
    #                        rasterized=True,
    #                            )
    #im2.set_norm(LogNorm())


    def plot_xi(row, xi):

        #im = axs[0,0].pcolormesh(chipa,chiph,xi[:,1:],
        #im = axs[0,0].pcolormesh(chipa,chiph,xi[:,:],
        #                    cmap=get_cmap('plasma'),
        #                    #shading='none',
        #                    linewidth=0,
        #                    rasterized=True,
        #                        )
        #im.set_norm(LogNorm())
        #cb = colorbar(im,ax=ax0)

        axs[row,0].contourf(chiph, chipa, xi[:,:],
                          cmap=get_cmap('plasma'),
                          levels=[1e-4,3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 0.3, 0.6, 0.8, 0.9, 1],
                          norm=LogNorm(),
                          vmin=1e-4,
                          vmax=1.0,
                          )

        #for i in range(0, size_particle_chi, 40):
        for i in [0, 40, 80, 120, 255]:
            i = int( (i/256)*len(chipa) )
            print("chipa", i, chipa[i])
            axs[row,1].plot(chipa, xi[:, i])

        #for i in [0, 40, 80, 120, 255]:
        for i in range(0,256,8):
            i = int( (i/256)*len(chipa) )
            print("chiph", i, chiph[i])
            axs[row,2].plot(chiph, xi[i, :])

    plot_xi(0, xi)



    def plot_error(row, d, **fmt):

        im = axs[row,0].pcolormesh(
                chiph,chipa,
                (xi[:,:] - d[:,:])/xi[:,:],
                cmap=get_cmap('RdBu'),
                #shading='none',
                linewidth=0,
                rasterized=True,
                vmin=-0.1,
                vmax=+0.1,
                )

        icol = 0
        for i in [0, 40, 80, 120, 255]:
            i = int( (i/256)*len(chipa) )

            #print("chipa", i, chipa[i])
            axs[row,1].plot(chipa, xi[:, i], color='C'+str(icol))
            axs[row,1].plot(chipa,  d[:, i], color='C'+str(icol), linestyle='dotted')
            icol += 1

        icol = 0
        #for i in [0, 40, 80, 120, 255]:
        for i in range(0,256,8):
            i = int( (i/256)*len(chipa) )

            #print("chiph", i, chiph[i])
            axs[row,2].plot(chiph, xi[i, :], color='C'+str(icol))
            axs[row,2].plot(chiph,  d[i, :], color='C'+str(icol), linestyle='dotted')
            icol += 1

        loss = np.sum( (xi[:,:] - d[:,:])**2 )
        rloss = np.sum( abs(xi[:,:] - d[:,:])/xi[:,:] )
        print('LOSS:', loss, rloss)

    #axs[0,0].plot(chi,d - norm_fun(chi), lw=1, color='C2', linestyle='dashed')
    #axs[0,0].plot(chi,norm_fun(chi), lw=1, color='C3', linestyle='dotted')
    #plot_error(norm_fun(chi), color='C3', linestyle='dotted', lw=1)



    #--------------------------------------------------
    #--------------------------------------------------
    #--------------------------------------------------
    if False:
        from pysr import PySRRegressor
        
        model = PySRRegressor(
            procs=4,
            populations=12,
            niterations=500,  # < Increase me for better results
            ncyclesperiteration=500,
            binary_operators=["+", "-", "*", "/"],
            unary_operators=[
                #"pow",
                "log",
                "exp",
                "square",
                "cube",
                #"sq3(x) = x^(-1/3)",
            ],
            #extra_sympy_mappings={"sq3": lambda x: x**(-1/3) },
            # ^ Define operator for SymPy as well
            loss="loss(prediction, target) = (prediction - target)^2",
            #loss="loss(prediction, target) = abs(prediction - target)/target",
            # ^ Custom loss function (julia syntax)
            #niterations=10000,  
            #early_stop_condition=(
            #    "stop_if(loss, complexity) = loss < 1e-6 && complexity < 10"
            #    # Stop early if we find a good and simple equation
            #),
            timeout_in_seconds=60 * 60 * 1,
            maxsize=35,
            maxdepth=5,
            constraints={
                #'pow': (9, 1),
                'exp': -1, #(-1, -1),
                'log': -1, #(-1, -1),
                'square': -1, #(-1, -1),
                'cube': -1, #(-1, -1),
                #'sq3': (9, 1),
                },
            # ^ Limit the complexity within each argument.
            # "inv": (-1, 9) states that the numerator has no constraint,
            # but the denominator has a max complexity of 9.
            # "exp": 9 simply states that `exp` can only have
            # an expression of complexity 9 as input.
            nested_constraints={
                "square": {"square": 0, "cube": 0, "exp": 0},
                "cube": {"square": 0, "cube": 0, "exp": 0},
                "exp": {"square": 1, "cube": 1, "exp": 0},
                "log": {"square": 1, "cube": 1, "exp": 0},
            },
            # ^ Nesting constraints on operators. For example,
            # "square(exp(x))" is not allowed, since "square": {"exp": 0}.
            #complexity_of_operators={"/": 2, "exp": 3},
            precision=64,
        )

        #X = chi.reshape(-1,1) # reshape to have 1 feature and N samples

        X = np.zeros(( (size_photon_chi//4)**2, 2))
        y = np.zeros(( (size_photon_chi//4)**2))

        for i in range(0, size_photon_chi, 4):
            for j in range(0, size_photon_chi, 4):
                X[i+j,0] = chipa[i]
                X[i+j,1] = chiph[j]
                y[i+j] = xi[j,i]

        model.fit(X, y)
        #print(model)


    def estimate_sr(x0,x1):
        def square(x): return x**2
        def cube(x): return x**3

        def tr(x):
            return 1e-1*x**(-1/4)*exp(-x)

        def tr_loglin(xx, a, b, xx0, xx1):
            x0 = log10(xx0)
            x1 = log10(xx1)
            x  = log10(xx)

            if x > x1: x = x1
            if x < x0: x = x0

            k = (b - a)/(x1 - x0) # slope
            c = (a*x1 - b*x0)/(x1 - x0)
            return k*x + c
        

        kA = tr_loglin(x0, -0.59, -1.29, 1e-4, 3e0)
        ks = tr_loglin(x0,  3/5,    3/7, 1e-4, 3e0)

        kd = tr_loglin(x0,     12,   500, 1e-4, 3e0)

        t = 10**kA
        t *= x1**ks 
        t *= exp(-(x1/kd))
        t += (1-exp(-(x1/(1.4*kd))))

        return np.clip(0.0, 1.0, t)

        # ver0
        #tlow   = 10**(-3/5  )*x1**(3/5)*exp(-x1/10) + (1-exp(-x1/13))    # at x = 1e-4
        #thigh  = 10**(-3*3/7)*x1**(3/7) # at x = 1e3
        #t = tlow*tr(x0) + (1-tr(x0))*thigh
        #return t

        #return 10**(-3*3/7)*x1**(3/7) # at x = 1e3
        #return 10**(-3/5)*x1**(3/5)*exp(-x1/10) + (1-exp(-x1/13))    # at x = 1e-4


    xi_sr = np.zeros_like(xi)
    for i in range(0, size_photon_chi):
        for j in range(0, size_photon_chi):
            xi_sr[j, i] = estimate_sr( chiph[i], chipa[j] )

    #d_sr = t1(chi)
    #axs[0,0].plot(chi, d_sr, color='C4', linestyle='dotted')
    #plot_error(d_sr, color='C4', linestyle='solid', lw=1)

    #plot_xi(1, xi_sr)
    #plot_error(2, xi_sr)

    #--------------------------------------------------
    axleft    = 0.10
    axbottom  = 0.15
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

    fname = 'integ_T_xi.pdf' 
    plt.savefig(fname)

    #fname = 'pml.png' 
    #plt.savefig(fname, dpi=300)
















