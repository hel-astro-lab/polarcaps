import numpy as np
import sys, os
import h5py as h5py

import matplotlib

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import get_cmap
from matplotlib import colorbar

from numpy import pi, log, log10, exp, sin, cos 


#-------------------------------------------------- 
if __name__ == "__main__":

    #--------------------------------------------------
    fig = plt.figure(1, figsize=(3.25, 4.0)) # single figure
    #fig = plt.figure(1, figsize=(7.0,  5.5)) # two-column figure

    plt.rc('xtick', top   = True)
    plt.rc('ytick', right = True)

    plt.rc('font',  family='serif',)
    plt.rc('text',  usetex=False)

    plt.rc('xtick', labelsize=7)
    plt.rc('ytick', labelsize=7)
    plt.rc('axes',  labelsize=8)
    plt.rc('legend',  handlelength=4.0)

    nrow_fig = 2
    ncol_fig = 1

    gs = plt.GridSpec(nrow_fig, ncol_fig)

    gs.update(wspace = 0.25)
    gs.update(hspace = 0.0)

    axs = np.empty( (nrow_fig,ncol_fig), dtype=object)

    for j in range(ncol_fig):
        for i in range(nrow_fig):
            axs[i,j] = plt.subplot(gs[i,j])
            axs[i,j].minorticks_on()

    axs[0,0].set_xlabel(r"$x$")
    axs[1,0].set_xlabel(r"$x$")

    axs[0,0].set_ylabel(r"$K(x)$")
    axs[1,0].set_ylabel(r"$\Delta(x)$")

    #axs[0,0].set_xlabel(r"$x~(\mathrm{cm})$")
    #axs[0,0].set_ylabel(r"$y_\mathrm{c}~(\mathrm{cm}\,\mathrm{s}^{-1})$")

    for j in range(ncol_fig):
        for i in range(nrow_fig):
            axs[i,j].set_xscale('log')
            axs[i,j].set_xlim((1e-4, 1e3))


    axs[0,0].set_ylim((-2, 6))
    axs[1,0].set_ylim((-0.03, 0.03))
    #axs[1,0].set_ylim((-1.0, 1.0))


    #--------------------------------------------------


    def norm_fun(x):
        #return 5.236*(1-exp(-x))*x**(-1/3)
        #return 5.236*exp(-0.5*x)
        return 5.236*( exp(-x) + (1-exp(-x))*x**(-1/3) ) # asymptotics
        #return  5.236*x**(-1/3) # asymptotics
        #return np.zeros_like(x)
        #return 5.236*exp(-1.5*x) + 5.235*(1-exp(-0.5*x))*x**(-0.332)  # asymptotics


    if True:
        path = "radiation_tables_256.h5"
        f = h5py.File(path, "r")
        if 'integfochi' in f:
            dataset = f['integfochi']

            size_particle_chi = dataset.attrs["size_particle_chi"]
            chipa_min = dataset.attrs["min_particle_chi"]
            chipa_max = dataset.attrs["max_particle_chi"]

            print("")
            print("Table 'integfochi'")
            print("Size: {}".format(size_particle_chi))
            print("Min particle chi: {}".format(chipa_min))
            print("Max particle chi: {}".format(chipa_max))

            chi = np.logspace(np.log10(chipa_min),np.log10(chipa_max),size_particle_chi)

            d = dataset[()]
            axs[0,0].plot(chi,d,lw=1,label=r'$\int{F/\chi d\chi}$')

            axs[0,0].set_xscale('log')
            axs[0,0].set_xlabel(r'$\chi_\pm$')
            axs[0,0].set_ylabel(r'$\int{F/\chi d\chi}$')
    else:
        chi = np.logspace(-4, 3, 512)


    def plot_error(y, **fmt):
        axs[1,0].plot(chi, y/d - 1, **fmt)

    axs[0,0].plot(chi,d - norm_fun(chi), lw=1, color='C2', linestyle='dashed')

    axs[0,0].plot(chi,norm_fun(chi), lw=1, color='C3', linestyle='dotted')
    plot_error(norm_fun(chi), color='C3', linestyle='dotted', lw=1)


    #--------------------------------------------------
    #--------------------------------------------------
    #--------------------------------------------------
    if False: # test g2

        def g(x):
            t1 = 1 + 4.8*(1 + x)*np.log(1 + 1.7*x) + 2.44*x**2
            return t1**(-2/3)

        d2 = g(chi)
        axs[0,0].plot(chi, d2, color='C1', lw=1, linestyle='dashed')


        def g2(x):
            def square(x):
                return x**2
            return log(1.5935 + (x / 1.3803)) / ((log(1.5935) + (x + x)) + (x + square(x)))
        d3 = g2(chi)
        axs[0,0].plot(chi, d3, color='C2', lw=1, linestyle='dashed')




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
            extra_sympy_mappings={"sq3": lambda x: x**(-1/3) },
            # ^ Define operator for SymPy as well
            loss="loss(prediction, target) = (prediction - target)^2",
            # ^ Custom loss function (julia syntax)
            #niterations=10000,  
            #early_stop_condition=(
            #    "stop_if(loss, complexity) = loss < 1e-6 && complexity < 10"
            #    # Stop early if we find a good and simple equation
            #),
            timeout_in_seconds=60 * 60 * 1,
            maxsize=20,
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

        X = chi.reshape(-1,1) # reshape to have 1 feature and N samples

        model.fit(X, d/norm_fun(chi) )

        #print(model)


    def t1(x):
        def square(x): return x**2
        def cube(x): return x**3

        #return ((x + x) / ((-1.5358 / 0.063052) - x)) - (cube(-1.5358) + ((x - 1.0116) / (0.63344 + x)))
        #return ((exp(0.148) - (0.14754 * x)) / (square(0.6816) + (x * 0.52895))) + exp(exp(x * -0.0155))

        #return -3.0194 / (((0.19373 + x) / x) + x) + norm_fun(x)
        #return (-3.0194 - (x/0.65976)) / ((0.81925 + (x*x)) + ((x + 0.14575) / x)) + norm_fun(x)
        #return  -3.4647 / (((1.8851 * x) + (0.0905 / x)) - log(x / 1.367)) + norm_fun(x)
        #return -(x/0.34219)/square(0.45255 + x) + norm_fun(x)
        return square(-1.9522) / (log(x/1.5746) - ((x + x) - (-0.11058 / x))) + norm_fun(x)
        #return exp(1.2 / ((-2.4412 - x) + (-0.38945 / x))) * norm_fun(x)
        #return exp((1.1973 - (-0.0047227 / x)) / ((-2.3762 - x) + (-0.43761 / x)))*norm_fun(x)


    d_sr = t1(chi)
    axs[0,0].plot(chi, d_sr, color='C4', linestyle='dotted')
    plot_error(d_sr, color='C4', linestyle='solid', lw=1)



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

    fname = 'integ_tau.pdf' 
    plt.savefig(fname)

    #fname = 'pml.png' 
    #plt.savefig(fname, dpi=300)
















