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

    axs[0,0].set_xlabel(r"$\chi_x$")
    axs[1,0].set_xlabel(r"$\chi_x$")

    axs[0,0].set_ylabel(r"$T(\chi_x)$")
    axs[1,0].set_ylabel(r"$\Delta(x)$")
    axs[2,0].set_ylabel(r"$T(\chi_x) - \mathrm{normfun}$")

    #axs[0,0].set_xlabel(r"$x~(\mathrm{cm})$")
    #axs[0,0].set_ylabel(r"$y_\mathrm{c}~(\mathrm{cm}\,\mathrm{s}^{-1})$")

    for j in range(ncol_fig):
        for i in range(nrow_fig):
            axs[i,j].set_xscale('log')
            axs[i,j].set_xlim((1e-4, 1e3))

    axs[0,0].set_ylim((1e-15, 1e4))
    axs[0,0].set_yscale('log')
    #axs[0,0].set_ylim((-17, -6))
    axs[1,0].set_ylim((-0.01, 0.01))
    #axs[1,0].set_ylim((-0.1, 0.1))
    #axs[1,0].set_ylim((-1.0, 1.0))

    axs[2,0].set_ylim((-3.0, 1))

    #--------------------------------------------------

    if True:
        path = "multiphoton_Breit_Wheeler_tables.h5"
        f = h5py.File(path, "r")


        if 'integration_dt_dchi' in f:
            dataset = f['integration_dt_dchi']
            
            size_photon_chi = dataset.attrs["size_photon_chi"]
            chiph_min = dataset.attrs["min_photon_chi"]
            chiph_max = dataset.attrs["max_photon_chi"]

            print("Table T")
            print("Dimension: {}".format(size_photon_chi))
            print("Min particle chi: {}".format(chiph_min))
            print("Max particle chi: {}".format(chiph_max))

            chi = np.logspace(np.log10(chiph_min),np.log10(chiph_max),size_photon_chi)

            d = dataset[()]

            axs[0,0].plot(chi,d, lw=1.0, color='C0')

    else:
        chi = np.logspace(-4, 3, 512)


    def norm_fun(x):
        #return np.zeros_like(x)
        #return np.log10( x/1e10 )
        #return x/1e10 

        #def bkn_pow_smooth(x, A, x_b, a_1, a_2, delta=1):
        #    a_1 *= -1
        #    a_2 *= -1
        #    return A*(x/x_b)**(-a_1) * (0.5*(1+(x/x_b)**(1/delta)))**((a_1-a_2)*delta)
        #return bkn_pow_smooth(x, 3e-10, 1.0, 2, 1, delta=1.45)

        #return 1.9*exp(-2.8/x)*x**(5/3) # ??? asymptotic
        #return 1.68*exp(-2.8/x)*x**(1.7) # better fit

        #return 1.68*exp(-2.8/x)*x**(1.7)*(1.0 - 0.275*exp(-(log(x)-0.898)**2/7.5677))
        
        base = 1.68*exp(-2.8/x)*x**(1.7)

        a1 = 0.275
        mu1 = -0.89846
        sig1 = 7.6766
        return base*(1.0 - a1*exp(-(log(x)-mu1)**2/sig1))



    def plot_error(y, **fmt):
        #axs[1,0].plot(chi, y/d - 1, **fmt)
        axs[1,0].plot(chi, (y-d)/d, **fmt)
        #axs[1,0].plot(chi, y-d, **fmt)

    axs[2,0].plot(chi, d - norm_fun(chi), lw=1, color='C2', linestyle='dashed')

    axs[0,0].plot(chi, norm_fun(chi), lw=1, color='C1', linestyle='dotted')
    plot_error(norm_fun(chi), color='C1', linestyle='dotted', lw=1)


    #--------------------------------------------------
    # fit
    if True:
        from scipy.optimize import minimize

        def model_asy_lngauss(theta, x):
            a1, mu1,sig1 = theta

            base = 1.68*exp(-2.8/x)*x**(1.7)
            return base*(1.0 - a1*exp(-(log(x)-mu1)**2/sig1))

        def log_prior(theta):
            a1, mu1,sig1 = theta

            if a1 < 0.0 or a1 > 1.0: return -np.inf
            if mu1 < -5.0 or a1 > 3.0: return -np.inf
            if sig1 < 0.0 or sig1 > 20.0: return -np.inf

            return 0


        # Gaussian fitting function
        def fit_func_gaussian(theta, x, y, fun):
        
            lp = log_prior(theta) #prior
            if not np.isfinite(lp): return -np.inf

            x = x[180:300]
            y = y[180:300]
        
            ymodel = fun(theta, x)

            sigma2 = np.ones(len(x))/x**1
            return -0.5*np.sum((y-ymodel)**2/sigma2 + np.log(sigma2)) + lp
        

        fit_params = {
                'model':   model_asy_lngauss,
                'names':   ['a', 'mu', 'sig',],
                'initial': [0.5, 0.0,  3.0],
                }

        # now fit

        axs[1,0].axvline(chi[180])
        axs[1,0].axvline(chi[300])

        fun = fit_params['model']
        nll  = lambda *args: -fit_func_gaussian(*args)
        soln = minimize( 
                        nll, 
                        np.array(fit_params['initial']), 
                        args=(chi, d, fun), 
                        method='Nelder-Mead', 
                        tol=1.0e-5,
                        )

        print("Maximum likelihood estimates: {0:.5f}".format(np.log10(soln.fun)))
        print("model parameters: ", soln.x)
        for i in range(len(soln.x)):
            print("{} \t : {}".format(fit_params['names'][i], soln.x[i]))

        axs[0,0].plot(chi, model_asy_lngauss(soln.x, chi), lw=1.5, color='C5', linestyle='dashed')

        plot_error( model_asy_lngauss(soln.x, chi),
                   lw=1.5, color='C5', linestyle='dashed')

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
        #return square(-1.9522) / (log(x/1.5746) - ((x + x) - (-0.11058 / x))) + norm_fun(x)
        #return exp(1.2 / ((-2.4412 - x) + (-0.38945 / x))) * norm_fun(x)
        #return exp((1.1973 - (-0.0047227 / x)) / ((-2.3762 - x) + (-0.43761 / x)))*norm_fun(x)
        return x

    #d_sr = t1(chi)
    #axs[0,0].plot(chi, d_sr, color='C4', linestyle='dotted')
    #nlot_error(d_sr, color='C4', linestyle='solid', lw=1)


    #--------------------------------------------------
    # compare to code 

    if True:
        import pyrunko

        mpp = pyrunko.qed.MultiPhotAnn("ph")
        mpp.B_QED = 1.0

        N = 512
        d_runko = np.zeros(N)
        chi_runko = np.zeros(N)

        vals = np.logspace(-4, 3, N)

        for i in range(N):
            ux,uy,uz = 0.3, 0.3, 0.3
            ex,ey,ez = 0.0, 0.0, 0.0
            bx,by,bz = 0.0, 0.0, 1.0

            RAMP = vals[i]  # vary B field (to get evolving chi)

            tau = mpp.comp_optical_depth(
                    "e-",
                    ux,uy,uz,
                    ex,ey,ez,
                    bx,by,RAMP*bz)

            c = mpp.comp_chi( ux,uy,uz, ex,ey,ez, bx,by,RAMP*bz)

            d_runko[i] = tau
            chi_runko[i] = c
        
        axs[0,0].plot(chi_runko, d_runko, color='C6', linestyle='dotted')
        plot_error(d_runko, color='C6', linestyle='solid', lw=1)


    #--------------------------------------------------
    axleft    = 0.20
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

    fname = 'integ_T.pdf' 
    plt.savefig(fname)

    #fname = 'pml.png' 
    #plt.savefig(fname, dpi=300)
















