import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm, LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import h5py as h5
import sys, os

import pytools
from init_problem import Configuration_Turbulence as Configuration



def lap2time(lap, conf):
    #print("lap: {} t: {} t0: {} L: {} l_0: {}".format(lap,lap/t0, t0, L, l0))
    return lap/conf.t0


def draw_label(ax, xloc, yloc, tt, conf):
    ax.text(xloc, yloc, tt, 
            color='black', 
            fontsize=8.0,
            ha='center',
            va='center',
        bbox=dict(facecolor='lightgray', edgecolor='gray', boxstyle='round,pad=0.4', alpha=0.9, zorder=20),
        )





# trick to make nice colorbars
# see http://joseph-long.com/writing/colorbars/
#def colorbar(mappable, 
#        loc="right", 
#        orientation="vertical", 
#        size="1%", 
#        #pad=0.05, 
#        pad=0.1, 
#        ticklocation='right'):
#        #loc="top", 
#        #orientation="horizontal", 
#        #size="8%", 
#        #pad=0.03, 
#        #ticklocation='top'):
#
#    ax = mappable.axes
#    fig = ax.figure
#    divider = make_axes_locatable(ax)
#    cax = divider.append_axes(loc, size=size, pad=pad)
#    return cax, fig.colorbar(mappable, cax=cax, orientation=orientation, ticklocation=ticklocation)


default_values = {
    'cmap':"viridis",
    'vmin': None,
    'vmax': None,
    'clip': None,
    'aspect':1,
    'vsymmetric':None,
    'winsorize_min':0.005,
    'winsorize_max':0.005,
    'title':'',
    'derived':False,
    'file': 'flds',
    'log': False,
    'norm': None,
}

default_turbulence_values = {
        # DONE
        'ex': {'title': r"$E_r$",
               'cmap': "PRGn",
               'vsymmetric':True,
                'vmin':-0.1,
                'vmax': 0.1,
                'norm': SymLogNorm(1e-2, vmin=-1.0, vmax=1.0, base=10),
                },
        'exy': {'title': r"$E_{xy}$",
               'cmap': "Blues",
               'vsymmetric':True,
                'vmin':-0.1,
                'vmax': 0.1,
                'norm': LogNorm(vmin=1e-4, vmax=10.0),
               'derived':True,
                },
        'ey': {'title': r"$E_y$",
               'cmap': "PRGn",
               'vsymmetric':True,
                'vmin':-0.1,
                'vmax': 0.1,
                },
        'ez': {'title': r"$E_z$",
               'cmap': "PRGn",
               'vsymmetric':True,
                'vmin':-0.1,
                'vmax': 0.1,
                },
        'epar': {'title': r"$E_\parallel$",
               'cmap': "RdBu", #"PRGn",
               'vsymmetric':True,
               'derived':True,
                'vmin':-0.1,
                'vmax': 0.1,
                'norm': SymLogNorm(1e-2, vmin=-1.0, vmax=1.0, base=10),
                },
        'bx': {'title': r"$B_x$",
               'cmap': "RdBu",
               'vsymmetric':True,
                'vmin':-0.1,
                'vmax': 0.1,
                },
        'by': {'title': r"$B_z$",
               'cmap': "RdBu",
               'vsymmetric':True,
                'vmin':-0.1,
                'vmax': 0.1,
                },
        'bz': {'title': r"$B_z$",
               'cmap': "RdBu",
               'vsymmetric':True,
                'vmin':-1.2,
                'vmax': 1.2,
                },
        'rho': {'title': r"$n/n_\mathrm{GJ}$",
                'vmin': 000000.0,
                'vmax': 4.0,
                },
        'logrho': {'title': r"$\log_{10}(n/n_\mathrm{GJ})$",
                'vmin':-3.0,
                'vmax': 1.0,
                'derived':True,
                'log':True,
                },
        'densx':  {'title': r"$n_x/n_{x,0}$",
                'vmin': 0.0,
                'vmax': 4.0,
                'file': 'moms',
                },
        'lognx':  {'title': r"$\log_{10}(n_x/n_{x,0})$",
                'vmin':-2.0,
                'vmax': 3.0,
                'file': 'moms',
                'derived':True,
                'cmap':'inferno',
                'log':True,
                },
        'je': {'title': r"$\mathbf{J} \cdot \mathbf{E}$",
               'cmap': "BrBG",
               'vmin': -5.00000,
               'vmax':  5.00000,
               'vsymmetric':True,
               'derived':True,
                },
        'jz': {'title': r"$J_z$",
               'cmap': "RdBu",
               'vsymmetric':True,
               'vmin': -2.0000,
               'vmax':  2.0000,
                },
        'jy': {'title': r"$J_y$",
               'cmap': "RdBu",
               'vsymmetric':True,
               'vmin': -2.0000,
               'vmax':  2.0000,
                },
        'S':  {'title': r"$\log_{10}(S)$",
                'vmin':-4.0,
                'vmax': 0.0,
                'derived':True,
                'cmap':'inferno',
                'log':True,
                },
        'logS':  {'title': r"$(S_z - \langle S \rangle_z)/S_0 $",
                'vmin':  -1,
                'vmax':  +1,
                'derived':True,
                'cmap':'RdBu',
                'norm': SymLogNorm(1e-4, vmin=-1.0, vmax=1.0, base=10),
                #'norm': LogNorm(vmin=-4, vmax=1.0),
                },
}


#def read_var(f5F, var):
#    try:
#        val = f5F[var][:,:,0]
#    except:
#        nx = f5F['Nx'][()]
#        ny = f5F['Ny'][()]
#        nz = f5F['Nz'][()]
#        #print("reshaping 1D array into multiD with {} {} {}".format(nx,ny,nz))
#        val = f5F[var][:]
#        val = np.reshape(val, (nx, ny, nz))
#        val = val[:,:,0]
#    return val


def build_bgradb(f5F):

    bx = read_var(f5F, "bx")
    by = read_var(f5F, "by")
    bz = read_var(f5F, "bz")

    # Compute partial derivatives of bx and by
    dbx_dx, dbx_dy = np.gradient(bx)
    dby_dx, dby_dy = np.gradient(by)
    dbz_dx, dbz_dy = np.gradient(bz)
    bgradbx = bx*dbx_dx + by*dbx_dy
    bgradby = bx*dby_dx + by*dby_dy
    bgradbz = bx*dbz_dx + by*dbz_dy
    # Compute the divergence
    return np.sqrt(bgradbx*bgradbx + bgradby*bgradby + bgradbz*bgradbz)/(bx*bx + by*by + bz*bz)

def build_epar(f5F, conf):

    ex = pytools.read_h5_array(f5F, "ex", stride=conf.stride)
    ey = pytools.read_h5_array(f5F, "ey", stride=conf.stride)
    ez = pytools.read_h5_array(f5F, "ez", stride=conf.stride)

    bx = pytools.read_h5_array(f5F, "bx", stride=conf.stride)
    by = pytools.read_h5_array(f5F, "by", stride=conf.stride)
    bz = pytools.read_h5_array(f5F, "bz", stride=conf.stride)
    b  = np.sqrt(bx**2 + by**2 + bz**2)

    epar = (ex*bx + ey*by + ez*bz)/b

    return epar


# Poynting flux
def build_exy(f5F, conf):

    ex = pytools.read_h5_array(f5F, "ex", stride=conf.stride)
    ey = pytools.read_h5_array(f5F, "ey", stride=conf.stride)
    #ez = read_var(f5F, "ez")

    return np.sqrt(ex*ex + ey*ey)



# Poynting flux
def build_S(f5F):

    ex = read_var(f5F, "ex")
    ey = read_var(f5F, "ey")
    ez = read_var(f5F, "ez")

    bx = read_var(f5F, "bx")
    by = read_var(f5F, "by")
    bz = read_var(f5F, "bz")

    e  = np.sqrt(ex**2 + ey**2 + ez**2)
    b  = np.sqrt(bx**2 + by**2 + bz**2)

    # remove y slice averages
    #Sy[:, 50] = 10.0
    #bx0 = np.mean(bx, axis=0)

    #ex -= np.mean(ex, axis=0)
    #ey -= np.mean(ey, axis=0)
    #ez -= np.mean(ez, axis=0)

    #bx -= np.mean(bx, axis=0)
    #by -= np.mean(by, axis=0)
    #bz -= np.mean(bz, axis=0)

    # determinant
    #x  y  z
    #ex ey ez
    #bx by bz

    # cross product E x B
    Sx =  (ey*bz - ez*by)
    Sy = -(ex*bz - ez*bx)
    Sz =  (ex*by - ey*bx)

    #Sy -= np.mean(Sy, axis=0)

    #--------------------------------------------------
    def gaussian( x , s):
        return 1./np.sqrt( 2. * np.pi * s**2 ) * np.exp( -x**2 / ( 2. * s**2 ) )

    sig = 5.0
    w = int(3*sig)
    kernel = np.fromiter( (gaussian( x , sig ) for x in range( -w, w+1, 1 ) ), float )
    #print('kernel', kernel)
    #--------------------------------------------------

    S0 = np.zeros_like(ex) # convolved Poynting flux

    # 1D smearing along the height direction
    for i in range(np.shape(ex)[1]): 
        S0[:,i] = np.convolve( Sy[:,i], kernel, mode='same' )

    # and then in x (does not seem to give good results)
    #for j in range(np.shape(ex)[0]): 
    #    #S0[j,:] = np.convolve( Sy[j,:], kernel, mode='same' )
    #    S0[j,:] = np.convolve( S0[j,:], kernel, mode='same' )

    Sr = Sy - S0

    #ind = np.where(Sr < 0)
    #ind = np.where(Sr < 0)
    #Sr[ ind ] = 0.0

    return Sr


def plot2dpcap_single(
        ax, 
        var,
        info, 
        args_cli,
        title= None,
        vmin = None,
        vmax = None,
        cmap = None,
        clip = None,
        ):

    #--------------------------------------------------
    # unpack incoming arguments that modify defaults
    args = {}

    #general defaults
    for key in default_values:
        args[key] = default_values[key]

    #overwrite with turbulence defaults
    try:
        for key in default_turbulence_values[var]:
            args[key] = default_turbulence_values[var][key]
    except:
        pass

    #finally, overwrite with user given values
    for key in args:
        try:
            user_val = eval(key)
            if not(user_val == None):
                args[key] = user_val
                print("overwriting {} key with {}".format(key, user_val))
        except:
            pass

    print('--------------------------------------------------')
    print("             lap: {}".format(info['lap']))

    if args['file'] == 'flds':
        print("reading {}".format(info['fields_file']))
        f5F = h5.File(info['fields_file'],'r')

    if args['file'] == 'moms':
        print("reading {}".format(info['moms_file']))
        f5F = h5.File(info['moms_file'],'r')


    # normal singular variables
    if not(args['derived']):
        val = pytools.read_h5_array(f5F, var, stride=conf.stride)

        #val /= 0.2056 #TODO automatize #n_0
        #val /= 0.2056*0.2056 #TODO automatize #curr

    # composite quantities
    else:
        print("building composite variable")

        if var == "je":      val = build_je(f5F)
        #elif var == "ve":    val, vex, vey, vez = build_ve(f5F)
        elif var == "epar":  val = build_epar(f5F, conf)
        elif var == 'logrho':val = np.log10( pytools.read_h5_array(f5F, 'rho', stride=conf.stride) )
        elif var == 'lognx': val = np.log10( pytools.read_h5_array(f5F, 'densx', stride=conf.stride) )
        elif var == 'S':     val = np.log10(abs(build_S(f5F)))
        elif var == 'logS':  val = build_S(f5F) 
        elif var == 'exy':   val = build_exy(f5F, conf) 


    #--------------------------------------------------
    if do_fieldlines and args['file'] == 'flds':
        bx = pytools.read_h5_array(f5F, "bx", stride=conf.stride)
        by = pytools.read_h5_array(f5F, "by", stride=conf.stride)
        bz = pytools.read_h5_array(f5F, "bz", stride=conf.stride)

        ex = pytools.read_h5_array(f5F, "ex", stride=conf.stride)
        ey = pytools.read_h5_array(f5F, "ey", stride=conf.stride)
        ez = pytools.read_h5_array(f5F, "ez", stride=conf.stride)

        if conf.twoD:
            bz[:] = 0
            ez[:] = 0
            

    elif do_fieldlines and not(args['file'] == 'flds'):
        f5F1 = h5.File(info['fields_file'],'r')

        sys.exit() # TODO

        if conf.twoD:
            bx = read_var(f5F1, "bx")
            by = read_var(f5F1, "by")
            bz = 0
            ex = read_var(f5F1, "ex")
            ey = read_var(f5F1, "ey")
            bz = 0
        if conf.threeD:
            bx = read_var(f5F1, "bx")
            by = read_var(f5F1, "by")
            bz = read_var(f5F1, "bz")
            ex = read_var(f5F1, "ex")
            ey = read_var(f5F1, "ey")
            ez = read_var(f5F1, "ez")


    #--------------------------------------------------
    # get shape
    nx, ny, nz = np.shape(val)
    print("nx={} ny={}".format(nx, ny))

    # strip out the 3rd dimension
    val = val[:,:,0]


    # typical image extent
    #xmin = 0.0
    #ymin = 0.0
    #xmax = nx/info['skindepth']
    #ymax = ny/info['skindepth']

    # star coordinates
    if conf.twoD:
        xmin = -conf.Lx//2
        ymin = 0
        xmax = conf.Lx//2
        ymax = conf.Ly
    if conf.threeD and args_cli.view == 'side':
        xmin = -conf.Lx//2
        ymin = 0
        xmax = conf.Lx//2
        ymax = conf.Lz
    if conf.threeD and args_cli.view == 'top':
        xmin = -conf.Lx//2
        xmax = +conf.Lx//2

        ymin = -conf.Ly//2
        ymax = +conf.Ly//2

    xmin /= conf.rad_pcap
    ymin /= conf.rad_pcap
    xmax /= conf.rad_pcap
    ymax /= conf.rad_pcap


    #if winsorize
    if not(args['winsorize_min'] == None) or not(args['winsorize_max'] == None):
        wvmin, wvmax = np.quantile(val, [args['winsorize_min'], 1.0-args['winsorize_max']])

        if args['vmin'] == None:
            args['vmin'] = wvmin
        if args['vmax'] == None:
            args['vmax'] = wvmax
    else:
        # else set vmin and vmax using normal min/max
        if args['vmin'] == None:
            args['vmin'] = np.min(val)
        if args['vmax'] == None:
            args['vmax'] = np.max(val)

    # make color limits symmetric
    if args['vsymmetric']:
        vminmax = np.maximum( np.abs(args['vmin']), np.abs(args['vmax']) )
        args['vmin'] = -vminmax
        args['vmax'] =  vminmax

    # finally, re-check that user did not give vmin/vmax
    args['vmin'] = vmin if not(vmin == None) else args['vmin']
    args['vmax'] = vmax if not(vmax == None) else args['vmax']

    #nor the project default
    args['vmin'] = default_turbulence_values[var]['vmin'] if not(vmin == None) else args['vmin']
    args['vmax'] = default_turbulence_values[var]['vmax'] if not(vmax == None) else args['vmax']


    #--------------------------------------------------
    # mask bad values
    if args['log']:
        inds = np.where(val < args['vmin'])
        val[inds] = args['vmin']

        inds = np.isnan(val)
        val[inds] = args['vmin']


    #--------------------------------------------------
    # normalization
    norm = 1.0
    n0 = conf.ppc*2 #*conf.stride**2 #number density per pixel in n_0 
    qe = np.abs(conf.qe)
    me_per_qe = np.abs(conf.me) / qe #for electrons = 1
    deltax = 1.0/conf.c_omp #\Delta x in units of skin depth

    lenscale = conf.rad_pcap #conf.Nx*conf.NxMesh*deltax/conf.max_mode #(large-eddy size in units of skin depth)

    if var in ['bx', 'by', 'bz']:
        norm = conf.binit #(me_per_qe*conf.cfl**2)/deltax
    if var in ['ex','ey','ez', 'epar', 'exy']:
        #norm = conf.binit #(me_per_qe*conf.cfl**2)/deltax
        norm = qe*conf.nGJ*conf.rad_pcap
    if var in ['rho', 'logrho']:
        #print(qe,conf.stride)
        #norm = n0
        norm = conf.nGJ
    if var in ['jx', 'jy', 'jz']:
        norm = qe*n0*conf.cfl*conf.cfl
    if var == 'je':
        norm_E = (me_per_qe*conf.cfl**2)/deltax/lenscale
        norm_J = qe*n0*conf.cfl*conf.cfl
        norm = norm_E * norm_J

        #correct for stride size in e/b fields
        norm /= conf.stride**2
        norm /= 1.0e3
    if var == 'logS':
        norm = conf.cfl*(qe*conf.nGJ*conf.rad_pcap)**2

    if not(args['log']):
        val = val / norm
    else:
        val = val - np.log10( norm )

    print("norm factor: {}".format( norm ))
    print("value at the corner {} / mean val {}".format( val[0,0], np.mean(val)) )
    print("time of ", lap2time(info['lap'], conf))


    #--------------------------------------------------
    # print all options
    print("--- {} ---".format(var))
    for key in args:
        if not(key == None):
            print(" setting {}: {}".format(key, args[key]))

    #--------------------------------------------------
    #if conf.threeD:
    #    val = val.T #  TODO transpose of transpose to neutralize the flips; done only when slice files are used

    im = pytools.visualize.imshow(
           ax, 
           val,
           xmin, xmax, ymin, ymax,
           cmap = args['cmap'],
           vmin = args['vmin'],
           vmax = args['vmax'],
           clip = args['clip'],
           aspect=args['aspect'],
           norm = args['norm'],
           )


    #-------------------------------------------------- 
    if do_fieldlines:

        nx, ny, nz = np.shape(bx)
        xx = np.linspace(xmin,xmax, nx)
        yy = np.linspace(ymin,ymax, ny)
        #X,Y = np.meshgrid(xx, yy, indexing='ij')

        bx = bx[:,:,0].T
        by = by[:,:,0].T
        bz = bz[:,:,0].T

        ex = ex[:,:,0].T
        ey = ey[:,:,0].T
        ez = ez[:,:,0].T

        # 2D setup
        bperp = bx
        bpara = by

        # 3D options
        if args_cli.view == 'side' and conf.threeD: # on side view we plot bx and bz
            bperp = bx
            bpara = bz
        if args_cli.view == 'top' and conf.threeD: # on top view we plot bx and by
            bperp = bx
            bpara = by

        norm_epar = abs(conf.qe)*conf.nGJ*conf.rad_pcap
        epar = (ex*bx + ey*by + ez*bz)/np.sqrt(bx**2 + by**2 + bz**2)/norm_epar

        print("epar", np.min(epar), np.max(epar) )

        lw = np.maximum( 1.0*np.ones_like(epar), 20.0*abs(epar))

        #print("streamlines")
        #print(np.shape(X), np.shape(Y))
        #print(np.shape(xx), np.shape(yy))
        #print(np.shape(bx), np.shape(by))
        #print(np.shape(lw))

        #print(bx)#
        #print(by)

        splt = ax.streamplot(
                xx, yy,
                bperp, bpara,
                density = 2.0,
                color = 'w',
                linewidth= lw,
                arrowsize=0.3,
                #broken_streamlines=False,
                )

        splt.lines.set_alpha(0.2)
        splt.arrows.set_alpha(0.2)


    #--------------------------------------------------
    # draw surface 

    if conf.twoD:
        cenx   = conf.Lx//2 + 0.5
        ceny   = -conf.rad_star + conf.rad_curv_shift
        cenz   = 0 
    elif conf.threeD:
        cenx   = conf.Lx//2 + 0.5
        #ceny   = conf.Ly//2 + 0.5
        ceny   = -conf.rad_star + conf.rad_curv_shift

    #--------------------------------------------------
    if args_cli.view == 'side':

        xx = np.linspace(0, conf.Nx*conf.NxMesh, 100) - cenx
        #xx**2 + yy**2 = rad**2
        yy = np.sqrt( conf.rad_star**2 - xx**2 ) + ceny
        ax.plot(xx/conf.rad_pcap, yy/conf.rad_pcap, "k-", lw=0.8)

        yy2 = np.sqrt( (conf.rad_star + conf.height_atms)**2 - xx**2 ) + ceny
        ax.plot(xx/conf.rad_pcap, yy2/conf.rad_pcap, "k-", alpha=1.0, lw=0.8)

        # draw polar cap limits
        sint = conf.rad_pcap/conf.rad_star
        Rbc = conf.rad_star/sint/sint

        #eta = sint*sint # require this to hold
        #rad/Rbc = s**2
        #rad     = Rbc*s**2
        #sqrt(x^2 + y^2 +) + Rbc*s^2
        #x^2 + y^2 = (Rbc*s^2)^2
        #yy3 = np.sqrt( (Rbc*sint*sint)**2 - xx**2 )
        #print(yy3)
        #ax.plot(xx/conf.rad_pcap, yy3/conf.rad_pcap, "w-", alpha=0.3)


    #--------------------------------------------------
    if do_dark:
        plt.subplots_adjust(left=0.10, bottom=0.05, right=0.90, top=0.97)
        ax.set_xlabel(r"$r/R_{\mathrm{pc}}$")
        ax.set_ylabel(r"$z/R_{\mathrm{pc}}$")

    else:
        #ax.set_xlabel(r"$x$ $(c/\omega_p)$")
        #ax.set_ylabel(r"$y$ $(c/\omega_p)$")

        ax.set_xlabel(r"$r/R_{\mathrm{pc}}$")
        ax.set_ylabel(r"$z/R_{\mathrm{pc}}$")

        plt.subplots_adjust(left=0.15, bottom=0.10, right=0.87, top=0.97)
        draw_label(ax, 850, 970, "a)", conf)


    wskip = 0.2
    pad = 0.01
    pos = ax.get_position()
    #print(pos)

    axleft   = pos.x0
    axbottom = pos.y0
    axright  = pos.x0 + pos.width
    axtop    = pos.y0 + pos.height

    cax = plt.fig.add_axes([axright+pad, axbottom+wskip, 0.01, axtop-axbottom-2*wskip])

    cb = plt.fig.colorbar(
            im, 
            cax=cax, 
            orientation='vertical',
            ticklocation='right',
            norm=args['norm'],
            )
    #cb.set_label(args['title']) 
    cax.text(1.0, 1.09, args['title'], transform=cax.transAxes)

    #if do_dark:
    #    for axis in ['top', 'bottom', 'left', 'right']:
    #        cax.spines[axis].set_visible(False)
    #        #cax.spines[axis].set_linewidth(0)
    #        #cax.spines[axis].set_color('red')
    #ax.set_title(title)



    slap = str(info['lap']).rjust(5, '0')
    if do_dark:
        if args_cli.view == 'side':
            fname = fdir + var +'_{}.png'.format(slap)
        elif args_cli.view == 'top':
            fname = fdir + var +'_top_{}.png'.format(slap)
        plt.savefig(fname)
    else:
        if args_cli.view == 'side':
            fname = fdir + var +'_{}.pdf'.format(slap)
        elif args_cli.view == 'top':
            fname = fdir + var +'_top_{}.pdf'.format(slap)
        plt.savefig(fname)
    cb.remove()


#--------------------------------------------------
    
def build_info(fdir, lap, conf, args):

    info = {}
    info['lap'] = lap
    info['fields_file']   = fdir + 'flds_'+str(lap)+'.h5'
    info['moms_file']     = fdir + 'moms_'+str(lap)+'.h5'
    info['analysis_file'] = fdir + 'analysis_'+str(lap)+'.h5'
    info['particle_file'] = fdir + 'test-prtcls'
     
    if conf.threeD and args.view == 'side':
        info['fields_file'] = fdir + 'slices-xz_' + str(lap) + '.h5' 
    if conf.threeD and args.view == 'top':
        info['fields_file'] = fdir + 'slices-xy_' + str(lap) + '.h5'
    
    return info



do_dark = True
do_fieldlines = True

if __name__ == "__main__":

    if do_dark:
        plt.fig = plt.figure(1, figsize=(8,7), dpi=300)
        
        plt.rc('font', family='serif', size=7)
        plt.rc('xtick')
        plt.rc('ytick')
    else:
        plt.fig = plt.figure(1, figsize=(4,3.5), dpi=200)
        plt.rc('font',  family='sans-serif')
        #plt.rc('text',  usetex=True)
        plt.rc('xtick', labelsize=8)
        plt.rc('ytick', labelsize=8)
        plt.rc('axes',  labelsize=8)
    

    gs = plt.GridSpec(1, 1)
    gs.update(hspace = 0.0)
    gs.update(wspace = 0.0)
    
    axs = []
    axs.append( plt.subplot(gs[0,0]) )
    
    #ax.set_xlabel(r"$x$ $(c/\omega_p)$")
    #ax.set_ylabel(r"$y$ $(c/\omega_p)$")

    #--------------------------------------------------
    # command line driven version
    args = pytools.parse_args()
    conf = Configuration(args.conf_filename, do_print=False)
    fdir = conf.outdir + "/"

    print("plotting {}".format(fdir))
    
    fname_F = "flds"
    fname_A = "analysis"
    fname_P = "test-prtcls"

    if do_dark: plt.style.use('dark_background')


    # if lap is defined, only plot that one individual round
    if not(args.lap == None):

        #info = build_info(fdir, args.lap)
        info = build_info(fdir, args.lap, conf, args)
        info['skindepth'] = conf.c_omp/conf.stride

        plot2dpcap_single(axs[0], args.var, info, args)

    # else plot every file that there is
    # FIXME
    else:
        #files_F = get_file_list(fdir, fname_F)
        #files_A = get_file_list(fdir, fname_A)
        #files_P = get_file_list(fdir, fname_P)
        #for lap, f in enumerate(files_F):
        for lap in range(0, conf.Nt + 1, conf.interval):

            info = build_info(fdir, lap, conf, args)
            info['skindepth'] = conf.c_omp/conf.stride
            plot2dpcap_single(axs[0], args.var, info, args)


    #varp = args.var + '_p'
    #figg.savefig(fname = fdir + varp+ '_dlngamma.pdf')

    
