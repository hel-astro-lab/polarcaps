import numpy as np
import matplotlib.pyplot as plt
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
def colorbar(mappable, 
        loc="right", 
        orientation="vertical", 
        size="1%", 
        #pad=0.05, 
        pad=0.1, 
        ticklocation='right'):
        #loc="top", 
        #orientation="horizontal", 
        #size="8%", 
        #pad=0.03, 
        #ticklocation='top'):

    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes(loc, size=size, pad=pad)
    return cax, fig.colorbar(mappable, cax=cax, orientation=orientation, ticklocation=ticklocation)



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
}

default_turbulence_values = {
        # DONE
        'by': {'title': r"$B_z$",
               'cmap': "RdBu",
               'vsymmetric':True,
                },
        'ex': {'title': r"$E_r$",
               'cmap': "PRGn",
               'vsymmetric':True,
                },
        'epar': {'title': r"$E_\parallel$",
               'cmap': "PRGn",
               'vsymmetric':True,
               'derived':True,
                'vmin':-0.01,
                'vmax': 0.01,
                },
        'bz': {'title': r"$B_z$",
               'cmap': "RdBu",
               'vsymmetric':True,
                },
        'rho': {'title': r"$n/n_0$",
                'vmin': 000000.0,
                'vmax': 4.0,
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
}



def read_var(f5F, var):
    try:
        val = f5F[var][:,:,0]
    except:
        nx = f5F['Nx'][()]
        ny = f5F['Ny'][()]
        nz = f5F['Nz'][()]
        #print("reshaping 1D array into multiD with {} {} {}".format(nx,ny,nz))

        val = f5F[var][:]
        val = np.reshape(val, (nx, ny, nz))
        val = val[:,:,0]

    return val


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

def build_epar(f5F):

    ex = read_var(f5F, "ex")
    ey = read_var(f5F, "ey")
    ez = read_var(f5F, "ez")

    bx = read_var(f5F, "bx")
    by = read_var(f5F, "by")
    bz = read_var(f5F, "bz")
    b  = np.sqrt(bx**2 + by**2 + bz**2)

    epar = (ex*bx + ey*by + ez*bz)/b

    return epar



def plot2dpcap_single(
        ax, 
        var,
        info, 
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
    print("reading {}".format(info['fields_file']))

    f5F = h5.File(info['fields_file'],'r')

    # normal singular variables
    if not(args['derived']):
        val = read_var(f5F, var)

        #val /= 0.2056 #TODO automatize #n_0
        #val /= 0.2056*0.2056 #TODO automatize #curr

    # composite quantities
    else:
        print("building composite variable")

        if var == "je":     val = build_je(f5F)
        elif var == "ve":   val, vex, vey, vez = build_ve(f5F)
        elif var == "epar": val = build_epar(f5F)


    #--------------------------------------------------
    if do_fieldlines:
        bx = read_var(f5F, "bx")
        by = read_var(f5F, "by")


    #--------------------------------------------------
    # get shape
    nx, ny = np.shape(val)
    print("nx={} ny={}".format(nx, ny))

    #xmin = 0.0
    #ymin = 0.0
    #xmax = nx/info['skindepth']
    #ymax = ny/info['skindepth']

    xmin = -conf.Lx//2
    ymin = 0
    xmax = conf.Lx//2
    ymax = conf.Ly

    xmin /= conf.Rpc
    ymin /= conf.Rpc
    xmax /= conf.Rpc
    ymax /= conf.Rpc


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
    # normalization
    norm = 1.0
    n0 = conf.ppc*2 #*conf.stride**2 #number density per pixel in n_0 
    qe = np.abs(conf.qe)
    me_per_qe = np.abs(conf.me) / qe #for electrons = 1
    deltax = 1.0/conf.c_omp #\Delta x in units of skin depth

    lenscale = conf.Rpc #conf.Nx*conf.NxMesh*deltax/conf.max_mode #(large-eddy size in units of skin depth)

    if var in ['bx', 'by', 'bz']:
        norm = conf.binit #(me_per_qe*conf.cfl**2)/deltax
    if var in ['ex','ey','ez', 'epar']:
        norm = conf.binit #(me_per_qe*conf.cfl**2)/deltax

    if var == 'rho':
        #print(qe,conf.stride)
        norm = n0
    if var == 'jz':
        norm = qe*n0*conf.cfl*conf.cfl
    if var == 'je':
        norm_E = (me_per_qe*conf.cfl**2)/deltax/lenscale
        norm_J = qe*n0*conf.cfl*conf.cfl
        norm = norm_E * norm_J

        #correct for stride size in e/b fields
        norm /= conf.stride**2
        norm /= 1.0e3

    val = val / norm

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

    im = pytools.visualize.imshow(
           ax, 
           np.transpose(val), 
           xmin, xmax, ymin, ymax,
           cmap = args['cmap'],
           vmin = args['vmin'],
           vmax = args['vmax'],
           clip = args['clip'],
           aspect=args['aspect'],
           )


    #-------------------------------------------------- 
    if do_fieldlines:
        xx = np.linspace(xmin,xmax, nx)
        yy = np.linspace(ymin,ymax, ny)
        X, Y = np.meshgrid(xx, yy, indexing='ij')

        val2 = np.sqrt(bx**2 + by**2)
        lw = 1 #1.5*val2/val2.max()

        #print("streamlines")
        #print(np.shape(X), np.shape(Y))
        #print(np.shape(xx), np.shape(yy))
        #print(np.shape(bx), np.shape(by))
        #print(np.shape(val2))

        #print(bx)#
        #print(by)

        splt = ax.streamplot(
                #X,Y,
                xx, yy,
                bx, by,
                density = 1.0,
                color = 'r',
                linewidth= lw,
                arrowsize=0.3
                )

        splt.lines.set_alpha(0.3)
        splt.arrows.set_alpha(0.3)

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

    cb = plt.fig.colorbar(im, cax=cax, 
            orientation='vertical',
            ticklocation='right')
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
        fname = fdir + var +'_{}.png'.format(slap)
        plt.savefig(fname)
    else:
        fname = fdir + var +'_{}.pdf'.format(slap)
        plt.savefig(fname)
    cb.remove()


#--------------------------------------------------

def build_info(fdir, lap):
    info = {}
    info['lap'] = lap
    info['fields_file']   = fdir + 'fields_'+str(lap)+'.h5'
    info['analysis_file'] = fdir + 'analysis'+str(lap)+'.h5'
    
    return info
    
def quick_build_info(fdir, lap):
    info = {}
    info['lap'] = lap
    info['fields_file']   = fdir + 'flds_'+str(lap)+'.h5'
    info['analysis_file'] = fdir + 'analysis_'+str(lap)+'.h5'
    info['particle_file'] = fdir + 'test-prtcls'
    
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
        info = quick_build_info(fdir, args.lap)
        info['skindepth'] = conf.c_omp/conf.stride

        plot2dpcap_single(axs[0], args.var, info,)

    # else plot every file that there is
    # FIXME
    else:
        #files_F = get_file_list(fdir, fname_F)
        #files_A = get_file_list(fdir, fname_A)
        #files_P = get_file_list(fdir, fname_P)
        #for lap, f in enumerate(files_F):
        for lap in range(0, conf.Nt + 1, conf.interval):

            info = quick_build_info(fdir, lap)
            info['skindepth'] = conf.c_omp/conf.stride
            plot2dpcap_single(axs[0], args.var, info,)


    #varp = args.var + '_p'
    #figg.savefig(fname = fdir + varp+ '_dlngamma.pdf')

    
