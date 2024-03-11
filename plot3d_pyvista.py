import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import h5py as h5
import sys, os
import matplotlib.ticker as ticker
from matplotlib import cm
import argparse

# 3D
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable

# pytools bindings
import pytools
from init_problem import Configuration_Turbulence as Configuration

from pytools.pybox.box import Box
from pytools.pybox.box import axisEqual3D

import pyvista as pv


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
    'file':'flds',
    'log':False,
    'vmin':-1,
    'vmax':+1,
}


default_turbulence_values = {
        'rho': {'title': r"$n/n_0$",
                'vmin': 0.0,
                'vmax': 4.0,
                },
        'jz': {'title': r"$J_z$",
               'cmap': "RdBu",
               'vsymmetric':True,
               'vmin': -1.0000,
               'vmax':  1.0000,
                },
        'bz': {'title': r"$B_z$",
               'cmap': "RdBu",
               'vsymmetric':True,
                },
        'bperp': {'title': r"$B_\perp$",
               'cmap': "magma",
               'vmin': -1.0000,
               'vmax':  1.0000,
               'derived':True,
                },
        'bvec': {'title': r"$B$",
               'cmap': "RdBu",
               'vsymmetric':True,
               'vmin': -1.0000,
               'vmax':  1.0000,
               'derived':True,
                },
}


#--------------------------------------------------
# normalization
def get_normalization(var, conf):
    norm = 1.0
    n0 = conf.ppc*2 #*conf.stride**2 #number density per pixel in n_0 
    qe = np.abs(conf.qe)
    me_per_qe = np.abs(conf.me) / qe #for electrons = 1
    deltax = 1.0/conf.c_omp #\Delta x in units of skin depth

    lenscale = conf.rad_pcap #conf.Nx*conf.NxMesh*deltax/conf.max_mode #(large-eddy size in units of skin depth)

    if var in ['bx', 'by', 'bz']:
        norm = conf.binit #(me_per_qe*conf.cfl**2)/deltax
    if var in ['ex','ey','ez', 'epar']:
        #norm = conf.binit #(me_per_qe*conf.cfl**2)/deltax
        norm = qe*conf.nGJ*conf.rad_pcap
    if var in ['rho', 'logrho']:
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

    return norm


def read_full_box(outdir, fname_fld, var, lap, stride=1):
    fields_file = outdir + "/" + fname_fld + "_" + str(lap) + ".h5"
    f5 = h5.File(fields_file, "r")
    #return pytools.read_h5_array(f5, var, stride=stride)

    nx = f5['Nx'][()]
    ny = f5['Ny'][()]
    nz = f5['Nz'][()]

    val = f5[var][()]

    #print("reshaping 1D array of {} into multiD with {} {} {}".format(len(val), nx,ny,nz))

    val = np.reshape(val, (nz, ny, nx))
    #val = val.ravel(order='F').reshape((nx,ny,nz))

    # pyvista expects C ordering of arrays so have to fiddle 
    # here and manually reshape the array

    val = val.reshape((nx,ny,nz)) 

    return val


do_print = True

if __name__ == "__main__":

    args_cli = pytools.parse_args()
    conf = Configuration(args_cli.conf_filename, do_print=do_print)
    var = args_cli.var # NOTE uncomment for basic functionality 
    #var = 'bvec' # manually set the plotted variable

    # additional args
    #parser = argparse.ArgumentParser(description="3D visualization")
    #args_script = parser.parse_args()
    #print('args s: ', args_script)


    #general defaults
    args = {}
    for key in default_values:
        args[key] = default_values[key]

    #overwrite with turbulence defaults
    try:
        for key in default_turbulence_values[var]:
            args[key] = default_turbulence_values[var][key]
    except:
        pass


    print(conf.outdir)
    print("plotting {}".format(var))


    fname_fld = args['file']
    fname_prtcls = "test-prtcls"

    lap = args_cli.lap
    #lap = 2000

    # reading funciton for data
    def read_h5(outdir, fname, var, lap, stride=1):
        return read_full_box(outdir, fname, var, lap, stride=stride)

    #--------------------------------------------------
    print(fname_fld)
    s = 1
    rho= read_h5(conf.outdir, fname_fld, "rho", lap,stride=s)

    jx = read_h5(conf.outdir, fname_fld, "jx", lap, stride=s)
    jy = read_h5(conf.outdir, fname_fld, "jy", lap, stride=s)
    jz = read_h5(conf.outdir, fname_fld, "jz", lap, stride=s)

    ex = read_h5(conf.outdir, fname_fld, "ex", lap, stride=s)
    ey = read_h5(conf.outdir, fname_fld, "ey", lap, stride=s)
    ez = read_h5(conf.outdir, fname_fld, "ez", lap, stride=s)

    bx = read_h5(conf.outdir, fname_fld, "bx", lap, stride=s)
    by = read_h5(conf.outdir, fname_fld, "by", lap, stride=s)
    bz = read_h5(conf.outdir, fname_fld, "bz", lap, stride=s)

    # normalize
    rho/= get_normalization('rho',conf)

    jx /= get_normalization('jx', conf)
    jy /= get_normalization('jy', conf)
    jz /= get_normalization('jz', conf)

    bx /= get_normalization('bx', conf)
    by /= get_normalization('by', conf)
    bz /= get_normalization('bz', conf)

    ex /= get_normalization('ex', conf)
    ey /= get_normalization('ey', conf)
    ez /= get_normalization('ez', conf)

    print(np.shape(rho))
    nx, ny, nz = np.shape(rho)


    #--------------------------------------------------
    # create pyvista object

    dx = conf.stride/conf.c_omp # skindepth resolution
    origin = 0,0,0

    Lx = nx*dx
    Ly = ny*dx
    Lz = nz*dx

    #mesh = pv.ImageData(
    #        dimensions=(nx, ny, nz), 
    #        spacing=(dx, dx, dx), 
    #        origin=origin)
    #x = mesh.points[:, 0]
    #y = mesh.points[:, 1]
    #z = mesh.points[:, 2]

    pv.global_theme.background = 'black'

    pv.global_theme.allow_empty_mesh=True

    p = pv.Plotter(
            #lighting='three lights'
            off_screen=True,
            window_size = [1280, 2560],
            )


    #--------------------------------------------------
    # star

    #--------------------------------------------------
    if True:
        cenx = conf.Lx//2 + 0.5
        ceny = conf.Ly//2 + 0.5
        cenz = -conf.rad_star + conf.rad_curv_shift

        sphere = pv.Sphere(
                radius=conf.rad_star,
                center=(cenx, ceny, cenz),
                theta_resolution=100,
                phi_resolution=100,
                )

        p.add_mesh(sphere, 
                   color='lightblue', 
                   #color='grey', 
                   #show_edges=True,
                   )


    if False: # jpar slices / projections

        mesh = pv.ImageData(
            dimensions=(nx, ny, nz), 
            spacing=(dx, dx, dx), 
            origin=origin)

        #--------------------------------------------------
        # j scalar field
        #jabs = np.sqrt(jx**2 + jy**2 + jz**2)
        #jabs = np.sqrt(jz**2)
        var = jz

        mesh['scalar'] = np.ravel(var)

        single_slice = mesh.slice(
                normal=[0, 0, 1],
                origin=[0, 0, 0.2*conf.Lz],
                )

        #surface = pv.Disc(
        #        center=(Lx//2,Ly//2,0.95*Lz), 
        #        inner=0,
        #        outer=conf.Lx//2,
        #        r_res=5,
        #        c_res=30,
        #        normal=(0,0,1),
        #        )
        #mesh.compute_implicit_distance(surface, inplace=True)
        #inner = mesh.threshold(0.0, scalars="implicit_distance", invert=True)

        p.add_mesh(
                single_slice, 
                #inner,
                show_scalar_bar=False,
                cmap='RdBu',
                clim=(-1.0, 1.0),
                opacity=0.9,
                   )

        #for z in np.linspace(0.1, 0.9, 10):
        #    single_slice = mesh.slice(
        #        normal=[0, 0, 1],
        #        origin=[0, 0, z*conf.Lz],
        #        )
        #    p.add_mesh(single_slice, 
        #        cmap='RdBu',
        #        clim=(-1.0, 1.0),
        #        opacity=0.8,
        #           )




    #--------------------------------------------------
    # b field
    if True: # streamlines

        m1 = pv.ImageData(
            dimensions=(nx, ny, nz), 
            spacing=(dx, dx, dx), 
            origin=origin)

        vectors = np.empty((m1.n_points, 3))
        vectors[:,0] = np.ravel(bx)
        vectors[:,1] = np.ravel(by)
        vectors[:,2] = np.ravel(bz)

        m1['vectors'] = vectors

        print(m1)

        #p.add_mesh(m1.outline(), color="k")

        #mesh = pv.ImageData(
        #    dimensions=(nx, ny, nz), 
        #    spacing=(dx, dx, dx), 
        #    origin=origin)
        b2 = np.sqrt( bx**2 + by**2 + bz**2 )
        print('bvec min/max', np.min(b2), np.max(b2), np.mean(b2))
        m1['scalar'] = 0.1*np.ravel(b2)
        #m1['scalar'] = 1.0*np.ones_like(b2).ravel() 


        midx = Lx//2
        midy = Ly//2
        midz = 0 #Lz//2

        #seed = pv.Plane(center=(midx,midy,midz), 
        #                i_size=conf.rad_pcap, 
        #                j_size=conf.rad_pcap,
        #                i_resolution=8,
        #                j_resolution=8,
        #                direction=(0,0,1),
        #                )

        seed = pv.Disc(
                center=(midx,midy,midz), 
                inner=0,
                outer=conf.rad_pcap,
                r_res=5,
                c_res=30,
                normal=(0,0,1),
                )

        #p.add_mesh(seed.outline(), color='r')

        stream = m1.streamlines_from_source(
            seed,
            vectors='vectors',
            max_time=2e3,
            initial_step_length=0.1,
            integration_direction='both',
            #integration_direction='forward',
            )


        #stream, src = m1.streamlines(
        #    'vectors', 
        #    return_source=True, 
        #    terminal_speed=0.0, 
        #    n_points=200, 
        #    source_radius=200,
        #    source_center=(0, 0, 0),
        #)

        #p.add_mesh(src)
        p.add_mesh(
                stream.tube(
                    radius=0.1,
                    scalars="scalar",
                    radius_factor=10.0,
                    ),
                show_scalar_bar=False,
                cmap='reds',
                #cmap='inferno',
                ambient=0.5,
                opacity=0.2,
                clim=(0.0, 2.0),
                   )


    #--------------------------------------------------
    # e field
    if False: # streamlines

        m1 = pv.ImageData(
            dimensions=(nx, ny, nz), 
            spacing=(dx, dx, dx), 
            origin=origin)

        vectors = np.empty((m1.n_points, 3))
        vectors[:,0] = np.ravel(ex)
        vectors[:,1] = np.ravel(ey)
        vectors[:,2] = np.ravel(ez)

        m1['vectors'] = vectors

        print(m1)

        #p.add_mesh(m1.outline(), color="k")

        e2 = np.sqrt( ex**2 + ey**2 + ez**2 )
        print('evec min/max', np.min(e2), np.max(e2), np.mean(e2))
        m1['scalar'] = 0.1*np.ravel(e2)

        Lx = nx*dx
        Ly = ny*dx
        Lz = nz*dx

        # xy plane
        midx = Lx//2
        midy = Ly//2
        midz = int(0.1*Lz)
        seed = pv.Plane(center=(midx,midy,midz), 
                        i_size=Lx, 
                        j_size=Ly,
                        i_resolution=8,
                        j_resolution=8,
                        direction=(1,1,0),
                        )
        # xz plane
        #midx = Lx//2
        #midy = Ly//2
        #midz = Lz//2
        #seed = pv.Plane(center=(midx,midy,midz), 
        #                i_size=Lx, 
        #                j_size=Ly,
        #                i_resolution=15,
        #                j_resolution=15,
        #                direction=(0,1,0),
        #                )

        stream = m1.streamlines_from_source(
            seed,
            vectors='vectors',
            max_time=2e3,
            initial_step_length=0.1,
            integration_direction='both',
            #integration_direction='forward',
            )

        #p.add_mesh(src)
        p.add_mesh(
                stream.tube(
                    radius=0.1,
                    scalars="scalar",
                    radius_factor=10.0,
                    ),
                show_scalar_bar=False,
                cmap='Blues',
                #cmap='inferno',
                ambient=0.5,
                clim=(0.0, 1.0),
                opacity=0.3,
                   )


    if True: # volume rendering
        mesh = pv.ImageData(
            dimensions=(nx, ny, nz), 
            spacing=(dx, dx, dx), 
            origin=origin)


        #--------------------------------------------------
        # j scalar field
        clim = (-4, 0.3)

        #ind = np.where(rho < 2.0)
        #rho[ind] = 0.0

        print('rho min/max', np.min(rho), np.max(rho), ' mean', np.mean(rho))

        mesh['rho'] = np.log10( np.ravel(rho) )
        #p.add_mesh( mesh.outline_corners(), color="k" )

        #mesh = mesh.gaussian_smooth(std_dev=0.5)

        ops = np.zeros(256)
        ops[0:192] = 2.0
        ops[192:]  = np.linspace(2, 128, 64)
        #ops[64:]  = np.linspace(0.0, 256, 192)

        p.add_volume(mesh,
                     scalars='rho',
                     #opacity=ops, #'geom',
                     opacity='geom',
                     clim=clim,
                     show_scalar_bar=False,
                     cmap='inferno',
                     shade=False,
                    )


    if False: # volume rendering
        mesh = pv.ImageData(
            dimensions=(nx, ny, nz), 
            spacing=(dx, dx, dx), 
            origin=origin)

        #--------------------------------------------------
        # e scalar field
        clim = (-3, -0.3)

        #b = np.sqrt(bx**2 + by**2 + bz**2)
        #epar = (ex*bx + ey*by + ez*bz)/b
        #v = abs(epar)

        b = np.sqrt(bx**2 + by**2 + bz**2)
        jpar = (jx*bx + jy*by + jz*bz)/b
        v = abs(jpar)

        #e = np.sqrt(ex**2 + ey**2 + ez**2)
        #v = abs(e)

        print('e min/max', np.min(v), np.max(v), ' mean', np.mean(v))

        mesh['v'] = np.log10( np.ravel(v) )
        #mesh = mesh.gaussian_smooth(std_dev=0.5)

        ops = np.zeros(256)
        ops[0:192] = 2.0
        ops[192:]  = np.linspace(2, 128, 64)
        #ops[64:]  = np.linspace(0.0, 256, 192)

        p.add_volume(mesh,
                     scalars='v',
                     opacity=ops, #'geom',
                     #opacity='linear',
                     clim=clim,
                     show_scalar_bar=False,
                     cmap='Blues',
                     shade=False,
                    )


    if False: # rho contours rendering
        mesh = pv.ImageData(
            dimensions=(nx, ny, nz), 
            spacing=(dx, dx, dx), 
            origin=origin)

        #--------------------------------------------------
        # rho
        print('min/max', np.min(rho), np.max(rho), ' mean', np.mean(rho))

        mesh['scalar'] = np.log10( np.ravel(rho) )
        mesh = mesh.gaussian_smooth(std_dev=0.6)

        #p.add_mesh( mesh.outline_corners(), color="k")

        contours = mesh.contour(np.linspace(-2, 1, 5))
        #contours = mesh.contour(np.logspace(-2, 1.0, 3))
        #contours = contours.smooth_taubin(50)

        p.add_mesh(contours,
                    opacity=0.1,
                    clim=(-4, 2),
                    show_scalar_bar=False,
                     cmap='inferno',
                    )

        hc  = mesh.contour([1.0])
        ##hc = hc.smooth_taubin(10, pass_band=10.0)
        p.add_mesh(hc, 
                   opacity=0.8,
                   clim=(-4, 2),
                   show_scalar_bar=False,
                   cmap='inferno',
                   )


    if False: # try isocountours rendering
        mesh = pv.ImageData(
            dimensions=(nx, ny, nz), 
            spacing=(dx, dx, dx), 
            origin=origin)

        #--------------------------------------------------
        # j scalar field
        jz = np.sign(jz)*np.sqrt( jx**2 + jy**2 + jz**2 )

        clim = (-1.5, 1.5)
        jz = np.sign(jz)*np.sqrt( jx**2 + jy**2 + jz**2 )
        mesh['jz'] = np.ravel(jz)
        p.add_mesh( mesh.outline_corners(), color="k" )

        ops = np.zeros(256)
        ops[0:128] = np.linspace(128, 0, 128)
        ops[128:]  = np.linspace(0, 128, 128)

        #ops = np.zeros(256)
        #ops[0:64]     = np.linspace(128.0, 0.0,   64)
        #ops[64:128]   = 0.0
        #ops[128:192:] = 0.0
        #ops[192:]     = np.linspace(0.0,   128.0, 64)

        print(ops)

        p.add_volume(mesh,
                     scalars='jz',
                     opacity=ops, #'linear',
                     clim=clim,
                     show_scalar_bar=False,
                     cmap='seismic',
                     #cmap='RdBu',
                     shade=False,
                    )



    if False: # try isocountours rendering
        mesh = pv.ImageData(
            dimensions=(nx, ny, nz), 
            spacing=(dx, dx, dx), 
            origin=origin)

        #--------------------------------------------------
        # j scalar field
        jz = np.sign(jz)*np.sqrt( jx**2 + jy**2 + jz**2 )
        mesh['scalar'] = np.ravel(jz)
        #mesh = mesh.gaussian_smooth(std_dev=0.8)

        p.add_mesh( mesh.outline_corners(), color="k")

        vv = 0.8
        clim = (-1.0, 1.0)

        contours = mesh.contour([0.5*vv, 0.7*vv, 0.9*vv])
        p.add_mesh(contours,
                    opacity=0.3,
                    clim=clim,
                    show_scalar_bar=False,
                    cmap='RdBu',
                    )

        contours = mesh.contour([-0.5*vv, -0.7*vv, -0.9*vv])
        p.add_mesh(contours,
                    opacity=0.3,
                    clim=clim,
                    show_scalar_bar=False,
                    cmap='RdBu',
                    )

        hcp = mesh.contour([+vv])
        p.add_mesh(hcp, 
                   opacity=1.0,
                    clim=clim,
                    show_scalar_bar=False,
                    cmap='RdBu',
                   )

        hcm = mesh.contour([-vv])
        p.add_mesh(hcm, 
                   opacity=1.0,
                    clim=clim,
                   show_scalar_bar=False,
                   cmap='RdBu',
                   )


    if False: # try isocountours rendering
        mesh = pv.ImageData(
            dimensions=(nx, ny, nz), 
            spacing=(dx, dx, dx), 
            origin=origin)

        #--------------------------------------------------
        # j scalar field
        jabs = np.sqrt(jx**2 + jy**2 + jz**2)
        mesh['scalar'] = np.ravel(jabs)

        #mesh = mesh.gaussian_smooth(std_dev=0.8)

        p.add_mesh( mesh.outline_corners(), color="k")

        contours = mesh.contour(np.linspace(0.2, 1.0, 5))
        #contours = contours.smooth_taubin(50)
        p.add_mesh(contours,
                    opacity=0.1,
                    clim=(0.3,2.0),
                    show_scalar_bar=False,
                    )

        hc  = mesh.contour([1.0])
        #hc = hc.smooth_taubin(10, pass_band=10.0)

        p.add_mesh(hc, 
                   opacity=1.0,
                   clim=(0.3,2.0),
                   show_scalar_bar=False,
                   )


        #contours = mesh.contour([1.0], method='marching_cubes')
        ##co = grid.contour([0], values, method='flying_edges')
        #dist = np.linalg.norm(mesh.points, axis=1)

        #mesh.plot(
        #        scalars=dist, 
        #        smooth_shading=True, 
        #        specular=1, 
        #        cmap="plasma", 
        #        show_scalar_bar=False)




    #cpos = [(380, 380, 380),
    #        (319.27461767196655, 319.5708401799202, 319.71321415901184),
    #        (0.0, 0.0, 1.0)]


    if args_cli.view == 'side': 
        cpos = [(400, 400, 278),
                (360, 360, 250),
               (-0.3, -0.3, 0.9)]

    if args_cli.view == 'topside': 
        cpos = [(250, 243, 444),
                (225, 218, 392),
                (-0.59, -0.578, 0.56)]

    if args_cli.view == 'top': # top
        cpos = [(123, 91, 479),
            (120, 90, 440),
            (-0.704, -0.704, 0.092)]


    p.camera_position = cpos

    p.enable_depth_peeling(100)
    #p.enable_anti_aliasing('msaa')
    p.enable_ssao(kernel_size=32)

    #p.save_graphic("3d_jz_b.png", title="", raster=True, painter=True)

    #cpos = p.show(return_cpos=True, auto_close=False)
    slap = str(lap).rjust(5, '0')
    p.screenshot(conf.outdir + "/" + "3d_" + args_cli.var + "_" + str(args_cli.view) + "_" + slap + ".png", scale=1)


    #cpos = p.camera.position
    print('camera pos:', cpos)

    print('azimuth', p.camera.azimuth)
    print('elevation', p.camera.elevation)
    print('view_angle', p.camera.view_angle)
    print('up', p.camera.up)

    if False:
        p.open_movie(conf.outdir + "/" + "3d_jz_b.mp4", quality=8)

        for az in np.linspace(0, 360, 360):
            p.camera.azimuth = az
            #p.screenshot("3d_jz_b_{}.png".format(int(az)) )
            p.write_frame()

        #for el in np.linspace(0, 50, 100):
        #    p.camera.elevation = el
        #    p.write_frame()

        ## zoom in
        #for va in np.linspace(30, 10, 100):
        #    p.camera.view_angle = va
        #    p.write_frame()

        #for pause in range(20):
        #    p.write_frame()

        ## zoom out
        #for va in np.linspace(10, 30, 100):
        #    p.camera.view_angle = va
        #    p.write_frame()

        #for el in np.linspace(50, 0, 100):
        #    p.camera.elevation = el
        #    p.write_frame()

        p.close()












