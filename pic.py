# -*- coding: utf-8 -*-

from mpi4py import MPI
import numpy as np
import sys, os
from numpy import log10, sin, cos, tan, exp

# runko + auxiliary modules
import pytools  # runko python tools
import pyrunko

# problem specific modules
from init_problem import Configuration_Turbulence as Configuration
from init_problem import velocity_profile
from init_problem import density_profile
from init_problem import weigth_profile

from qed_toolset import QEDToolset


#--------------------------------------------------
live_plot = True
rnd_seed_default = 1
np.random.seed(rnd_seed_default)  # global simulation seed

import warnings
#warnings.filterwarnings("once") # suppress warnings after reporting them once
warnings.filterwarnings('ignore')


#--------------------------------------------------
# Field initialization (guide field)
def insert_em_fields(grid, conf, do_initialization=True):

    for tile in pytools.tiles_all(grid):
        g = tile.get_grids(0)

        ii,jj,kk = tile.index if conf.threeD else (*tile.index, 0)


        # insert values into Yee lattices; includes halos from -3 to n+3
        if do_initialization:
            for n in range(-3, conf.NzMesh + 3):
                for m in range(-3, conf.NyMesh + 3):
                    for l in range(-3, conf.NxMesh + 3):

                        if not(conf.use_maxwell_split): # if no static component

                            # get global coordinates
                            #iglob, jglob, kglob = pytools.ind2loc((ii, jj, kk), (l, m, n), conf)
                            #r = np.sqrt(iglob ** 2 + jglob ** 2 + kglob ** 2)

                            g.ex[l,m,n] = 0.0
                            g.ey[l,m,n] = 0.0
                            g.ez[l,m,n] = 0.0

                            g.bx[l,m,n] = 0.0
                            g.by[l,m,n] = 0.0 
                            g.bz[l,m,n] = 0.0 #conf.binit

                        elif conf.use_maxwell_split: # static component
                            1
                            # TODO
    return
        
#-------------------------------------------------- 
#-------------------------------------------------- 
#-------------------------------------------------- 
if __name__ == "__main__":


    # --------------------------------------------------
    # initialize auxiliary tools
    sch  = pytools.Scheduler() # Set MPI aware task scheduler
    args = pytools.parse_args() # parse command line arguments
    tplt = pytools.TerminalPlot(13, 13) # terminal plotting tool


    # create conf object with simulation parameters based on them
    conf = Configuration(args.conf_filename, do_print=sch.is_master)
    sch.conf = conf # remember to update scheduler

    if conf.mpi_task_mode: sch.switch_to_task_mode()

    # --------------------------------------------------
    # Timer for profiling
    timer = pytools.Timer()
    timer.start("total")
    timer.start("init")
    timer.do_print = sch.is_example_worker
    timer.verbose = 0  # 0 normal; 1 - debug mode

    sch.timer = timer # remember to update scheduler


    #--------------------------------------------------
    toolset = QEDToolset(conf)

    if live_plot and sch.is_master:
        import matplotlib
        import matplotlib.pyplot as plt
        from matplotlib.colors import Normalize
        from matplotlib.cm import get_cmap
        from matplotlib import colorbar

        fig = plt.figure(1, figsize=(12.0, 12.0))

        #plt.style.use('classic')    
        plt.rc('xtick', top = True)
        plt.rc('ytick', right = True)

        plt.rc('font',  family='sans-serif', size=8)
        plt.rc('text',  usetex=False)
        plt.rc('xtick', labelsize=7)
        plt.rc('ytick', labelsize=7)
        plt.rc('axes',  labelsize=8)
        plt.rc('legend',  handlelength=4.0)

        nrow_fig = 5
        ncol_fig = 4

        gs = plt.GridSpec(nrow_fig, ncol_fig)
        gs.update(hspace = 0.3)
        gs.update(wspace = 0.3)
        
        axs = np.empty( (nrow_fig,ncol_fig), dtype=object)

        for j in range(ncol_fig):
            for i in range(nrow_fig):
                axs[i,j] = plt.subplot(gs[i,j])
                axs[i,j].minorticks_on()
                axs[i,j].set_yscale('log')

        axleft    = 0.10
        axbottom  = 0.13
        axright   = 0.97
        axtop     = 0.97

        pos1 = axs[0,0].get_position()
        axwidth  = axright - axleft
        axheight = (axtop - axbottom)*0.02
        axpad = 0.01

        # colormap
        Nsteps = conf.Nt + 1
        norm = matplotlib.colors.Normalize(vmin=0.0, vmax=Nsteps)
        cmap = matplotlib.colormaps['turbo_r']


        # photon spectra
        for (i,j) in [ (0,0), (4,0), (4,2) ]:
            axs[i,j].set_title('photons')
            axs[i,j].set_xlabel(r'$\log x$')
            axs[i,j].set_ylabel(r'$x \ell_x \propto x d n_p/d \log x$') # lx
            axs[i,j].set_xlim(toolset.xxlims)
            axs[i,j].set_ylim(toolset.xylims)

        #axs[0,2].set_title('esc photons')
        #axs[0,2].set_xlabel(r'$\log x$')
        #axs[0,2].set_ylabel(r'$x \ell_x \propto x d n_p/d \log x$') # lx
        #axs[0,2].set_xlim(toolset.xxlims)
        #axs[0,2].set_ylim(toolset.xylims)

        # pair spectra
        for (i,j) in [ (0,1), (0,2), (4,1),]:
            axs[i,j].set_xlabel(r'$\log p$')
            axs[i,j].set_ylabel(r'$p d\tau/d p$ $\propto f_e$')
            axs[i,j].set_ylim(toolset.pylims)
            axs[i,j].set_xlim(toolset.pxlims)
            axs[i,j].set_title('pairs')

        axs[0,1].set_title('e-')
        axs[0,2].set_title('e+')

        # photon LP distribution
        #axs[1,0].set_title('photon weights')
        axs[1,0].set_xlabel(r'$\log x$')
        axs[1,0].set_ylabel(r'$d N_x/d \log x$') 
        axs[1,0].set_xlim(toolset.xxlims)

        # pair LP distribution
        #axs[1,1].set_title('pair weights')
        axs[1,1].set_xlabel(r'$\log p$')
        axs[1,1].set_ylabel(r'$d N_p/d \log p$') 
        axs[1,1].set_xlim(toolset.pxlims)

        # photon w distribution
        axs[2,0].set_xlabel(r'$\log w$')
        axs[2,0].set_ylabel(r'$d N_x/d \log w$') 
        axs[2,0].set_xlim(toolset.wxlims)

        # pair w distribution
        axs[2,1].set_xlabel(r'$\log w$')
        axs[2,1].set_ylabel(r'$d N_p/d \log w$') 
        axs[2,1].set_xlim(toolset.wxlims)

        # optical depth
        axs[0,3].set_xlabel(r'lap')
        axs[0,3].set_ylabel(r'$\tau$') 
        #axs[0,3].set_yscale('linear')
        axs[0,3].set_yscale('log')
        axs[0,3].set_ylim((1e-4, 1.0))

        axs[1,2].set_xlabel(r'$t$ ($H/c$)')
        axs[1,2].set_ylabel(r'$U_x/U_\pm$') 
        axs[1,2].set_ylim((1e-2, 1e0))

        axs[1,3].set_xlabel(r'$t$ ($H/c$)')
        axs[1,3].set_ylabel(r'energy $\ell$') 
        axs[1,3].set_ylim((1e-2, 1e4))

        axs[2,2].set_xlabel(r'lap')
        axs[2,2].set_ylabel(r'$n_p/n_0$') 
        axs[2,2].set_ylim((1e-3, 1e1))

        axs[2,3].set_xlabel(r'lap')
        axs[2,3].set_ylabel(r'$N_p/N_0$') 
        axs[2,3].set_ylim((1e-3, 1e1))

        # 2d histogram
        axs[3,0].set_xscale("linear")
        axs[3,0].set_yscale("linear")
        axs[3,0].set_xlabel("log w")
        axs[3,0].set_ylabel("log x")
        axs[3,0].set_xlim(toolset.wxlims)
        axs[3,0].set_ylim(toolset.xxlims)
        im30 = axs[3,0].imshow( np.zeros((toolset.Nhist, toolset.Nhist)),
                         extent=[toolset.wxlims[0], toolset.wxlims[1], 
                                 toolset.xxlims[0], toolset.xxlims[1]],
                         origin='lower',
                         cmap='turbo',
                         aspect='auto',
                         interpolation='nearest',
                         vmin=-1,
                         vmax=6,
                         )


        # height vs ene; photons
        axs[3,1].set_title("ph height")
        axs[3,1].set_xscale("linear")
        axs[3,1].set_yscale("linear")
        axs[3,1].set_xlabel("height")
        axs[3,1].set_ylabel("log x")
        axs[3,1].set_xlim(toolset.hhlims)
        #axs[3,1].set_ylim(toolset.xxlims) # TODO
        axs[3,1].set_ylim((-toolset.Nhist, toolset.Nhist))
        im31 = axs[3,1].imshow( np.zeros((toolset.Nhist, 2*toolset.Nhist)),
                         extent=[toolset.hhlims[0], toolset.hhlims[1], 
                                 -toolset.Nhist, toolset.Nhist], # TODO
                               #toolset.xxlims[0], toolset.xxlims[1]],
                         origin='lower',
                         cmap='turbo',
                         aspect='auto',
                         interpolation='nearest',
                         vmin=np.log10( toolset.xylims[0]),
                         vmax=np.log10( toolset.xylims[1]),)

        # height vs ene; photons
        axs[3,2].set_title("e- height")
        axs[3,2].set_xscale("linear")
        axs[3,2].set_yscale("linear")
        axs[3,2].set_xlabel("height")
        axs[3,2].set_ylabel("log p")
        axs[3,2].set_xlim(toolset.hhlims)
        #axs[3,2].set_ylim(toolset.pxlims) # TODO
        axs[3,2].set_ylim((-toolset.Nhist, toolset.Nhist))
        im32 = axs[3,2].imshow( np.zeros((toolset.Nhist, 2*toolset.Nhist)),
                         extent=[toolset.hhlims[0], toolset.hhlims[1], 
                                 -toolset.Nhist, toolset.Nhist], # TODO
                               #toolset.pxlims[0], toolset.pxlims[1]],
                         origin='lower',
                         cmap='turbo',
                         aspect='auto',
                         interpolation='nearest',
                         vmin=np.log10(toolset.pylims[0]),
                         vmax=np.log10(toolset.pylims[1]),)

        axs[3,3].set_title("e+ height")
        axs[3,3].set_xscale("linear")
        axs[3,3].set_yscale("linear")
        axs[3,3].set_xlabel("height")
        axs[3,3].set_ylabel("log p")
        axs[3,3].set_xlim(toolset.hhlims)
        axs[3,3].set_ylim((-toolset.Nhist, toolset.Nhist))
        #axs[3,3].set_ylim(toolset.pxlims) # TODO
        im33 = axs[3,3].imshow( np.zeros((toolset.Nhist, 2*toolset.Nhist)),
                         extent=[toolset.hhlims[0], toolset.hhlims[1], 
                                 -toolset.Nhist, toolset.Nhist], # TODO
                                 #toolset.pxlims[0], toolset.pxlims[1]],
                         origin='lower',
                         cmap='turbo',
                         aspect='auto',
                         interpolation='nearest',
                         vmin=np.log10(toolset.pylims[0]),
                         vmax=np.log10(toolset.pylims[1]),)


        #axs[3,3].set_xlabel(r'lap')
        #axs[3,3].set_ylabel(r'$\ell_\mathrm{inj}$, $\ell_\mathrm{esc}$') 
        #axs[3,3].set_ylim((1e0, 1e3))

        axs[4,3].set_xlabel(r'lap')
        axs[4,3].set_ylabel(r'$\ell_\mathrm{inj}/\ell_\mathrm{out}$') 
        axs[4,3].set_yscale('linear')
        axs[4,3].set_ylim((0.0, 4.0))



    # --------------------------------------------------
    # create output folders
    if sch.is_master: pytools.create_output_folders(conf)

    # --------------------------------------------------
    # load runko
    if conf.threeD: # 3D modules
        import pycorgi.threeD as pycorgi   # corgi   bindings
        import pyrunko.pic.threeD as pypic # pic c++ bindings
        import pyrunko.emf.threeD as pyfld # fld c++ bindings
        import pyrunko.qed.threeD as pyqed # qed c++ bindings
    elif conf.twoD: # 2D modules
        import pycorgi.twoD as pycorgi     # corgi   bindings
        import pyrunko.pic.twoD as pypic   # pic c++ bindings
        import pyrunko.emf.twoD as pyfld   # fld c++ bindings
        import pyrunko.qed.twoD as pyqed   # qed c++ bindings



    # --------------------------------------------------
    # setup grid
    grid = pycorgi.Grid(conf.Nx, conf.Ny, conf.Nz)
    grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax, conf.zmin, conf.zmax)
    sch.grid = grid # remember to update scheduler


    # compute initial mpi ranks using Hilbert's curve partitioning
    if conf.use_injector:
        pytools.load_catepillar_track_mpi(grid, conf.mpi_track, conf) 
    else:
        #pytools.load_mpi_y_strides(grid, conf) #equal x stripes

        #if grid.size() > 24: # big run
        #    pytools.balance_mpi_3D_rootmem(grid, 4) # advanced Hilbert curve optimization (leave 4 first ranks empty)
        #else:
        pytools.balance_mpi(grid, conf) #Hilbert curve optimization 


    # load pic tiles into grid
    pytools.pic.load_tiles(grid, conf)
    #load_damping_tiles(grid, conf)

    # --------------------------------------------------
    # simulation restart

    # get current restart file status
    io_stat = pytools.check_for_restart(conf)


    # no restart file; initialize simulation
    if io_stat["do_initialization"]:
        if sch.is_master: print("initializing simulation...")
        lap = 0

        rseed = MPI.COMM_WORLD.Get_rank()
        np.random.seed(rseed + rnd_seed_default)  # sync rnd generator seed for different mpi ranks

        # injecting plasma particles
        # NOTE: we still need to call injector to set the container types; density_profile ensure = 0
        if not(conf.use_injector): 
            prtcl_stat = pytools.pic.inject(grid, velocity_profile, density_profile, conf)
            if sch.is_example_worker: 
                print("injected:")
                print("     e- prtcls: {}".format(prtcl_stat[0]))
                print("     e+ prtcls: {}".format(prtcl_stat[1]))

        # inserting em grid
        insert_em_fields(grid, conf, do_initialization=True)

        # save a snapshot of the state to disk
        pytools.save_mpi_grid_to_disk(conf.outdir, 0, grid, conf)

    else:
        if sch.is_master:
            print("restarting simulation from lap {}...".format(io_stat["lap"]))

        # insert damping reference values
        insert_em_fields(grid, conf, do_initialization=False)

        # read restart files
        pyfld.read_grids(grid, io_stat["read_lap"], io_stat["read_dir"])
        pypic.read_particles(grid, io_stat["read_lap"], io_stat["read_dir"])

        # set particle types
        for tile in pytools.tiles_all(grid):
            for ispcs in range(conf.Nspecies):
                container = tile.get_container(ispcs)
                container.type = conf.prtcl_types[ispcs] # name container

        # step one step ahead
        lap = io_stat["lap"] + 1

    # --------------------------------------------------
    # static load balancing setup; communicate neighborhood info once

    if sch.is_master: print("load balancing grid..."); sys.stdout.flush()


    # update boundaries
    grid.analyze_boundaries()
    grid.send_tiles()
    grid.recv_tiles()
    MPI.COMM_WORLD.barrier()

    if sch.is_master: print("loading virtual tiles..."); sys.stdout.flush()

    # load virtual mpi halo tiles
    pytools.pic.load_virtual_tiles(grid, conf)


    # --------------------------------------------------
    # load physics solvers

    if sch.is_master: print("loading solvers..."); sys.stdout.flush()


    #sch.fldpropE = pyfld.FDTD2_pml()
    #sch.fldpropB = pyfld.FDTD2_pml()

    sch.fldpropE = pyfld.FDTD2()
    sch.fldpropB = pyfld.FDTD2()

    #sch.fldpropE = pyfld.FDTD4()
    #sch.fldpropB = pyfld.FDTD4()

    # enhance numerical speed of light slightly to suppress numerical Cherenkov instability
    sch.fldpropE.corr = conf.c_corr
    sch.fldpropB.corr = conf.c_corr

    # --------------------------------------------------
    # set PML absorption region into the pusher
 
    # setup perfectly matching layer; vacuum outside boundaries
    #for prop in [ sch.fldpropE, sch.fldpropB]:
    #    prop.cenx = conf.Lx//2
    #    prop.ceny = conf.Ly//2 
    #    prop.cenz = 0

    #    prop.radx = conf.Lx//2 # box x half length
    #    prop.rady = conf.Ly//2 # box y half length
    #    prop.radz = 1;        # box z half length

    #    prop.norm_abs = conf.cfl/3.0 # PML coefficient
    #    prop.rad_lim = 0.80 # PML sphere damping radius


    # --------------------------------------------------
    #sch.pusher = pypic.BorisPusher()
    #sch.pusher = pypic.VayPusher()
    #sch.pusher = pypic.HigueraCaryPusher()
    #sch.pusher  = pypic.rGCAPusher()
    sch.pusher  = pypic.PulsarPusher()

    #if conf.gammarad > 0:
    #    sch.pusher   = pypic.BorisDragPusher()
    #    sch.pusher.drag = conf.drag_amplitude
    #    sch.pusher.temp = 0.0

    sch.pusherx= pypic.PhotonPusher() # photon pusher

    # background field from external pusher components
    if conf.use_maxwell_split:
        sch.pusher.bx_ext = conf.bx_ext 
        sch.pusher.by_ext = conf.by_ext
        sch.pusher.bz_ext = conf.bz_ext

        sch.pusher.ex_ext = conf.ex_ext 
        sch.pusher.ey_ext = conf.ey_ext
        sch.pusher.ez_ext = conf.ez_ext


    sch.fintp = pypic.LinearInterpolator()
    #sch.fintp = pypic.QuadraticInterpolator() #2nd order quadratic
    #sch.fintp = pypic.CubicInterpolator() #3rd order cubic 3d
    #sch.fintp = pypic.QuarticInterpolator() #4th order quartic; 3d

    sch.currint = pypic.ZigZag()
    #sch.currint = pypic.ZigZag_2nd()
    #sch.currint = pypic.ZigZag_3rd()
    #sch.currint = pypic.ZigZag_4th()
    #sch.currint = pypic.Esikerpov_2nd() # 3d
    #sch.currint = pypic.Esikerpov_4th() # 3d

    sch.flt = pyfld.Binomial2(conf.NxMesh, conf.NyMesh, conf.NzMesh)


    # QED
    mc = pyqed.Pairing()
    mc.prob_norm = 1/(conf.N_time*conf.N_wgt*conf.N_qdt) # units of [per prtcl per time]

    # histogram edges of leaking photons
    mc.update_hist_lims(toolset.xxlims[0], toolset.xxlims[1], toolset.Nhist)

    a = pyrunko.qed.PhotAnn("ph", "ph")
    b = pyrunko.qed.PairAnn("e-", "e+")
    c = pyrunko.qed.PairAnn("e+", "e-")
    d = pyrunko.qed.Compton("ph", "e-")
    e = pyrunko.qed.Compton("ph", "e+")
    f = pyrunko.qed.Compton("e-", "ph")
    g = pyrunko.qed.Compton("e+", "ph")

    #mc.add_interaction(a) # ON                      # phot-ann
    #mc.add_interaction(b) # ON                      # pair-ann
    #mc.add_interaction(c) # off for double counting # pair-ann
    #mc.add_interaction(d) # ON
    #mc.add_interaction(e) # ON
    #mc.add_interaction(f) # off for double counting
    #mc.add_interaction(g) # off for double counting


    #--------------------------------------------------
    mc.prob_norm_onebody = conf.N_onebody/conf.N_qdt # units of [TODO]

    a0 = pyrunko.qed.Synchrotron("e-")
    a1 = pyrunko.qed.Synchrotron("e+")
    a0.do_accumulate = True
    a1.do_accumulate = True

    b  = pyrunko.qed.MultiPhotAnn("ph")

    # set critical magnetic field 
    for intr in [a0, a1, b]:
        intr.B_QED = conf.binit/conf.B_QED

    mc.add_interaction(a0) 
    mc.add_interaction(a1) 
    mc.add_interaction(b ) # 



    # --------------------------------------------------
    # I/O objects
    if sch.is_master: print("loading IO objects..."); sys.stdout.flush()

    # quick field snapshots
    #fld_writer = pyfld.MasterFieldsWriter(
    fld_writer = pyfld.FieldsWriter(
        conf.outdir,
        conf.Nx,
        conf.NxMesh,
        conf.Ny,
        conf.NyMesh,
        conf.Nz,
        conf.NzMesh,
        conf.stride,
    )


    # test particles; only works with no injector (id's are messed up when coming out from injector)
    prtcl_writers = []

    #for ispc in [0]: #electrons
    #for ispc in [0, 1, 2]: #electrons & positrons
    #    mpi_comm_size = grid.size() if not(conf.mpi_task_mode) else grid.size() - 1
    #    n_local_tiles = int(conf.Nx*conf.Ny*conf.Nz/mpi_comm_size)
    #    #print("average n_local_tiles:", n_local_tiles, " / ", mpi_comm_size)

    #    prtcl_writer = pypic.TestPrtclWriter(
    #            conf.outdir,
    #            conf.Nx, conf.NxMesh, conf.Ny, conf.NyMesh, conf.Nz, conf.NzMesh,
    #            conf.ppc,
    #            n_local_tiles, #len(grid.get_local_tiles()),
    #            conf.n_test_prtcls,)
    #    prtcl_writer.ispc = ispc
    #    prtcl_writers.append(prtcl_writer)


    #mom_writer = pypic.MasterPicMomentsWriter(
    mom_writer = pypic.PicMomentsWriter(
        conf.outdir,
        conf.Nx,
        conf.NxMesh,
        conf.Ny,
        conf.NyMesh,
        conf.Nz,
        conf.NzMesh,
        conf.stride_mom,
    )

    # 3D box peripherals
    if conf.threeD:
        st = 1 # stride
        slice_xy_writer = pyfld.FieldSliceWriter( conf.outdir, 
                conf.Nx, conf.NxMesh, conf.Ny, conf.NyMesh, conf.Nz, conf.NzMesh, st, 0, 1)
        slice_xz_writer = pyfld.FieldSliceWriter( conf.outdir, 
                conf.Nx, conf.NxMesh, conf.Ny, conf.NyMesh, conf.Nz, conf.NzMesh, st, 1, 1)
        slice_yz_writer = pyfld.FieldSliceWriter( conf.outdir, 
                conf.Nx, conf.NxMesh, conf.Ny, conf.NyMesh, conf.Nz, conf.NzMesh, st, 2, 1)

        slice_xy_writer.ind = int(0.8*conf.Lz) # bottom slice; z = 8 slice through atmosphere
        slice_xz_writer.ind = int(0.5*conf.Ly) # mid slice 
        slice_yz_writer.ind = int(0.5*conf.Lx) # mid slice 


    # --------------------------------------------------
    #star = pyfld.Conductor()
    star = pypic.Star()
    for obj in [star, sch.pusher]:
        obj.radius    = conf.rad_star
        obj.radius_pc = conf.rad_pcap
        obj.period    = conf.period_star # NOTE should be period_star for normal runs 

        if conf.twoD:
            obj.cenx   = conf.Lx//2 + 0.5
            obj.ceny   = -conf.rad_star + conf.rad_curv_shift
            obj.cenz   = 0 
        elif conf.threeD:
            obj.cenx   = conf.Lx//2 #+ 0.5
            obj.ceny   = conf.Ly//2 #+ 0.5
            obj.cenz   = -conf.rad_star + conf.rad_curv_shift


    #sch.pusher.grav_const = 0.0 # gravitational constant
    height_atms = conf.height_atms #0.005*conf.rad_star
    sch.pusher.grav_const = conf.delgam**2/(2*height_atms)

    star.B0       = conf.b_dipole_norm    # TODO normalization here?
    star.chi_om   = np.deg2rad(conf.chi)  # rotation axis inclination
    star.phase_mu = 0.0
    star.phase_om = 0.0
    star.delta    = 0.5*conf.height_atms #4.0 #1  # in units of cells; radial smoothing function sharpness; 2x delta = about tanh limit

    star.Nx = conf.Lx
    star.Ny = conf.Ly
    star.Nz = conf.Lz
    
    star.delta_pc  = 2 # transverse polar cap smoothing (g(x) function)

    star.temp_pairs = conf.delgam
    star.temp_phots = conf.delgam_x
    star.ninj_pairs = conf.ninj_pairs #0 #0.05
    star.ninj_phots = conf.ninj_phots #0.0
    star.ninj_min_pairs = conf.ninj_min_pairs #0.1
    star.ninj_min_phots = conf.ninj_min_phots #0.0

    sch.lwall = star # add to scheduler


    # induce initial magnetic and electric field from the star
    for tile in pytools.tiles_all(grid):
        star.insert_em(tile)

    # --------------------------------------------------
    # --------------------------------------------------
    # --------------------------------------------------
    # end of initialization
    timer.stop("init")
    timer.stats("init")

    if sch.is_master: print('init: starting simulation...'); sys.stdout.flush()

    ##################################################
    # simulation time step loop

    # simulation loop
    time = lap * (conf.cfl / conf.c_omp)
    for lap in range(lap, conf.Nt + 1):

        # ramp up plate smoothly
        #ramp_up_laps = 1.0*conf.rad_pcap/conf.cfl # duration of the ramp up in polar cap light crossing times
        #pc_freq = min(max(1,lap)/ramp_up_laps, 1.0)*(1/conf.period_star) # polar cap rotation frequency
        #pc_freq = 1e-5/conf.period_star # polar cap rotation frequency
        #star.period = 1/pc_freq

        # rotate star's Omega vector with the period
        sch.lwall.phase_om += conf.Om_star

        # --------------------------------------------------
        # QED interaction loop
        if conf.qed_mode and lap % conf.qed_step == 0:

            #timer.start_comp("ph_inj")
            #if conf.zeta_xinj1 > 0 and conf.lum_ph1 > 0:
            #    for tile in pytools.tiles_local(grid):
            #        mc.inject_photons(tile, conf.kTbb1, conf.wph_inj1, conf.Nph_inj1 )

            #if conf.zeta_xinj2 > 0 and conf.lum_ph2 > 0:
            #    for tile in pytools.tiles_local(grid):
            #        mc.inject_photons(tile, conf.kTbb2, conf.wph_inj2, conf.Nph_inj2 )
            #timer.stop_comp("ph_inj")

            #timer.start_comp("ep_inj")
            #for tile in pytools.tiles_local(grid):
            #    if conf.zeta_pinj > 0 and conf.lum_ep > 0:
            #        mc.inject_plaw_pairs(tile, conf.pslope, conf.pmin, conf.pmax, conf.wep_inj, conf.Nep_inj )
            #timer.stop_comp("ep_inj")

            #--------------------------------------------------
            #timer.start_comp("qed")
            #for tile in pytools.tiles_local(grid):
            #    i,j,k = pytools.get_index(tile, conf)
            #    mc.solve_twobody(tile)
            #timer.stop_comp("qed")

            #--------------------------------------------------
            timer.start_comp("ph_esc")

            # ver2: local tau
            tau_tile_min = 1e7
            tau_tile_max = 0
            for tile in pytools.tiles_local(grid): 
                mc.tau_global = 0.0 # clear tau measure 
                mc.comp_tau(tile, conf.N_wgt) # sum over tiles

                # TODO
                #mc.leak_photons(tile, conf.t_c/conf.dt/conf.N_qdt, conf.tau_ext)  # apply photon escape

                tau_tile_min = mc.tau_global if mc.tau_global < tau_tile_min else tau_tile_min
                tau_tile_max = mc.tau_global if mc.tau_global > tau_tile_max else tau_tile_max

            # store values
            toolset.tau_tile_min = tau_tile_min
            toolset.tau_tile_max = tau_tile_max


            # ver1; global tau
            #mc.tau_global = 0.0 # clear tau measure 
            #for tile in pytools.tiles_local(grid): 
            #    mc.comp_tau(tile, conf.N_tau) # sum over tiles # TODO
            #
            # MPI sum over ranks TODO no tau reduction
            #tau_global = mc.tau_global
            #tau_global = MPI.COMM_WORLD.allreduce(tau_global, op=MPI.SUM)
            #mc.tau_global = tau_global

            #for tile in pytools.tiles_local(grid):
            #    mc.leak_photons(tile, conf.t_c/conf.dt, conf.tau_ext)  # apply photon escape

            timer.stop_comp("ph_esc")

            # TODO testing second particle comm after QED

            # local and global particle exchange 
            #sch.operate( dict(name='check_outg_prtcls',     solver='tile',  method='check_outgoing_particles',     nhood='local', ) )
            #sch.operate( dict(name='pack_outg_prtcls',      solver='tile',  method='pack_outgoing_particles',      nhood='boundary', ) )

            #sch.operate( dict(name='mpi_prtcls',            solver='mpi',   method='p1',                           nhood='all', ) )
            #sch.operate( dict(name='mpi_prtcls',            solver='mpi',   method='p2',                           nhood='all', ) )

            #sch.operate( dict(name='unpack_vir_prtcls',     solver='tile',  method='unpack_incoming_particles',    nhood='virtual', ) )
            #sch.operate( dict(name='check_outg_vir_prtcls', solver='tile',  method='check_outgoing_particles',     nhood='virtual', ) )
            #sch.operate( dict(name='get_inc_prtcls',        solver='tile',  method='get_incoming_particles',       nhood='local', args=[grid,]) )

            #sch.operate( dict(name='del_trnsfrd_prtcls',    solver='tile',  method='delete_transferred_particles', nhood='local', ) )
            #sch.operate( dict(name='del_vir_prtcls',        solver='tile',  method='delete_all_particles',         nhood='virtual', ) )


        # --------------------------------------------------
        # comm E and B
        sch.operate( dict(name='mpi_b0', solver='mpi', method='b', ) )
        sch.operate( dict(name='mpi_e0', solver='mpi', method='e', ) )
        sch.operate( dict(name='upd_bc', solver='tile',method='update_boundaries', args=[grid,[1,2] ], nhood='local', ) )

        # --------------------------------------------------
        # push B half
        sch.operate( dict(name='push_half_b1', solver='fldpropB', method='push_half_b', nhood='local',) )
        sch.operate( dict(name='wall_bc_b',    solver='lwall',    method='update_b',    nhood='local',) )

        # comm B
        sch.operate( dict(name='mpi_b1',    solver='mpi', method='b',                 ) )
        sch.operate( dict(name='upd_bc ',   solver='tile',method='update_boundaries',args=[grid, [2,] ], nhood='local',) )


        # interpolate fields and push particles in x and u
        #sch.operate( dict(name='interp_em', solver='fintp',  method='solve', nhood='local', ) )
        #sch.operate( dict(name='push',      solver='pusher', method='solve', nhood='local', ) )

        #--------------------------------------------------
        # single-body QED interactions
        if conf.qed_mode and lap % conf.qed_step == 0:
            timer.start_comp("qed_onebody")
            for tile in pytools.tiles_local(grid):
                #print("calling onebody")
                mc.solve_onebody(tile)
            timer.stop_comp("qed_onebody")

        # --------------------------------------------------
        # move particles (only locals tiles)

        # NOTE need to recalculate the interpolation step since new particles dont have up-to-date Bpart and Epart; 
        #      this is clearly a design error in mc.solve_onebody...
        # DONE new PulsarPusher calculates the fields inside the prtcl loop so this is fixed
        #      for every other pusher, need to uncomment this
        sch.operate( dict(name='interp_em', solver='fintp',  method='solve', nhood='local', ) )

        sch.operate( dict(name='push',      solver='pusher', method='solve', nhood='local', args=[0]) ) # e^-
        sch.operate( dict(name='push',      solver='pusher', method='solve', nhood='local', args=[1]) ) # e^+
        sch.operate( dict(name='push',      solver='pusherx',method='solve', nhood='local', args=[2]) ) # x

        # clear currents; need to call this before wall operations since they can deposit currents too 
        sch.operate( dict(name='clear_cur', solver='tile',   method='clear_current', nhood='all', ) )

        # apply moving/reflecting/injecting walls
        # TODO
        #if lap*conf.cfl > 1.0*conf.rad_pcap: # apply after a fraction of the disk light crossing time 
        sch.operate( dict(name='star',     solver='lwall', method='solve', nhood='local', ) )


        # --------------------------------------------------
        # advance half B 
        sch.operate( dict(name='push_half_b2', solver='fldpropB', method='push_half_b', nhood='local', ) )
        sch.operate( dict(name='wall_bc_b',    solver='lwall',    method='update_b',    nhood='local',) ) # dynamic B

        # comm B
        sch.operate( dict(name='mpi_b2', solver='mpi', method='b',                 ) )
        sch.operate( dict(name='upd_bc', solver='tile',method='update_boundaries', args=[grid, [2,] ], nhood='local',) )

        # --------------------------------------------------
        # push E
        sch.operate( dict(name='push_e',    solver='fldpropE', method='push_e',  nhood='local', ) )
        sch.operate( dict(name='wall_bc_e', solver='lwall',    method='update_e',nhood='local', ) )


        # TODO current deposit + MPI was here

        # --------------------------------------------------
        # particle communication (only local/boundary tiles)

        # local and global particle exchange 
        sch.operate( dict(name='check_outg_prtcls',     solver='tile',  method='check_outgoing_particles',     nhood='local', ) )
        sch.operate( dict(name='pack_outg_prtcls',      solver='tile',  method='pack_outgoing_particles',      nhood='boundary', ) )

        sch.operate( dict(name='mpi_prtcls',            solver='mpi',   method='p1',                           nhood='all', ) )
        sch.operate( dict(name='mpi_prtcls',            solver='mpi',   method='p2',                           nhood='all', ) )

        sch.operate( dict(name='unpack_vir_prtcls',     solver='tile',  method='unpack_incoming_particles',    nhood='virtual', ) )
        sch.operate( dict(name='check_outg_vir_prtcls', solver='tile',  method='check_outgoing_particles',     nhood='virtual', ) )
        sch.operate( dict(name='get_inc_prtcls',        solver='tile',  method='get_incoming_particles',       nhood='local', args=[grid,]) )

        sch.operate( dict(name='del_trnsfrd_prtcls',    solver='tile',  method='delete_transferred_particles', nhood='local', ) )
        sch.operate( dict(name='del_vir_prtcls',        solver='tile',  method='delete_all_particles',         nhood='virtual', ) )


        # --------------------------------------------------
        # current calculation; charge conserving current deposition
        # clear virtual current arrays for boundary addition after mpi, send currents, and exchange between tiles
        sch.operate( dict(name='comp_curr',     solver='currint', method='solve',             nhood='local', ) )
        sch.operate( dict(name='clear_vir_cur', solver='tile',    method='clear_current',     nhood='virtual', ) )
        sch.operate( dict(name='mpi_cur',       solver='mpi',     method='j',                 nhood='all', ) )
        sch.operate( dict(name='cur_exchange',  solver='tile',    method='exchange_currents', nhood='local', args=[grid,], ) )


        # --------------------------------------------------
        # filter
        for fj in range(conf.npasses):


            # flt uses halo=2 padding so only every 3rd (0,1,2) pass needs update
            if fj % 2 == 0:
                sch.operate( dict(name='mpi_cur_flt', solver='mpi', method='j', ) )
                sch.operate( dict(name='upd_bc',      solver='tile',method='update_boundaries',args=[grid, [0,] ], nhood='local', ) )
                MPI.COMM_WORLD.barrier()
            sch.operate( dict(name='filter', solver='flt', method='solve', nhood='local', ) )


        # --------------------------------------------------
        # add antenna
        #sch.antenna.update_rnd_phases()
        #antenna.get_brms(grid)
        #sch.operate( dict(name='add_antenna', solver='antenna', method='add_ext_cur', nhood='local', ) )

        # --------------------------------------------------
        # add current to E
        sch.operate( dict(name='add_cur',   solver='tile',  method='deposit_current',       nhood='local', ) )
        sch.operate( dict(name='wall_bc_e', solver='lwall', method='update_e',              nhood='local', ) )


        ##################################################
        # data reduction and I/O

        timer.lap("step")
        if lap % conf.interval == 0:
            if sch.is_master:
                print("--------------------------------------------------")
                print("------ lap: {} / t: {}".format(lap, time))

            timer.stats("step")
            timer.comp_stats()
            timer.purge_comps()

            # internal pairing timers
            if sch.is_example_worker: 
                mc.timer_stats()
            mc.timer_clear()


            # io/analyze (independent)
            timer.start("io")
            # shrink particle arrays
            for tile in pytools.tiles_all(grid):
                tile.shrink_to_fit_all_particles()

            # barrier for quick writers
            MPI.COMM_WORLD.barrier()

            # shallow IO
            # NOTE: do moms before other IOs to keep rho up-to-date
            mom_writer.write(grid, lap)  # pic distribution moments; 
            fld_writer.write(grid, lap)  # quick field snapshots

            for pw in prtcl_writers:
                pw.write(grid, lap)  # test particles
            
            #pytools.save_mpi_grid_to_disk(conf.outdir, lap, grid, conf)

            #box peripheries 
            if conf.threeD:
                slice_xy_writer.write(grid, lap)
                slice_xz_writer.write(grid, lap)
                slice_yz_writer.write(grid, lap)

            #--------------------------------------------------
            # terminal plot 
            if sch.is_master:

                #epar
                #d = fld_writer.get_slice(1) # use ex xy-slice as the image
                #d_norm = 0.1*conf.binit #

                # rho
                #d = fld_writer.get_slice(9) # rho
                #d_norm = conf.ppc
                #tplt.plot(np.abs(d)/d_norm)

                tplt.col_mode = False

                if conf.twoD:
                    tplt.plot_panels( (2,3),
                    dict(axs=(0,0), data=fld_writer.get_slice( 0)/conf.e_norm ,   name='ex', cmap='RdBu'   ,vmin=-1, vmax=1),
                    dict(axs=(0,1), data=fld_writer.get_slice( 1)/conf.e_norm ,   name='ey', cmap='RdBu'   ,vmin=-1, vmax=1),
                    dict(axs=(0,2), data=fld_writer.get_slice( 9)/conf.p_norm   , name='ne', cmap='viridis',vmin= 0, vmax=4),
                    dict(axs=(1,0), data=fld_writer.get_slice( 3)/conf.b_norm ,   name='bx', cmap='RdBu'   ,vmin=-1, vmax=1),
                    dict(axs=(1,1), data=fld_writer.get_slice( 4)/conf.b_norm ,   name='by', cmap='RdBu'   ,vmin=-1, vmax=1),
                    dict(axs=(1,2), data=mom_writer.get_slice(14)/conf.x_norm   , name='ph', cmap='viridis',vmin= 0, vmax=4),
                    )
                if conf.threeD:
                    tplt.plot_panels( (2,3),
                    dict(axs=(0,0), data=slice_xz_writer.get_slice( 0)/conf.e_norm ,   name='ex (xz)', cmap='RdBu'   ,vmin=-1, vmax=1),
                    dict(axs=(0,1), data=slice_xz_writer.get_slice( 2)/conf.e_norm ,   name='ez (xz)', cmap='RdBu'   ,vmin=-1, vmax=1),
                    dict(axs=(0,2), data=slice_xz_writer.get_slice( 9)/conf.p_norm   , name='ne (xz)', cmap='viridis',vmin= 0, vmax=1),
                    dict(axs=(1,0), data=slice_xz_writer.get_slice( 3)/conf.b_norm ,   name='bx (xz)', cmap='RdBu'   ,vmin=-1, vmax=1),
                    dict(axs=(1,1), data=slice_xy_writer.get_slice( 2)/conf.e_norm   , name='ez (top)',cmap='RdBu',   vmin=-1, vmax=1),
                    dict(axs=(1,2), data=slice_xy_writer.get_slice( 9)/conf.p_norm ,   name='ne (top)',cmap='viridis',vmin=-0, vmax=1),
                    )

            #--------------------------------------------------

            #print statistics
            if sch.is_master:
                print('simulation time    {:7d} ({:7.1f} omp)   {:5.1f}%'.format( int(lap), time, 100.0*lap/conf.Nt))
                #print('sim time {:7d} ({:7.1f} omp) ({:5.1f} H/c) ({:7.2f} l0/c) {:5.2f}%'.format(int(lap), 
                #                                                                         time, 
                #                                                                         lap*conf.dt/conf.t_c,
                #                                                                         sch.antenna.tcur,
                #                                                                         100.0*lap/conf.Nt))


            #--------------------------------------------------
            # deep IO
            if conf.full_interval > 0 and (lap % conf.full_interval == 0) and (lap > 0):
                pyfld.write_grids(grid, lap, conf.outdir + "/full_output/")
                pypic.write_particles(grid, lap, conf.outdir + "/full_output/")


            # restart IO (overwrites)
            if (lap % conf.restart == 0) and (lap > 0):

                # flip between two sets of files
                io_stat["deep_io_switch"] = 1 if io_stat["deep_io_switch"] == 0 else 0

                pyfld.write_grids(
                    grid, io_stat["deep_io_switch"] + io_stat['restart_num'],
                    conf.outdir + "/restart/"
                )

                pypic.write_particles(
                    grid, io_stat["deep_io_switch"] + io_stat['restart_num'],
                    conf.outdir + "/restart/"
                )

                # if successful adjust info file
                MPI.COMM_WORLD.barrier()  # sync everybody in case of failure before write
                if grid.rank() == 0:
                    with open(conf.outdir + "/restart/laps.txt", "a") as lapfile:
                        lapfile.write("{},{}\n".format(
                            lap, 
                            io_stat["deep_io_switch"]+io_stat['restart_num']))

            MPI.COMM_WORLD.barrier() # extra barrier to synch everybody after IOs

            timer.stop("io")
            timer.stats("io")


        ###################################################
        # QED plotting
        if conf.qed_mode and lap % conf.N_qdt == 0:
            toolset.update_stats(grid, lap, conf)


        if lap % conf.plot_interval == 0 and live_plot:

            timer.start("io2")

            if conf.qed_mode:
                toolset.update_hists(grid, lap, conf)
                toolset.update_esc_stats(mc, lap, conf) # needs to be called w/ plot_interval
                toolset.update_esc_hists(mc, lap, conf) # needs to be called w/ plot_interval
                mc.clear_hist() # remember to clear histogram after reading

            if conf.qed_mode and sch.is_master:

                col = cmap(norm(lap)) # NOTE overwrite color
                lw = 0.8
                ls = 'solid'

                axs[0,0].plot(toolset.lnxs, toolset.h1_enes['ph'], drawstyle='steps-pre', color=col, alpha=1.0, lw = lw, linestyle=ls,)

                #if lap > 1: # ignore first time step
                #    axs[0,2].plot(toolset.lnxs, toolset.h1_enes['esc'], drawstyle='steps-pre', color=col, alpha=1.0, lw = lw, linestyle=ls,)

                zs = toolset.zs
                axs[0,1].plot(toolset.lnzs, toolset.h1_enes['e-']*zs, drawstyle='steps-pre', color=col, alpha=1.0, lw = lw, linestyle=ls,)
                axs[0,2].plot(toolset.lnzs, toolset.h1_enes['e+']*zs, drawstyle='steps-pre', color=col, alpha=1.0, lw = lw, linestyle=ls,)

                # LP weight spectrum
                axs[1,0].plot(toolset.lnxs, toolset.h1_ws['ph'],    drawstyle='steps-pre', color=col, alpha=1.0, lw = lw, linestyle=ls,)
                axs[1,1].plot(toolset.lnzs, toolset.h1_ws['e-']*zs, drawstyle='steps-pre', color=col, alpha=1.0, lw = lw, linestyle=ls,)
                axs[1,1].plot(toolset.lnzs, toolset.h1_ws['e+']*zs, drawstyle='steps-pre', color=col, alpha=1.0, lw = lw, linestyle=ls,)

                #--------------------------------------------------
                # LP w spectrum
                axs[2,0].plot(toolset.lnws, toolset.h1_nums['ph'], drawstyle='steps-pre', color=col, alpha=1.0, lw = lw, linestyle=ls,)
                axs[2,1].plot(toolset.lnws, toolset.h1_nums['e-'], drawstyle='steps-pre', color=col, alpha=1.0, lw = lw, linestyle=ls,)
                axs[2,1].plot(toolset.lnws, toolset.h1_nums['e+'], drawstyle='steps-pre', color=col, alpha=1.0, lw = lw, linestyle=ls,)

                im30.set_data(np.log10(toolset.h2_nums.T)) # 2d weight - ene histogram
                im31.set_data(np.log10(toolset.h2_enes['ph'].T)) # 2d spatial histogram
                im32.set_data(np.log10(toolset.h2_enes['e-'].T)) # 2d spatial histogram
                im33.set_data(np.log10(toolset.h2_enes['e+'].T)) # 2d spatial histogram


                #--------------------------------------------------
                # line plots
                laps = toolset.storage.data['lap']
                laps_sparse = toolset.storage.data['lap_sparse']
                ts = (conf.dt/conf.t_c)*np.array(laps)

                # optical depth
                axs[0,3].plot(laps[-1],        toolset.storage.data['tau'][-1],      color='C0', marker='.')
                axs[0,3].plot(laps_sparse[-1], toolset.storage.data['tau_meas'][-1], color='C1', marker='.')
                axs[0,3].plot(laps_sparse[-1], toolset.storage.data['tau_min'][-1],  color='C1', marker='^')
                axs[0,3].plot(laps_sparse[-1], toolset.storage.data['tau_max'][-1],  color='C1', marker='v')

                # energy
                axs[1,3].plot(ts[-1], toolset.storage.data['ene_e-'][-1], color='C0', marker='.')
                axs[1,3].plot(ts[-1], toolset.storage.data['ene_e+'][-1], color='C1', marker='.')
                axs[1,3].plot(ts[-1], toolset.storage.data['ene_ph'][-1], color='C2', marker='.')

                # energy
                Ux = np.array(toolset.storage.data['ene_ph'])
                Ue = np.array(toolset.storage.data['ene_e-'])
                Up = np.array(toolset.storage.data['ene_e+'])

                axs[1,2].plot(ts[-1], Ux[-1]/(Ue[-1]+Up[-1]), color='C0', marker='.')

                # LP number
                #print('num_e', toolset.storage.data['num_e-'])
                axs[2,3].plot(laps[-1], toolset.storage.data['lp_num_e-'][-1], color='C0', marker='.')
                axs[2,3].plot(laps[-1], toolset.storage.data['lp_num_e+'][-1], color='C1', marker='.')
                axs[2,3].plot(laps[-1], toolset.storage.data['lp_num_ph'][-1], color='C2', marker='.')

                axs[2,3].plot(laps_sparse[-1], toolset.storage.data['lp_num_esc'][-1], color='C3', marker='.')

                lum_in = conf.lum_ep + conf.lum_ph1 + conf.lum_ph2
                #axs[3,3].axhline(y=lum_in, lw=0.5)
                #axs[3,3].plot(laps_sparse[-1], toolset.storage.data['ene_inj_ep'][-1], color='C0', marker='.')
                #axs[3,3].plot(laps_sparse[-1], toolset.storage.data['ene_inj_ph'][-1], color='C2', marker='.')
                #axs[3,3].plot(laps_sparse[-1], toolset.storage.data['ene_esc'][-1],    color='C3', marker='.')

                axs[4,3].plot(laps_sparse[-1], toolset.storage.data['lum_rat'][-1], color='C0', marker='.')
                axs[4,3].axhline(y=1.0, lw=0.5)

                axs[2,2].plot(laps[-1], toolset.storage.data['num_e-'][-1], color='C0', marker='.')
                axs[2,2].plot(laps[-1], toolset.storage.data['num_e+'][-1], color='C1', marker='.')
                axs[2,2].plot(laps[-1], toolset.storage.data['num_ph'][-1], color='C2', marker='.')
                axs[2,2].plot(laps_sparse[-1], toolset.storage.data['num_esc'][-1], color='C3', marker='.')


                # QED analysis
                toolset.save(lap, conf)

                #--------------------------------------------------    
                # close and save fig
                fig.subplots_adjust(left=axleft, bottom=axbottom, right=axright, top=axtop)
                fname = conf.outdir + '/qed_analysis.pdf'
                plt.savefig(fname)


                #--------------------------------------------------    
                # extra prints
                #print("--------------------------------------------------")
                #print(' tau:  {:6.3f} min/max ({:6.3f}, {:6.3f}'.format(
                #    toolset.storage.data['tau'][-1],
                #    toolset.storage.data['tau_min'][-1],
                #    toolset.storage.data['tau_max'][-1])
                #      )
                #print(' lin/lout: {:5.1f} / {:5.1f} = {:6.3f}'.format(
                #    toolset.storage.data['ene_inj_ep'][-1] + 
                #    toolset.storage.data['ene_inj_ph'][-1],
                #    toolset.storage.data['ene_esc'][-1],
                #    toolset.storage.data['lum_rat'][-1],
                #      ))
                #print(' N/N_0: ph {:6.1f} e- {:6.1f} e+ {:6.1f}'.format(
                #    toolset.storage.data['lp_num_ph'][-1],
                #    toolset.storage.data['lp_num_e-'][-1],
                #    toolset.storage.data['lp_num_e+'][-1]
                #      ))

            # erase histogram and start a new monitoring cycle
            if lap % conf.plot_interval == 0:
                mc.clear_hist()

            timer.stop("io2")
            timer.stats("io2")
            #--------------------------------------------------


            timer.start("step")  # refresh lap counter (avoids IO profiling)

            sys.stdout.flush()

        # next step
        time += conf.cfl / conf.c_omp
        #MPI.COMM_WORLD.barrier() # extra barrier to synch everybody
    # end of loop

    # --------------------------------------------------
    # end of simulation

    timer.stop("total")
    timer.stats("total")

    if conf.qed_mode and sch.is_master and live_plot:
        print('-------finishing plotting---------')

        laps = toolset.storage.data['lap']
        laps_sparse = toolset.storage.data['lap_sparse']

        ts = (conf.dt/conf.t_c)*np.array(laps)

        # optical depth
        axs[0,3].plot(laps, toolset.storage.data['tau'])
        axs[0,3].plot(laps_sparse, toolset.storage.data['tau_meas'], color='C1', ls='solid')
        axs[0,3].plot(laps_sparse, toolset.storage.data['tau_min'],  color='C1', ls='solid')
        axs[0,3].plot(laps_sparse, toolset.storage.data['tau_max'],  color='C1', ls='solid')

        # energy
        axs[1,3].plot(ts, toolset.storage.data['ene_e-'], color='C0')
        axs[1,3].plot(ts, toolset.storage.data['ene_e+'], color='C1')
        axs[1,3].plot(ts, toolset.storage.data['ene_ph'], color='C2')

        # energy
        Ux = np.array(toolset.storage.data['ene_ph'])
        Ue = np.array(toolset.storage.data['ene_e-'])
        Up = np.array(toolset.storage.data['ene_e+'])

        axs[1,2].plot(ts, Ux/(Ue+Up), color='C0')

        # LP number
        #print('num_e', toolset.storage.data['num_e-'])
        axs[2,3].plot(laps, toolset.storage.data['lp_num_e-'], color='C0')
        axs[2,3].plot(laps, toolset.storage.data['lp_num_e+'], color='C1')
        axs[2,3].plot(laps, toolset.storage.data['lp_num_ph'], color='C2')

        axs[2,3].plot(laps_sparse, toolset.storage.data['lp_num_esc'], color='C3', ls='dashed')


        lum_in = conf.lum_ep + conf.lum_ph1 + conf.lum_ph2
        #axs[3,3].axhline(y=lum_in, lw=0.5)
        #axs[3,3].plot(laps_sparse, toolset.storage.data['ene_inj_ep'], color='C0', ls='solid')
        #axs[3,3].plot(laps_sparse, toolset.storage.data['ene_inj_ph'], color='C2', ls='solid')
        #axs[3,3].plot(laps_sparse, toolset.storage.data['ene_esc'],    color='C3', ls='dashed')

        axs[4,3].plot(laps_sparse, toolset.storage.data['lum_rat'], color='C0', ls='solid')
        axs[4,3].axhline(y=1.0, lw=0.5)


        axs[2,2].plot(laps, toolset.storage.data['num_e-'], color='C0')
        axs[2,2].plot(laps, toolset.storage.data['num_e+'], color='C1')
        axs[2,2].plot(laps, toolset.storage.data['num_ph'], color='C2')
        axs[2,2].plot(laps_sparse, toolset.storage.data['num_esc'], color='C3', ls='dashed')

        # plot the final spectrum

        col = 'C0'
        axs[4,0].plot(toolset.lnxs, toolset.h1_enes['ph'], drawstyle='steps-pre', color=col, alpha=1.0, lw = lw, linestyle=ls,)
        axs[4,2].plot(toolset.lnxs, toolset.h1_enes['esc'],drawstyle='steps-pre', color=col, alpha=1.0, lw = lw, linestyle=ls,)
        axs[4,1].plot(toolset.lnzs, toolset.h1_enes['e-'], drawstyle='steps-pre', color=col, alpha=1.0, lw = lw, linestyle=ls,)
        axs[4,1].plot(toolset.lnzs, toolset.h1_enes['e+'], drawstyle='steps-pre', color=col, alpha=1.0, lw = lw, linestyle=ls,)

        #--------------------------------------------------    
        # close and save fig
        if True:
            axleft    = 0.10
            axbottom  = 0.13
            axright   = 0.97
            axtop     = 0.97

            pos1 = axs[0,0].get_position()
            axwidth  = axright - axleft
            axheight = (axtop - axbottom)*0.02
            axpad = 0.01

            fig.subplots_adjust(left=axleft, bottom=axbottom, right=axright, top=axtop)

        fname = conf.outdir + '/qed_analysis.pdf'
        plt.savefig(fname)


