import numpy as np
import h5py as h5

import pytools  # runko python tools
import sys, os

from numpy import sqrt, pi 
from mpi4py import MPI
from tools import Storage, plot_connect_points


class QEDToolset:

    def __init__(self, conf):

        self.nx = conf.Nx*conf.NxMesh
        self.ny = conf.Ny*conf.NyMesh
        self.nz = conf.Nz*conf.NzMesh

        if conf.oneD:
            self.nh = self.nx
        elif conf.twoD:
            self.nh = self.ny
        elif conf.threeD:
            self.nh = self.nz

        # limits for various histograms
        self.xxlims = (-2, 7) # photon xlim
        self.xylims = (1e1, 1e5) #(1e-2, 1e6) # photon ylim (1e-4, 1e2)

        self.pxlims = (-2, 10) # pair xlim
        self.pylims = (1e-1, 1e2) # pair ylim

        self.wxlims = (-1, 7) # pair wlim
        self.hhlims = (0, self.nh) # spatial box height limits

        self.Nhist = 256 #128

        self.zs = np.logspace(self.pxlims[0], self.pxlims[1], self.Nhist, endpoint=False) # momenta
        self.xs = np.logspace(self.xxlims[0], self.xxlims[1], self.Nhist, endpoint=False) # energy
        self.ws = np.logspace(self.wxlims[0], self.wxlims[1], self.Nhist, endpoint=False) # weight
        self.hs = np.linspace(self.hhlims[0], self.hhlims[1], self.Nhist, endpoint=False) # height

        # histogram bins have +1 elemenets
        self.zs_bin = np.append(self.zs, 10.0**self.pxlims[1] ) # momenta
        self.xs_bin = np.append(self.xs, 10.0**self.xxlims[1] ) # energy
        self.ws_bin = np.append(self.ws, 10.0**self.wxlims[1] ) # weight
        self.hs_bin = np.append(self.hs,       self.hhlims[1] ) # height

        self.lnzs = np.log10(self.zs)
        self.lnxs = np.log10(self.xs)
        self.lnws = np.log10(self.ws)

        self.dlnz = np.diff(np.log10(self.zs_bin))
        self.dlnx = np.diff(np.log10(self.xs_bin))
        self.dlnw = np.diff(np.log10(self.ws_bin))

        self.tau_tile_min = 1e7
        self.tau_tile_max = 0

        self.types = ["e-", "e+", "ph"]

        # normalizations
        #self.N_ene = conf.N_w*conf.t_c/conf.dt # units of luminosity compactness
        #self.N_ene *= 1.0/(conf.Nx*conf.Ny*conf.Nz)

        self.N_wgt   = conf.N_wgt
        self.N_box   = conf.N_box
        self.N_time  = conf.N_time

        #self.N_tau = conf.N_tau #(conf.Nx*conf.Ny*conf.Nz)**2 # mean tau in tiles
        self.N_npc = conf.npc_ref*conf.N_box # total number of ref LPs in sim box


        # mark root rank
        self.root = True if MPI.COMM_WORLD.Get_rank() == 0 else False

        # data storage
        if self.root: 
            self.storage = Storage()


    def update_stats(self, grid, lap, conf):

        lp_nums = {}
        nums = {}
        enes = {}
        taus = {}

        for t in self.types:
            nums[t] = 0
            enes[t] = 0
            taus[t] = 0
            lp_nums[t] = 0



        for tile in pytools.tiles_local(grid):
            i,j,k = pytools.get_index(tile, conf)
            #print('tile', i,j,k)
            for ispcs in range(conf.Nspecies):
                container = tile.get_container(ispcs)
                t = container.type

                # LP number
                lp_nums[t] += container.size()/self.N_npc

                # energy
                ux = np.array( container.vel(0) )
                uy = np.array( container.vel(1) )
                uz = np.array( container.vel(2) )
                w = np.array( container.wgt() )

                mass = 1.0 if t in ["e-", "e+"] else 0.0

                ene = np.sqrt(mass*mass + ux**2 + uy**2 + uz**2) 
                enes[t] += np.sum( ene*w )*self.N_wgt/self.N_time/self.N_box # total energy carried by LP # TODO

                # optical depth
                if t in ['e-', 'e+']:
                    beta = sqrt(1.0 - 1.0/ene**2 + 1e-20)
                    taus[t] += np.sum(w*beta)*self.N_wgt #/self.N_box # N_tau # TODO

                # particle number
                nums[t] += np.sum(w)/conf.wsum0


        #--------------------------------------------------
        # MPI variables to rank0
        taus["e-"]    = MPI.COMM_WORLD.reduce(taus['e-'], op=MPI.SUM, root=0)
        taus["e+"]    = MPI.COMM_WORLD.reduce(taus['e+'], op=MPI.SUM, root=0)

        lp_nums["e-"] = MPI.COMM_WORLD.reduce(lp_nums['e-'], op=MPI.SUM, root=0)
        lp_nums["e+"] = MPI.COMM_WORLD.reduce(lp_nums['e+'], op=MPI.SUM, root=0)
        lp_nums["ph"] = MPI.COMM_WORLD.reduce(lp_nums['ph'], op=MPI.SUM, root=0)

        enes["e-"]    = MPI.COMM_WORLD.reduce(enes['e-'], op=MPI.SUM, root=0)
        enes["e+"]    = MPI.COMM_WORLD.reduce(enes['e+'], op=MPI.SUM, root=0)
        enes["ph"]    = MPI.COMM_WORLD.reduce(enes['ph'], op=MPI.SUM, root=0)

        nums["e-"]    = MPI.COMM_WORLD.reduce(nums['e-'], op=MPI.SUM, root=0)
        nums["e+"]    = MPI.COMM_WORLD.reduce(nums['e+'], op=MPI.SUM, root=0)
        nums["ph"]    = MPI.COMM_WORLD.reduce(nums['ph'], op=MPI.SUM, root=0)

        #--------------------------------------------------
        if self.root:
            self.storage.append('lap', lap)

            tau_tot_mean = (taus['e-'] + taus['e+'])/self.N_box
            self.storage.append('tau', tau_tot_mean)

            self.storage.append('lp_num_e-', lp_nums['e-'])
            self.storage.append('lp_num_e+', lp_nums['e+'])
            self.storage.append('lp_num_ph', lp_nums['ph'])

            self.storage.append('ene_e-', enes['e-'])
            self.storage.append('ene_e+', enes['e+'])
            self.storage.append('ene_ph', enes['ph'])

            self.storage.append('num_e-', nums['e-'])
            self.storage.append('num_e+', nums['e+'])
            self.storage.append('num_ph', nums['ph'])

        return


    def update_esc_stats(self, mc, lap, conf):

        # these values are accumulated for plot_interval period; calculate current val

        # energy carried away by LPs
        ene_esc = mc.leaked_ene*self.N_wgt/self.N_time/conf.plot_interval/self.N_box

        # TODO normalize these?
        ene_inj_ph = (self.N_wgt/self.N_time/self.N_box)*mc.inj_ene_ph/conf.plot_interval
        ene_inj_ep = (self.N_wgt/self.N_time/self.N_box)*mc.inj_ene_ep/conf.plot_interval

        lp_num_esc = mc.leaked_pnum/self.N_npc # LPs lost

        num_esc = mc.leaked_wsum/conf.wsum0

        #--------------------------------------------------
        # MPI
        lp_num_esc = MPI.COMM_WORLD.reduce(lp_num_esc, op=MPI.SUM, root=0)
        num_esc    = MPI.COMM_WORLD.reduce(num_esc,    op=MPI.SUM, root=0)
        ene_esc    = MPI.COMM_WORLD.reduce(ene_esc,    op=MPI.SUM, root=0)
        ene_inj_ep = MPI.COMM_WORLD.reduce(ene_inj_ep, op=MPI.SUM, root=0)
        ene_inj_ph = MPI.COMM_WORLD.reduce(ene_inj_ph, op=MPI.SUM, root=0)

        tau_tile_min = MPI.COMM_WORLD.reduce(self.tau_tile_min, op=MPI.MIN, root=0)
        tau_tile_max = MPI.COMM_WORLD.reduce(self.tau_tile_max, op=MPI.MAX, root=0)

        #--------------------------------------------------
        if self.root:

            # ver1: NOTE root has the correct (already MPI reduced) tau_global

            # ver2: NOTE value varies based on what tile ran comp_tau the last
            tau_meas = mc.tau_global # tau used in escape prob calc; 

            lum_in  = ene_inj_ph + ene_inj_ep #+ conf.lum_ant 
            lum_out = ene_esc
            lum_rat = lum_in/(lum_out + 1e-7)

            self.storage.append('lap_sparse', lap) # sparse sampling of laps
            self.storage.append('lp_num_esc', lp_num_esc)
            self.storage.append('num_esc', num_esc)
            self.storage.append('ene_esc', ene_esc)
            self.storage.append('ene_inj_ph', ene_inj_ph)
            self.storage.append('ene_inj_ep', ene_inj_ep)

            self.storage.append('lum_rat', lum_rat)

            self.storage.append('tau_meas', tau_meas)
            self.storage.append('tau_min', tau_tile_min)
            self.storage.append('tau_max', tau_tile_max)

        return



    def update_esc_hists(self, mc, lap, conf):


        self.h1_enes['esc'] = np.zeros(len(self.lnzs), dtype='d')

        hx_cnts_rank = np.array(mc.get_hist_cnts(), dtype='d')
        hx_cnts      = np.zeros(len(hx_cnts_rank), dtype='d') # global array

        ncnts = len(hx_cnts_rank)

        # reduce globally
        MPI.COMM_WORLD.Reduce(
            [hx_cnts_rank, ncnts, MPI.DOUBLE],
            [hx_cnts,      ncnts, MPI.DOUBLE],
            op = MPI.SUM,
            root = 0,
            )

        #--------------------------------------------------
        if self.root:
            self.h1_enes['esc'][:] = hx_cnts

            # units
            # NOTE: we normalize by the sampling interval = conf.plot_interval
            dx = self.dlnx*self.xs
            #self.h1_enes['esc'][:] *= self.xs*self.N_ene/conf.plot_interval
            self.h1_enes['esc'][:] *= (self.N_wgt/self.N_time)*self.xs**2/dx/conf.plot_interval/self.N_box

        return



    def update_hists(self, grid, lap, conf):

        #--------------------------------------------------
        # prepare arrays

        self.h1_enes = {}
        self.h1_nums = {}
        self.h1_ws   = {}
        for t in ['e-', 'e+']:
            self.h1_enes[t] = np.zeros(len(self.lnzs), dtype='d')
            self.h1_nums[t] = np.zeros(len(self.lnzs), dtype='d')
            self.h1_ws[  t] = np.zeros(len(self.lnzs), dtype='d')

        for t in ['ph']:
            self.h1_enes[t] = np.zeros(len(self.lnxs), dtype='d')
            self.h1_nums[t] = np.zeros(len(self.lnxs), dtype='d')
            self.h1_ws[  t] = np.zeros(len(self.lnxs), dtype='d')

        self.h2_nums = np.zeros((self.Nhist, self.Nhist), dtype='d')

        #--------------------------------------------------
        # height dependent momentum histogram

        self.h2_enes = {}
        for t in ['e-', 'e+']:
            self.h2_enes[t] = np.zeros( (len(self.hs), 2*len(self.lnzs)), dtype='d')

        for t in ['ph']:
            self.h2_enes[t] = np.zeros( (len(self.hs), 2*len(self.lnxs)), dtype='d')

        #--------------------------------------------------
        # local arrays for MPI
        # NOTE danger zone: dont mix self.arr and arr
        h1_enes = {}
        h1_nums = {}
        h1_ws   = {}
        for t in ['e-', 'e+', 'ph']:
            h1_enes[t] = np.zeros_like( self.h1_enes[t] ) #(len(self.lnzs), dtype='d')
            h1_nums[t] = np.zeros_like( self.h1_nums[t] ) #(len(self.lnzs), dtype='d')
            h1_ws[  t] = np.zeros_like( self.h1_ws[t]   ) #(len(self.lnzs), dtype='d')

        h2_nums = np.zeros((self.Nhist, self.Nhist), dtype='d')

        h2_enes = {}
        for t in ['e-', 'e+', 'ph']:
            h2_enes[t] = np.zeros_like(self.h2_enes[t])

        #--------------------------------------------------
        # analyze tiles

        for tile in pytools.tiles_local(grid):
            for ispcs in range(conf.Nspecies):
                container = tile.get_container(ispcs)
                t = container.type

                # choose bins depending on type
                if t in ["e-", "e+"]:
                    bins = self.zs_bin
                elif t in ["ph"]:
                    bins = self.xs_bin

                bins_w = self.ws_bin

                if conf.oneD:
                    h  = np.array( container.loc(0) ) # y
                    hv = np.array( container.vel(0) ) # vy
                elif conf.twoD:
                    h  = np.array( container.loc(1) ) # y
                    hv = np.array( container.vel(1) ) # vy
                elif conf.threeD:
                    h = np.array( container.loc(2) ) # z
                    hv= np.array( container.vel(2) ) # vz

                # compute energy
                ux = np.array( container.vel(0) )
                uy = np.array( container.vel(1) )
                uz = np.array( container.vel(2) )
                w  = np.array( container.wgt() )

                mass = 1.0 if t in ["e-", "e+"] else 0.0

                #ene = self.N_ene*np.sqrt(mass*mass + ux**2 + uy**2 + uz**2) 
                mom = np.sqrt(ux**2 + uy**2 + uz**2) # p = \gamma \beta OR x 

                # update histograms

                # physical energy spectra
                h1, edges_e = np.histogram(mom, bins, weights=w)
                h1_enes[t][:] += h1[:] 

                # LP weight energy spectra
                h1, edges_l = np.histogram(mom, bins)
                h1_ws[t][:] += h1[:]

                # LP weight spectra (for numerical analysis only)
                h1, edges_w = np.histogram(w, bins_w)
                h1_nums[t][:] += h1[:]

                if t == "ph":
                    # 2d histogram of weight and energy
                    h2, edges_w, edges_e = np.histogram2d(
                            w, mom, 
                            bins=(self.ws_bin, self.xs_bin),
                            weights=w,
                            )
                    h2_nums[:,:] += h2[:,:]

                #--------------------------------------------------
                # 2D histograms with spatial + energy info

                n_mid_h2_enes = len(bins) - 1

                # particles going up
                ind = np.where(hv > 0)
                h2, edges_h, edges_e = np.histogram2d(
                            h[ind], mom[ind], 
                            bins=(self.hs_bin, bins),
                            weights=w[ind],)
                h2_enes[t][:, n_mid_h2_enes:] += h2[:,:] # upper half

                # particles going down
                ind = np.where(hv < 0)
                h2, edges_h, edges_e = np.histogram2d(
                            h[ind], mom[ind], 
                            bins=(self.hs_bin, bins),
                            weights=w[ind],)
                h2_enes[t][:, :n_mid_h2_enes] += np.fliplr(h2[:,:]) # lower half



        #--------------------------------------------------
        # MPI
        # reduce globally

        for t in ['e-', 'e+', 'ph']:
            ncnts = len(h1_enes[t])
            MPI.COMM_WORLD.Reduce(
                [     h1_enes[t], ncnts, MPI.DOUBLE],
                [self.h1_enes[t], ncnts, MPI.DOUBLE],
                op = MPI.SUM,
                root = 0,)

            ncnts = len(h1_ws[t])
            MPI.COMM_WORLD.Reduce(
                [     h1_ws[t], ncnts, MPI.DOUBLE],
                [self.h1_ws[t], ncnts, MPI.DOUBLE],
                op = MPI.SUM,
                root = 0,)

            ncnts = len(h1_nums[t])
            MPI.COMM_WORLD.Reduce(
                [     h1_nums[t], ncnts, MPI.DOUBLE],
                [self.h1_nums[t], ncnts, MPI.DOUBLE],
                op = MPI.SUM,
                root = 0,)

            ncnts = 2*self.Nhist*self.Nhist
            MPI.COMM_WORLD.Reduce(
                [     h2_enes[t], ncnts, MPI.DOUBLE],
                [self.h2_enes[t], ncnts, MPI.DOUBLE],
                op = MPI.SUM,
                root = 0,)


        ncnts = self.Nhist*self.Nhist
        MPI.COMM_WORLD.Reduce(
                [     h2_nums[:], ncnts, MPI.DOUBLE],
                [self.h2_nums[:], ncnts, MPI.DOUBLE],
                op = MPI.SUM,
                root = 0,)


        #--------------------------------------------------
        if self.root:
            # make into units of d\tau/dp from original d\tau/dlnp = p d\tau/dp
            dz = self.dlnz*self.zs
            self.h1_enes['e-'][:] *= self.N_wgt/dz/self.N_box
            self.h1_enes['e+'][:] *= self.N_wgt/dz/self.N_box

            # make into units of p d\tau/\dp (note that this differs from h1_ene)
            # multiply both halfs of the spatial array with the same unit conversion factor
            dz2 = np.array( [np.flip(dz/self.zs), dz/self.zs ] ).flatten()
            self.h2_enes['e+'][:,:] *= self.N_wgt/dz2/self.N_box
            self.h2_enes['e-'][:,:] *= self.N_wgt/dz2/self.N_box

            # make into units of compactness; x d(x n_x) / dx from original d(n_x)/dlnx
            dx = self.dlnx*self.xs
            self.h1_enes['ph'][:] *= (self.N_wgt/self.N_time)*self.xs**2/dx/self.N_box

            # make into units of compactness as well; note that x^2 is multiplied in a different place
            # multiply both halfs of the spatial array with unit conversion factor
            # TODO make into same units as regular h1_enes array
            #      this needs a bit of thinking since array needs to be flipped in middle
            #dx2 = self.dlnx[1]
            dx2 = np.array( [np.flip(dx/self.xs**2), dx/self.xs**2 ] ).flatten()
            self.h2_enes['ph'][:,:] *= (self.N_wgt/self.N_time)/dx2/self.N_box

            # make into units of reference particle number
            ppc_tot = conf.npc_ref*conf.Nx*conf.Ny*conf.Nz
            dw = self.dlnw*self.ws
            self.h1_nums['e-'][:] *= 1.0/ppc_tot/dw
            self.h1_nums['e+'][:] *= 1.0/ppc_tot/dw 
            self.h1_nums['ph'][:] *= 1.0/ppc_tot/dw 

            # normalize 2d histogram of photons as well
            xx = (self.N_wgt/self.N_time)*self.xs**2/dx/self.N_box
            yy = 1.0/dw
            XX, YY = np.meshgrid(xx, yy, indexing='xy')
            self.h2_nums[:,:] *= XX*YY

            # same slope as regular energy specs; normalization to total prtcl num
            self.h1_ws['e-'][:] *= 1.0/ppc_tot/dz 
            self.h1_ws['e+'][:] *= 1.0/ppc_tot/dz 
            self.h1_ws['ph'][:] *= ppc_tot*self.xs**2/dx

        return


    def save(self, lap, conf):

        if self.root:
            fname = conf.outdir + '/qed_{}.h5'.format(str(lap)) #.rjust(4,'0'))
            #print(fname)
            f5 = h5.File(fname,'w')

            f5.create_dataset('h1_ene_e-', data=self.h1_enes['e-' ][:] )
            f5.create_dataset('h1_ene_e+', data=self.h1_enes['e+' ][:] )
            f5.create_dataset('h1_ene_ph', data=self.h1_enes['ph' ][:] )
            f5.create_dataset('h1_ene_esc',data=self.h1_enes['esc'][:] )

            f5.create_dataset('h1_nums_e-', data=self.h1_nums['e-' ][:] )
            f5.create_dataset('h1_nums_e+', data=self.h1_nums['e+' ][:] )
            f5.create_dataset('h1_nums_ph', data=self.h1_nums['ph' ][:] )

            f5.create_dataset('h1_ws_e-', data=self.h1_ws['e-' ][:] )
            f5.create_dataset('h1_ws_e+', data=self.h1_ws['e+' ][:] )
            f5.create_dataset('h1_ws_ph', data=self.h1_ws['ph' ][:] )

            f5.create_dataset('h2_ene_e-', data=self.h2_enes['e-' ][:,:] )
            f5.create_dataset('h2_ene_e+', data=self.h2_enes['e+' ][:,:] )
            f5.create_dataset('h2_ene_ph', data=self.h2_enes['ph' ][:,:] )

            #--------------------------------------------------
            for key in ['lap_sparse',
                        'lp_num_esc',
                        'num_esc',
                        'ene_esc',
                        'ene_inj_ph',
                        'ene_inj_ep',
                        'lum_rat',
                        'tau_meas',
                        'tau_min',
                        'tau_max',
                        'lap',
                        'tau',
                        'num_e-',
                        'num_e+',
                        'num_ph',
                        'lp_num_e-',
                        'lp_num_e+',
                        'lp_num_ph',
                        'ene_e-',
                        'ene_e+',
                        'ene_ph',
                        ]:
                f5.create_dataset(key, data=self.storage.data[key][:] )
            #--------------------------------------------------

            f5.close()

        return


