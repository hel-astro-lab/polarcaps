
import numpy as np
import pytools  # runko python tools
import sys


from numpy import pi, sqrt, exp
from mpi4py import MPI
from numpy import newaxis as nax


# 3D Langevin antenna that drives a current
class Antenna:
    def __init__(self, min_mode, max_mode, conf):
        # print("initializing...")

        self.cfl = conf.cfl

        self.Nx = conf.Nx
        self.Ny = conf.Ny
        self.Nz = conf.Nz

        self.NxMesh = conf.NxMesh
        self.NyMesh = conf.NyMesh
        self.NzMesh = conf.NzMesh

        # local current storage
        self.jx = np.zeros((self.NxMesh, self.NyMesh, self.NzMesh))
        self.jy = np.zeros((self.NxMesh, self.NyMesh, self.NzMesh))
        self.jz = np.zeros((self.NxMesh, self.NyMesh, self.NzMesh))

        lx = conf.Nx*self.NxMesh #/ conf.c_omp
        ly = conf.Ny*self.NyMesh #/ conf.c_omp
        lz = conf.Nz*self.NzMesh #/ conf.c_omp

        self.jnorm = conf.jm_scaling*conf.qe*conf.ppc*conf.cfl

        # physical timer inside the antenna to track measured time
        self.tcur = 0.0

        return


    def add_ext_cur(self, tile):
        g = tile.get_grids(0)
        
        for r in range(self.NzMesh):
            for s in range(self.NyMesh):
                for q in range(self.NxMesh):
                    g.jx[q,s,r] += self.jnorm                
                    #g.jx[q,s,r] += -self.jnorm
                    #g.jy[q,s,r] += self.jy[q,s,r]
                    #g.jz[q,s,r] += self.jz[q,s,r]

        return


