
#from configSetup import Configuration #old conf
from pytools import Configuration # new conf
import pytools  # runko python tools

from pytools import simplify_string, simplify_large_num

from numpy import sqrt, pi
import numpy as np



class Configuration_Pulsar(Configuration):

    def __init__(self, *file_names, do_print=False):
        Configuration.__init__(self, *file_names)

        #--------------------------------------------------
        # problem specific initializations
        if do_print:
            print("init:--------------------------------------------------")
            print("init: Initializing pulsar setup...")
            print("init: conf file:", file_names[0])

        #--------------------------------------------------
        # create directory naming scheme
        if self.outdir == "auto":

            # base name
            self.outdir = self.prefix if "prefix" in self.__dict__ else ""

            # grid size
            str_x = str(int(self.Nx))
            str_m = str(int(self.NxMesh))
            self.outdir += "x" + str_x + "m" + str_m + "_"

            #--------------------------------------------------
            # particle info
            self.outdir += "p" + simplify_string(self.ppc)
            self.outdir += "np" + simplify_string(self.npasses)
            self.outdir += "_"

            #--------------------------------------------------
            # pulsar
            self.outdir += "v" + simplify_string(self.vrot)
            self.outdir += "inje" + simplify_string(self.ninj_pairs)
            self.outdir += "injx" + simplify_string(self.ninj_phots)
            self.outdir += "_"
            self.outdir += "jm" + simplify_string(self.jm_scaling)
            self.outdir += "rc" + simplify_string(self.rad_curv)
            self.outdir += "_"

            #--------------------------------------------------
            # QED
            self.outdir += "b" + simplify_string(self.bratio)
            self.outdir += "N" + simplify_large_num(self.Nmp)
            self.outdir += "_"

            #--------------------------------------------------
            # postfix string
            if "postfix" in self.__dict__:
                self.outdir += self.postfix


        #--------------------------------------------------
        # end of directory naming

        if do_print:
            print("init: output dir", self.outdir)

        #--------------------------------------------------
        # turn off explicitly some advanced functionality
        self.use_injector = False
        self.c_corr = 1.0 # no speed of light correction
        self.use_maxwell_split = False

        # grid sizes
        self.Lx = self.Nx*self.NxMesh
        self.Ly = self.Ny*self.NyMesh
        self.Lz = self.Nz*self.NzMesh

        # re-name particle types
        self.prtcl_types = ['e-', 'e+', 'ph'] 

        #--------------------------------------------------
        #required cgs constants
        alphaf = 1/137

        # numerical Compton wavelength lambda_C
        self.lamC = 1.0/(4.0*pi*self.cfl**2*self.Nmp*alphaf)

        # 1-body QED process normalization constant = t_free/dt
        self.N_onebody = self.lamC
        self.N_onebody *= 1.0/alphaf/self.cfl #additionally, countering the normalized units in interact (where c=alpha=1)

        # 2-body QED process normalization constant = t_free/dt
        self.N_twobody = 6.0*pi*self.cfl**3*self.Nmp*self.NxMesh*self.NyMesh*self.NzMesh

        # super-stepping for the QED reactions 
        self.N_qdt = self.qed_step

        # gyroradius = lamC/b
        self.rg = self.lamC/self.bratio 

        # B normalization, r_g = lamC/b = mc^2/eB = c^2/B (in code units)
        self.bstar = self.cfl**2/self.rg
        self.binit = self.bstar/2.0 # B_{*,r} = radial magnetic field component

        #Schwinger field: b = B_*/B_Q
        self.B_QED = self.bstar/self.bratio

        # set polarcap and star size in cells
        if self.oneD:
            #self.rad_pcap = 0.6*self.Lx//4 #same size as would be for 3D sim (since there Nz = Nx/2)
            #self.rad_star = 10*self.rad_pcap

            # medium gap setup
            #self.rad_pcap = 0.8*self.Lx

            # large gap setup
            self.rad_pcap = 1.0*self.Lx
            self.rad_star = 2*self.rad_pcap #3*self.rad_pcap #10*self.rad_pcap
        else:
            sys.error() # not implemented


        # get electron charge from Goldreich-Julian number density
        # n_GJ = v B_*/e H = nppc 
        # with this definition, user-given ppc per species corresponds to n_GJ
        self.qe = -self.vrot*self.bstar/(self.ppc*self.rad_pcap)
        self.qi = -self.qe # positron charge

        # NOTE: local definitions of me and mi, the code initialiation uses a different definition
        me = np.abs(self.me)
        mi = np.abs(self.mi)
        me *= abs(self.qe) # update electron mass to right units
        mi *= abs(self.qe) # update electron mass to right units

        # photons
        self.qp = 0.0 # dummy charge; NOTE: needs to be 0 for photons to avoid charge deposit
        self.mp = 0.0 # dummy mass

        # temperatures
        self.delgam_e = self.delgam
        self.delgam_i = self.temp_ratio*self.delgam_e

        # magnetization corresponding to GJ number density
        self.sigma = self.bstar**2/(self.ppc*abs(me)*self.cfl**2)

        # pulsar parameters from user-given parameters
        self.period_star = 2*pi*self.rad_pcap/(self.vrot*self.cfl) # P_*
        self.Om_star = 2.0*pi/self.period_star # Omega_*
        self.rad_lcyl = self.cfl/self.Om_star # light cylinder distance in cells; not self consistent in PC setup


        # the height of the star's curving surface (in cells)
        #self.rad_curv_shift = self.rad_star - np.sqrt(self.rad_star**2 - self.rad_pcap**2) # height of the curved atmosphere
        #self.rad_curv_shift += 6 # pad with some cells to avoid boundary effects close to grid limit
        self.rad_curv_shift = 5 # NOTE using flat surface 

        self.height_atms = 1 #3 # add padding; this is the height of the atmosphere at r=Rpc in units of cells
        self.wph = self.wph


        # normalization coefficient for internal dipole coordinate system
        # defined so that we have B_* at the star's surface
        self.b_dipole_norm = self.binit*self.rad_star**3

        # rotation velocity (in units of c) of the polar cap disk; given by user
        #vrot = self.Om_star*self.rad_pcap/cfl

        # Goldreich-Julian density; given by user but recalculated
        self.nGJ = self.vrot*self.bstar/(abs(self.qe)*self.rad_pcap) # total ppc for both species

        # omega_p plasma frequency
        self.omp = sqrt(self.ppc*self.qe**2/abs(me))

        # skindepth in units of cell size
        self.c_omp = self.cfl/self.omp


        #--------------------------------------------------
        # gap parameters

        #polar cap height in code units, set to equal the polarcap radius
        self.h_pcap = self.rad_pcap 

        # maximum energy from the gap
        self.gam_gap = self.vrot*self.bstar*self.rad_pcap/(2*self.cfl**2) 

        # 1D curvature radius of a field line
        #self.rad_curv = self.rad_star**2/self.rad_pcap
        self.rad_curv = self.rad_curv*self.rad_pcap

        # gamma_rad, radiation reaction limit where gap gains equals radiation losses
        self.gam_rad_synch  = self.gam_gap**0.25
        self.gam_rad_synch *= (self.rg/self.rad_curv)**-0.5
        self.gam_rad_synch *= self.bratio**-0.5
        self.gam_rad_synch *= ( 1.5*self.lamC/self.h_pcap/alphaf)**0.25

        ninj_phots_per_cell = self.ninj_phots*self.xpc*self.wph
        self.gam_rad_comp = self.gam_gap**0.5
        self.gam_rad_comp *= self.delgam_x**-0.5 #m_e c^2 / kT = 1.0 / delgam_x
        self.gam_rad_comp *= (0.28*6.0*np.pi*self.cfl**5*self.Nmp/(max(ninj_phots_per_cell,1)*self.h_pcap))**0.5

        # radiation length (distance that the particle travels before reaching the radiation limit)
        self.len_rad = (self.gam_rad_synch/self.gam_gap)*self.h_pcap # simpler v2

        # characteristic synchrotron photon energy
        self.xsyn = 1.5*self.bratio*(self.rg/self.rad_curv)*self.gam_rad_synch**3

        # 1-photon pair creation mean free path

        lmfp1_per_h  = self.rad_curv/(self.xsyn*self.bratio) # prefactor
        lmfp1_per_h *= 1.0/( np.log(0.333*(alphaf*self.rad_curv/( self.lamC*self.bratio*self.xsyn**2 ) )**(3/8) )
                            - np.log( np.log(0.333*(alphaf*self.rad_curv/( self.lamC*self.bratio*self.xsyn**2 ) )**(3/8) )) )
        lmfp1_per_h *= 1.0/self.rad_pcap # into units of gap height

        #--------------------------------------------------
        # Calculation of 2-photon pair creation mean free path

        hplanck = 6.6261e-27 #erg s
        c = 2.9979e10 #cm/s
        kev2erg = 1.602176633e-9
        mec2_inkev = 511.0
        Tkev = 0.5
        kTperhc = kev2erg*Tkev/(hplanck*c)
        nBB = 16.0*np.pi*1.202*kTperhc**3
        x2_av_ene = 2.7*Tkev/mec2_inkev

        comp_scale = 6.0*np.pi*self.cfl**5*self.Nmp/(max(ninj_phots_per_cell,1)*self.h_pcap)
        x1 = 1e5
        x2 = x2_av_ene
        xpr = x1*x2
        efac = 0.652*(xpr**2-1.0)*np.log(xpr)*np.heaviside(xpr-1.0,1.0)/xpr**3
        efac = 1.0/efac
        lmfp2_per_h = comp_scale*efac
        #print("comp_scale: ", comp_scale)
        #print("lmfp2_per_h:", lmfp2_per_h)

        #--------------------------------------------------
        # extra undefined parameters
        # TO BE REMOVED
        H = 1.0e5
        c = 3e10
        self.t_c = H/c # light crossing time across the system
        self.xi = 1.0 #self.t0 # TODO setting this automatically to match PIC 
        self.dt = self.xi*self.t_c # time step; we resolve mean interaction time by xi
        ppc0 = self.ppc # FIXME free to choose this
        self.npc_ref = self.NxMesh*self.NyMesh*self.NzMesh*ppc0
        self.N_box = self.Nx*self.Ny*self.Nz
        self.N_time = self.dt/self.t_c # time step in units of light crossing time
        self.N_wgt = 1.0

        self.Nref = {}
        self.Nref["ph"] = self.npc_ref # reference photon num dens
        self.Nref["e-"] = self.npc_ref # reference electron num dens 
        self.wsum0 = self.Nref["e-"] # reference weight


        #--------------------------------------------------
        if do_print:
            print("conf:              b:", self.bratio)
            print("conf:            Nmp:", self.Nmp)
            print("conf:            cfl:", self.cfl)
            print("conf:            ppc:", self.ppc)

            print("star:")
            print("star:          v_rot:", self.vrot)
            print("star:            P_*:", self.period_star)
            print("star:           Om_*:", self.Om_star)
            print("star:           R_LC:", self.rad_lcyl)
            print("star:           R_pc:", self.rad_pcap)
            print("star:            R_*:", self.rad_star)
            print("star:         B_norm:", self.b_dipole_norm)
            print("star: rad_curv_shift:", self.rad_curv_shift)
            print("star:    height_atms:", self.height_atms)
            print("star:     r_g/R_curv:", self.rg/self.rad_curv)
            print("star:         R_curv:", self.rad_curv)

            print("init:")
            print("init:           lamC:", self.lamC)
            print("init:             rg:", self.rg)
            print("init:          sigma:", self.sigma)
            print("init:          binit:", self.binit)
            print("init:            B_*:", self.bstar)
            print("init:            nGJ:", self.nGJ)
            print("init:             qe:", self.qe, me)
            print("init:             qi:", self.qi, mi)
            print("init:          m-/me:", self.me)
            print("init:          m+-me:", self.mi)
            print("init:        omega_p:", self.omp)
            print("init:    R = c/omega:", self.c_omp, " dx")
            print("init: R_pc / c/omega:", self.rad_pcap/self.c_omp)

            print("time")
            print("time:         t_esc:", self.rad_pcap/self.cfl, " laps")

            print("phys:")
            print("phys:        gam_gap:", self.gam_gap)
            print("phys:    gam_rad_syn:", self.gam_rad_synch)
            print("phys:   gam_rad_comp:", self.gam_rad_comp)
            print("phys:    g_syn/g_gap:", self.gam_rad_synch/self.gam_gap)
            print("phys:    g_com/g_gap:", self.gam_rad_comp/self.gam_gap)

            print("phys:  1-photon mfp/H:", lmfp1_per_h)
            print("phys:  2-photon mfp/H:", lmfp2_per_h)

            print("phys:         len_rad:", self.len_rad)
            print("phys:           xcurv:", self.xsyn)


        #if(self.gam_rad > self.gam_rad_comp):
        #    self.gam_rad = self.gam_rad_comp
        self.gam_rad = self.gam_rad_comp
        #--------------------------------------------------
        # default normalization

        self.b_norm = self.bstar 
        self.e_norm = self.bstar*self.vrot
        self.j_norm = abs(self.qe)*self.nGJ*self.cfl**2 # j_m \Delta t # abs(self.qe)*self.ppc*2*self.cfl**2
        self.p_norm = self.ppc 
        self.x_norm = max(self.xpc,1)
        self.t_norm = self.rad_pcap/self.cfl  # t_pc = R_pc/c; lightcrossing time across polar cap


        # DONE
        if do_print:
            print("init:--------------------------------------------------")

# end of class


# Prtcl velocity (and location modulation inside cell)
#
# NOTE: Cell extents are from xloc to xloc + 1, i.e., dx = 1 in all dims.
#       Therefore, typically we use position as x0 + RUnif[0,1).
#
# NOTE: Default injector changes odd ispcs's loc to (ispcs-1)'s prtcl loc.
#       This means that positrons/ions are on top of electrons to guarantee
#       charge conservation (because zero charge density initially).
#
def velocity_profile(xloc, ispcs, conf):

    # electrons
    if ispcs == 0:
        delgam = conf.delgam_e  # * np.abs(conf.mi / conf.me) * conf.temp_ratio

    # positrons/ions/second species
    elif ispcs == 1:
        delgam = conf.delgam_i

    # photons
    elif ispcs == 2:
        delgam = conf.delgam_x

    # perturb position between x0 + RUnif[0,1)
    xx = xloc[0] + np.random.rand()
    yy = xloc[1] + np.random.rand()
    zz = xloc[2] + np.random.rand()

    # velocity sampling from Maxwell-Juttner
    if ispcs in [0,1]:
        gamma = 0 # no bulk motion
        direction = +1
        ux, uy, uz, uu = pytools.sample_boosted_maxwellian(
            delgam, gamma, direction=direction, dims=3
        )
    elif ispcs == 2: # photons
        ux, uy, uz, uu = pytools.sample_blackbody(delgam)

        #if np.isnan(ux) or np.isnan(uy) or np.isnan(uz) or np.isnan(uu):
        #    sys.exit()

    # TODO debug
    #ux = 1e3

    x0 = [xx, yy, zz]
    u0 = [ux, uy, uz]
    return x0, u0



# number of prtcls of species 'ispcs' to be added to cell at location 'xloc'
#
# NOTE: Plasma frequency is adjusted based on conf.ppc (and prtcl charge conf.qe/qi
#       and mass conf.me/mi are computed accordingly) so modulate density in respect
#       to conf.ppc only.
#
def density_profile(xloc, ispcs, conf):

    return 0 # NOTE no injection in the beginning of the simulation
    
    # TODO debug
    #if xloc[0] > conf.height_atms + 11:
    #    return 0
    #if xloc[0] < conf.height_atms + 10:
    #    return 0

    if ispcs in [0,1]:
        return conf.ppc
    if ispcs == 2:
        return conf.xpc

# Particle weight that is added to cell at location 'xloc'
def weigth_profile(xloc, ispcs, conf):
    return 1.0


