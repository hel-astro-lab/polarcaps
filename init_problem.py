
#from configSetup import Configuration #old conf
from pytools import Configuration # new conf
import pytools  # runko python tools

from pytools import simplify_string, simplify_large_num

from numpy import sqrt, pi
import numpy as np



class Configuration_Turbulence(Configuration):

    def __init__(self, *file_names, do_print=False):
        Configuration.__init__(self, *file_names)

        #--------------------------------------------------
        # problem specific initializations
        if do_print:
            print("init:--------------------------------------------------")
            print("init: Initializing turbulence setup...")
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
            # ppc
            str_ppc  = simplify_string(self.ppc)
            str_np   = simplify_string(self.npasses)
            str_comp = simplify_string(self.c_omp)
            self.outdir += "p" + str_ppc + "np" + str_np + "c" + str_comp + "_"

            #--------------------------------------------------
            # sigma
            str_sig = simplify_large_num(self.sigma)
            str_delg = simplify_string(self.delgam)
            self.outdir += "s" + str_sig + "d" + str_delg + "_"


            #--------------------------------------------------
            # radiation 
            if "gammarad" in self.__dict__:
                if self.gammarad > 0.0:
                    str_r = str(int(self.gammarad))
                    self.outdir += "r" + str_r + "_"

            #--------------------------------------------------
            # driving  amplitude
            if "drive_ampl" in self.__dict__:
                str_a = simplify_string(self.drive_ampl)
                self.outdir += "a" + str_a + "_"

            if "drive_freq" in self.__dict__:
                str_f = simplify_string(self.drive_freq)
                self.outdir += "w" + str_f + "_"

            if "decorr_time" in self.__dict__:
                str_d = simplify_string(self.decorr_time)
                self.outdir += "g" + str_d + "_"

            #--------------------------------------------------
            if self.qed_mode:
                if "lum_ph1" in self.__dict__ and "lum_ph2" in self.__dict__:
                    str_lumx1 = simplify_string(self.lum_ph1)
                    str_lumx2 = simplify_string(self.lum_ph2)

                    self.outdir += "lx" + str_lumx1 + "t" + simplify_large_num(self.kTbb1) 
                    self.outdir += "_" 

                    self.outdir += "ld" + str_lumx2 + "t" + simplify_large_num(self.kTbb2) 
                    self.outdir += "_"

                if "tau_ext" in self.__dict__ and "tau_ini" in self.__dict__:
                    str_t1 = simplify_string(self.tau_ini)
                    str_t2 = simplify_string(self.tau_ext)
                    self.outdir += "taui" + str_t1 + "_taue" + str_t2 + "_"

                if "H" in self.__dict__:
                    self.outdir += "H" + simplify_large_num(self.H) 


            #--------------------------------------------------
            # postfix string
            if "postfix" in self.__dict__:
                self.outdir += self.postfix

        if do_print:
            print("init: output dir", self.outdir)


        #--------------------------------------------------
        # turn off some advanced functionality
        self.use_injector = False
        self.c_corr = 1.0 # no speed of light correction
        self.use_maxwell_split = False

        # local variables just for easier/cleaner syntax
        me = np.abs(self.me)
        mi = np.abs(self.mi)
        cfl = self.cfl
        ppc = self.ppc * 2.0  # multiply x2 to account for 2 species/pair plasma

        # plasma reaction & subsequent normalization
        self.gamma = 1.0
        self.omp = cfl / self.c_omp
        self.qe = -(self.omp ** 2.0 * self.gamma) / ((ppc * 0.5) * (1.0 + me / mi))
        self.qi = -self.qe

        me *= abs(self.qi)
        mi *= abs(self.qi)

        # temperatures
        self.delgam_e = self.delgam
        self.delgam_i = self.temp_ratio*self.delgam_e

        if do_print:
            print("init: Positron thermal spread: ", self.delgam_i)
            print("init: Electron thermal spread: ", self.delgam_e)

        sigmaeff = self.sigma #* temperature corrections

        if do_print:
            print("init: Alfven vel:    ",sqrt(sigmaeff/(1.+sigmaeff)))
            print("init: Ion beta:      ",     2.*self.delgam_i/(self.sigma*(mi/me+1.)/(mi/me)) )
            print("init: Electron beta: ",2.*self.delgam_e/(self.sigma*(mi/me+1.)) )

        #--------------------------------------------------
        # photons
        self.qp = 0.0 # dummy charge; NOTE: needs to be 0 for photons to avoid charge deposit
        self.mp = 0.0 # dummy mass
        #self.delgam_x = 0.1 #  black body temperature

        # re-name particle types
        self.prtcl_types = ['e-', 'e+', 'ph'] 


        #--------------------------------------------------
        # get magnetic field from sigma: various definitions exist

        # no corrections; cold sigma
        self.binit_nc = sqrt(ppc*cfl**2.*self.sigma*me)

        #self.gammath = 1.0 + 3.0*delgam_e # approximative \gamma_th = 1 + 3\theta
        self.gammath = 1.0 + (3.0/2.0)*self.delgam_e # another approximation which is more accurate at \delta ~ 1

        #self.gammath = 1.55 #manual value for theta=0.3
        self.binit_approx = sqrt(self.gammath*ppc*me*cfl**2.*self.sigma)

        #self.binit = sqrt((self.gamma)*ppc*.5*c**2.*(me*(1.+me/mi))*self.sigma)

        self.binit = self.binit_approx # NOTE: selecting this as our sigma definitions

        self.omB = self.omp*sqrt(self.sigma) # nominal gyrofrequency

        if do_print:
            print("init: sigma:            ", self.sigma)
            print("init: mass term:        ", sqrt(mi+me))
            print("init: q_e:              ", self.qe)
            #print("init: warm corr:: ", self.warm_corr)
            print("init: gamma_th:         ", self.gammath)
            print("init: B_guide (no corr):", self.binit_nc)
            print("init: B_guide (approx): ", self.binit_approx)
            print("init: B_guide (used):   ", self.binit)


        #forcing scale in units of skin depth
        #self.l0 = self.Nx*self.NxMesh/self.max_mode/self.c_omp  
        
        #thermal larmor radius
        lth = self.gammath / np.sqrt(self.sigma)*np.sqrt(self.gammath) #gamma_0

        #reconnection radius; gam ~ sigma
        lsig = self.sigma / np.sqrt(self.sigma)*np.sqrt(self.gammath) 

        # maximum attainable particle energy; from gyro-scale resonance with largest scale
        #self.g0 = self.l0*np.sqrt(self.sigma)*np.sqrt(self.gammath) 


        if do_print:
            #print("init: driving scale l_0:  ", self.l0)
            print("init: larmor radi r_g  :  ", lth)
            print("init: larmor rad post rec:", lsig)
            #print("init: max gamma:          ", self.g0)


        #--------------------------------------------------
        # adding of static background component to particle pusher
        self.bpar = 1
        self.bperp = 1
        self.bplan = 1

        self.Lx = self.Nx*self.NxMesh
        self.Ly = self.Ny*self.NyMesh
        self.Lz = self.Nz*self.NzMesh

        # external fields
        if self.use_maxwell_split:
            self.bx_ext = self.binit*self.bpar  #np.cos(bphi)
            self.by_ext = self.binit*self.bplan #self.binit * bperp #np.sin(bphi) * np.sin(btheta)
            self.bz_ext = self.binit*self.bperp #np.sin(bphi) * np.cos(btheta)

            self.ex_ext = 0.0
            self.ey_ext = -self.beta * self.binit*self.bperp
            self.ez_ext = +self.beta * self.binit*self.bplan
        else:
            self.bx_ext = 0.0
            self.by_ext = 0.0
            self.bz_ext = 0.0

            self.ex_ext = 0.0
            self.ey_ext = 0.0
            self.ez_ext = 0.0

        #-------------------------------------------------- 
        # radiative cooling
        if "gammarad" in self.__dict__:
            if not(self.gammarad == 0.0):
                self.drag_amplitude = 0.1*self.binit/(self.gammarad**2.0)

                sigma_perp = self.sigma*self.drive_ampl**2
                self.g0 = self.l0*np.sqrt(sigma_perp)*np.sqrt(self.gammath) #gamma_0
                A = 0.1*self.g0/self.gammarad**2 #definition we use in Runko

                if do_print:
                    print("init: using radiation drag...")
                    print("init: drag amplitude: {} with gamma_rad: {}".format(self.drag_amplitude, self.gammarad))
                    print("init: rad A:        ", A)


        #-------------------------------------------------- 
        # driving
        lx = self.Nx*self.NxMesh 
        #self.t0 = 1.0/(lx/self.cfl/self.max_mode) # time steps in units of light-crossing 

        #if do_print:
        #    print("init:--------------------------------------------------")
        #    print("init: using Langevin antenna...")
        #    print("init: ampl: ", self.drive_ampl)
        #    print("init: w0:   ", self.drive_freq)
        #    print("init: g0:   ", self.decorr_time)
        #    print("init: t0:   ", self.t0, " laps")


        #-------------------------------------------------- 
        # polar cap setup

        if self.oneD:
            self.rad_pcap = 0.6*self.Lx//4 #same size as woudld be for 3D sim (since there Nz = Nx/2)
            self.rad_star = 10*self.rad_pcap

        else:
            if False: # radio pulsar setup
                self.rad_pcap = 0.6*self.Lx//2 #0.65 is exactly at the corner for 2d domain when R_* = 10 R_pc
                self.rad_star = 10*self.rad_pcap

            else: # millisecond pulsar setup
                self.rad_pcap = 0.4*self.Lx//2 
                self.rad_star = 4*self.rad_pcap


        # v1 based on polar cap size we get the star's spin 
        #self.rad_lcyl = self.rad_star*(self.rad_star/self.rad_pcap)**2
        #self.Om_star = self.cfl/self.rad_lcyl # \Omega_*
        #self.period_star = 2*np.pi/self.Om_star # P_*
        #self.vrot  = self.Om_star*self.rad_pcap/self.cfl # v_\phi of the polar cap edge

        # v2 based on given vrot infer other parameters
        self.period_star = 2*pi*self.rad_pcap/(self.vrot*self.cfl)
        self.Om_star = 2.0*np.pi/self.period_star
        self.rad_lcyl = self.cfl/self.Om_star # light cylinder distance in cells; not self consistent in PC setup

        self.t0 = self.cfl/self.rad_pcap # time steps in units of light-crossing across the polar cap

        self.rad_curv_shift = self.rad_star - np.sqrt(self.rad_star**2 - self.rad_pcap**2) # height of the curved atmosphere
        self.rad_curv_shift += 6 # pad with some cells to avoid boundary effects close to grid limit

        if self.oneD:
            self.rad_curv_shift = 5 # NOTE using flat surface in 3D
        if self.threeD:
            self.rad_curv_shift = 5 # NOTE using flat surface in 3D

        self.height_atms = 3 # add padding; this is the height of the atmosphere at r=Rpc in units of cells

        #self.chi = 10 # magnetic obliquity angle
        phase = 0.0 #global rotator phase

        self.b_dipole_norm = self.binit*self.rad_star**3
        self.bstar = 2*self.b_dipole_norm*self.rad_star**-3 # B_{*,r} = radial magnetic field component 
                                                       # at the star's surface; 
                                                       # factor of 2 comes from dipole coordinate system
        # NOTE: using constant B in 1D
        #if self.oneD: self.b_dipole_norm = self.binit

        vrot = self.Om_star*self.rad_pcap/cfl

        #omB = sqrt(self.sigma)*self.omp
        #dV_gap = omB*(self.rad_star/self.cfl)*(self.rad_star/self.rad_lcyl)**2 # DONE
        #print('init: dV_gap:', dV_gap)

        #self.nGJ = self.Om_star*bstar/(cfl*abs(self.qe))
        self.nGJ = vrot*self.bstar/(abs(self.qe)*self.rad_pcap) # total ppc for both species

        #nGJ = self.Om_star*bstar/(2*pi*self.cfl*abs(self.qe))

        deGJ = self.c_omp*(self.nGJ/self.ppc)**(-0.5)

        gam_gap  = vrot*self.bstar*self.rad_pcap/(2*cfl**2) # OK
        #gam_gap  = vrot**2*(bstar*self.rad_pcap/(cfl**2))
        #gam_gap2 = vrot**2*(bstar*self.rad_star/(cfl**2))
        #gam_gap3 = (self.omB/self.Om_star)*(self.rad_star/self.rad_lcyl)**3

        if do_print:
            print(' pulsar initialization...')
            print('init: P_*:   ', self.period_star)
            print('init: Om_*:  ', self.Om_star, " f_*:",self.Om_star/(2*pi) )
            print('init: R_*:   ', self.rad_star)
            print('init: R_pc:  ', self.rad_pcap)
            print('init: R_LC:  ', self.rad_lcyl, ' not consistent for PC setup')
            print('init: chi:   ', np.deg2rad(self.chi))
            print('init: B_*    ', self.bstar)
            #print('init: B_LC   ', bstar*(self.rad_lcyl/self.rad_star)**(-3.0))
            print('init: v_rot  ', vrot)
            print('init: n_GJ   ', self.nGJ)
            print('init: d_e GJ ', deGJ, 'dx')
            print('init: gam_gap', gam_gap)

        # add variables to the class
        self.gam_gap = gam_gap

            #sys.exit()


        #-------------------------------------------------- 
        # default normalization

        self.b_norm = 2.0*self.binit # B_0 = 2*binit; hence 2x
        self.e_norm = 2.0*self.binit*self.vrot
        self.j_norm = abs(self.qe)*self.nGJ*self.cfl**2 # j_m \Delta t # abs(self.qe)*self.ppc*2*self.cfl**2
        self.p_norm = self.nGJ # max(self.ppc*2,1)
        self.x_norm = max(self.xpc,1)
        self.t_norm = self.rad_pcap/self.cfl  # t_pc = R_pc/c; lightcrossing itme across polar cap


        #-------------------------------------------------- 
        # QED definitions
        if self.qed_mode:

            if do_print:
                print("init:--------------------------------------------------")

            me   = 9.1093826e-28   # electron mass g
            sigT = 6.65245873e-25  # Thomson optical depth; cm**2
            c    = 2.99792458e10   # speed of light cm/s
            kB   = 1.3806505e-16   # erg/K
            re   = 2.82e-13        # classical electron radius

            # QED reaction rate normalization
            #ppc0 = max(2*self.ppc, self.xpc)
            #ppc0 = max(8, ppc0) # min limiter
            #ppc0 = 1.0 # FIXME free to choose this
            ppc0 = self.ppc # FIXME free to choose this
            self.npc_ref = self.NxMesh*self.NyMesh*self.NzMesh*ppc0


            #--------------------------------------------------
            # inject

            # synchrotorn seed photon
            self.enebb1 = 2.7*self.kTbb1

            # disk photons
            self.enebb2 = 2.7*self.kTbb2

            self.Nref = {}
            self.Nref["ph"] = self.npc_ref # reference photon num dens
            self.Nref["e-"] = self.npc_ref # reference electron num dens 
            self.wsum0 = self.Nref["e-"] # reference weight

            #--------------------------------------------------
            H = self.H          # characteristic physical system size
            R = self.H          # characteristic physical system radius (assuming cube so not needed) 
            vol = R*R*H         # physical system volume
            dx_phys = H/self.Lx # physical size of \Delta x

            #self.nep0  = 1e17 #
            self.nz0  = self.tau_ini/(sigT*H) # number density required for optical depth unity
            #tau0       = sigT*self.nep0*self.H # reference optical depth
            self.t_c   = H/c # light crossing time across the system
            self.t_tau = self.t_c/self.tau_ext if self.tau_ext > 0 else self.t_c # minimum time between interactions 

            self.Nmp   = dx_phys/(cfl**2*ppc0)/(4.0*pi*re) # number of real electrons in a computational macro particle
            self.Nmp2  = self.nz0*dx_phys**3/ppc0 # required macro-particle number to model the given optical depth

            #self.xi = 1e-3 #0.01 # numerical safety factor for resolving QED interactions
            self.xi = 1 #self.t0 # TODO setting this automatically to match PIC 
            self.dt = self.xi*self.t_c # time step; we resolve mean interaction time by xi

            #self.dt *= self.qed_step # increase step size based on interval # TODO

            if do_print:
                print("init: t_0/t_c   : ", self.t0/self.t_c)
                print("init: xi        : ", self.xi)
                print("init: dt_QED    : ", self.dt)
                print("init: dt_PIC    : ", self.t0)
                print("init: dx_phys   : ", dx_phys)
                print("init: r_e (code): ", -self.qe/self.cfl**2)
                print("init: r_e (phys): ", re)
                print("init: r_eC/r_eP : ", -self.qe/self.cfl**2/re)


            # normalizations 
            self.N_wgt = self.nz0*sigT*H/self.wsum0 # unit conversion factor
            self.N_wgt *= 0.5 # normalize with two species pair plasma
            #self.N_w *= 1.0/(self.Nx*self.Ny*self.Nz) # FIXME normalize with total box size
            # NOTE N_wgt uses wsum0 tile ppc as a reference; hence it already has 1/N_box

            # tau is sum over tiles; we normalize to get mean tau per tile
            #self.N_tau = self.N_w
            #self.N_tau *= 1.0/(self.Nx*self.Ny*self.Nz) # FIXME

            self.N_time = self.dt/self.t_c # time step in units of light crossing time
            self.N_box = self.Nx*self.Ny*self.Nz
            self.N_qdt = self.qed_step     # QED reaction time steps to plasma time steps  

            # NOTE: it then follows that unit of luminosity is N_wgt / N_time

            # normaliation of onebody interaction; reduced Compton wavelength

            # this is inverse of tau_C = lamC/c = alpha_f r_e /c = e^2/mc^2 = e/c^2
            #self.N_lamC = 1.0e3 #H/self.Rpc # TODO
            #self.N_lamC = 137*self.cfl**2/abs(self.qe)
            #self.N_lamC2= 137*self.c_omp**2*self.ppc/self.cfl

            self.qed_vir_ampl = 1e12 # this becomes inefficient at >1e14; too rare for anything to happen
            self.N_onebody = re/(self.cfl*dx_phys) # reference normalization that would be self-consistent with the grid size
            self.N_onebody *= self.qed_vir_ampl # 1e4 # artificial amplification factor; = size of r_e in \Delta x

            self.lamC = self.N_onebody*self.cfl # \lam_C in units of \Delta x 
            #self.lamC *= 1/137 # (\alpha_f) NOTE: not needed (for some reason)

            #self.lamC = self.N_onebody*self.cfl*dx_phys*137
            #print('lamC', self.lamC)
            #print('lamC2', self.qe/(self.cfl**2)/137)

            #self.lamC2 = 2*pi*self.B_QED*self.c_omp/sqrt(self.sigma)
            #print('lamC1', self.lamC)
            #print('lamC2', self.lamC2)

            # radiation reaction limit \gamma_rad
            # v1: gives weird results
            #self.gam_rad  = (137*3/2)**0.25
            #self.gam_rad *= (self.vrot*self.B_QED)**(0.25)
            #self.gam_rad *= (self.rad_star**2/self.rad_pcap/self.lamC)**0.5

            # v2 wrong
            #self.gam_rad  = self.gam_gap**(0.25) 
            #self.gam_rad *= self.rad_star/self.rad_pcap
            #self.gam_rad *= ( 3*self.rad_pcap/(137*self.lamC) )**0.25

            #v3
            self.h_pcap = self.rad_pcap # polar cap height

            self.rad_curv = self.rad_star**2/self.rad_pcap
            rg = self.c_omp/sqrt(self.sigma)

            self.gam_rad  = self.gam_gap**0.25
            self.gam_rad *= (rg/self.rad_curv)**-0.5
            self.gam_rad *= self.B_QED**-0.5
            self.gam_rad *= ( 1.5*self.lamC/(137*self.h_pcap) )**0.25

            #self.gam_rad *= ( 1.5*self.lamC/(137*self.H) )**0.25
            #self.gam_rad *= ( 1.5*137*re/self.H )**0.25

            #self.gam_rad *= ( dx_phys*137 )**0.25 # arbitrary correction factor

            # synchrotorn radiation cooling happens over a distance (in units of polar cap size)
            self.len_rad = self.gam_rad*(1/self.B_QED/vrot)*self.lamC/self.rad_pcap

            #self.xsyn = 1.5*(self.lamC/self.rad_star)*(self.rad_pcap/self.rad_star)*self.gam_rad**3
            self.xsyn = 1.5*self.B_QED*(rg/self.rad_curv)*self.gam_rad**3

            #self.N_inj = self.N_w*self.t_c/self.dt # unit of injection luminosity (per dt) DONE
            #self.N_inj *= self.Nx*self.Ny*self.Nz # FIXME

            # how much to inject as an initial solution
            #self.ppc_ini = int( 0.25* ppc0 )
            #if do_print: print('ppc_ini', self.ppc_ini)

            #--------------------------------------------------
            #self.N_Q = (self.t_c/self.dt)/self.N_w   # TODO
            #self.N_Q *= 1.0/(self.Nx*self.Ny*self.Nz) # scale up
            #self.N_Q = self.N_wgt*self.N_time


            if do_print:
                print('init: n_pm   :', np.log10(self.nz0))
                print('init: H      :', np.log10(self.H))
                print('init: tau    :', self.tau_ext)
                print('init: t_c    :', self.t_c)
                print('init: t_tau  :', self.t_tau)
                print('init: N_mp   : {:.2e}'.format(self.Nmp))
                print('init: N_mp,t : {:.2e}'.format(self.Nmp2))
                print('init: dt     :', self.dt, 'dt/t_tau:', self.dt/self.t_tau, 'dt/tc:', self.dt/self.t_c)
                print('init: N_ph   :', self.Nref['ph'])
                print('init: N_ep   :', self.Nref['e-'])
                print('init: N_wgt  :', self.N_wgt)
                print('init: N_time :', self.N_time)
                print('init: N_1bQ  :', self.N_onebody)
                print('init: lamC/Hp: {:.2e}'.format(self.lamC/self.h_pcap))
                print('init: x_syn  : {:.2e}'.format(self.xsyn))
                print('init: rg/R_c : {:.2e}'.format(rg/self.rad_curv))
                print('init: gam_rad: {:.2e}'.format(self.gam_rad))
                print('init: len_rad: {:.2e}'.format(self.len_rad))
                #print('init: N_lamC :', self.N_lamC)
                #print('init: N_lamC2:', self.N_lamC2)
                #print('init: N_Q    :', self.N_Q)
                #print('init: N_inj  :', self.N_inj)
                print('init: wsum0  :', self.wsum0)
                #print('L:     :', lum_units)
                print('init: esc con1:', 'dt/tc < 1:', self.dt/self.t_c)
                print('init: esc con1:', 'dt/tt < 1:', self.dt/self.t_tau)


            #inject; total number of injected photons

            # required weight to get lum_ph
            if True:
                if self.zeta_xinj1 > 0 and self.lum_ph1 > 0:
                    self.Nph_inj1 = self.zeta_xinj1*self.wsum0 #self.Nref["ph"]
                    #self.wph_inj1 = self.lum_ph1/(self.Nph_inj1*self.enebb1*self.N_wgt/self.N_time*self.N_box)
                    self.wph_inj1 = self.lum_ph1/(self.Nph_inj1*self.enebb1*self.N_wgt/self.N_time)
                else:
                    self.Nph_inj1 = 0
                    self.wph_inj1 = 1

                if self.zeta_xinj2 > 0 and self.lum_ph2 > 0:
                    self.Nph_inj2 = self.zeta_xinj2*self.wsum0 #self.Nref["ph"]
                    #self.wph_inj2 = self.lum_ph2/(self.Nph_inj2*self.enebb2*self.N_wgt/self.N_time*self.N_box)
                    self.wph_inj2 = self.lum_ph2/(self.Nph_inj2*self.enebb2*self.N_wgt/self.N_time)
                else:
                    self.Nph_inj2 = 0
                    self.wph_inj2 = 1

                # TODO increase weight 
                self.wph_inj1 *= self.N_qdt 
                self.wph_inj2 *= self.N_qdt

            else: # inject uni weight photons (but lots of 'em)

                if do_print:
                    print("WARN: forcing photons to have w=1")
                self.wph_inj1 = 1.0
                self.wph_inj2 = 1.0

                #self.zeta_xinj1 = self.lum_ph1/( self.Nref['ph']*self.wph_inj1*self.enebb1*self.N_wgt/self.N_time*self.N_box)
                self.zeta_xinj1 = self.lum_ph1/( self.Nref['ph']*self.wph_inj1*self.enebb1*self.N_wgt/self.N_time)
                #self.zeta_xinj2 = self.lum_ph2/( self.Nref['ph']*self.wph_inj2*self.enebb2*self.N_wgt/self.N_time*self.N_box)
                self.zeta_xinj2 = self.lum_ph2/( self.Nref['ph']*self.wph_inj2*self.enebb2*self.N_wgt/self.N_time)

                self.Nph_inj1 = self.zeta_xinj1*self.Nref['ph']
                self.Nph_inj2 = self.zeta_xinj2*self.Nref['ph']


            #gplaw = ( pmax**(pslope+2) - pmin**(pslope+2) )/(pslope + 2)

            #integral_a^b (b^p x^(-p) (a x^p - x a^p))/(a b^p - b a^p) dx = 
            #(b^2 a^p - a^2 (p - 1) b^p + a (p - 2) b^(p + 1))/((p - 2) (a b^p - b a^p)) for 0<a<b

            def int_plaw(a, b, p):
                #(b^2 a^p - a^2 (p-1) b^p + a(p-2) b^(p + 1))/((p - 2) (a b^p - b a^p)) 
                return (b*b*a**p - a*a*(p-1)*b**p + a*(p-2)*b**(p+1))/( (p-2)*(a*b**p - b*a**p)) 
            gplaw = int_plaw(self.pmin, self.pmax, self.pslope)

            z = np.random.rand(10000000)
            v = ( self.pmin**(1+self.pslope) + z*( self.pmax**(1+self.pslope) - self.pmin**(1+self.pslope) ))**(1/(1+self.pslope))
            if do_print: print('init: mean brute force', np.mean(v), gplaw)
            gplaw = np.mean(v) # FIXME

            # v1
            #Nep_inj = zeta_pinj*Nref["e-"]
            #wep_inj = lum_ep/(Nep_inj*gplaw*N_inj) 

            #v2 # NOTE: assume wep_inj = 1 and calc zeta_pinj
            self.wep_inj = 1.0 
            # NOTE pair luminosity to electron luminosity so 0.5x
            self.zeta_pinj = 0.5*self.lum_ep/( self.Nref['e-']*self.wep_inj*gplaw*self.N_wgt/self.N_time )  
            self.Nep_inj  = self.zeta_pinj*self.Nref["e-"]   # injecting pair plasma so x2
            self.Nep_inj *= self.N_qdt # TODO increase by qed activity ratio

            Lp1_ergs = self.Nph_inj1*(self.wph_inj1/self.wsum0)*self.enebb1/self.dt
            Lp1_ergs *= 1/self.N_qdt # TODO normalize by weight since wph has N_qdt
            Lp1_ergs *= vol*self.nz0*me*c**2

            if do_print:
                print("init:-------------------- ph inj 1---------------------")
                print('init: lum_ph   :', self.lum_ph1)
                print('init: kTbb     :', self.kTbb1)
                print('init: enebb    :', self.enebb1)
                print('init: zeta_xinj:', self.zeta_xinj1, ' N:', self.Nph_inj1)
                print('init: wph_inj  :', self.wph_inj1, np.log10(self.wph_inj1) + 1e-100)
                print('init: Lp1      :', np.log10(Lp1_ergs + 1e-100), ' erg/s')

            Lp2_ergs = self.Nph_inj2*(self.wph_inj2/self.wsum0)*self.enebb2/self.dt
            Lp2_ergs *= 1/self.N_qdt # TODO normalize by weight since wph has N_qdt
            Lp2_ergs *= vol*self.nz0*me*c**2

            if do_print:
                print("init:-------------------- ph inj 2---------------------")
                print('init: lum_ph   :', self.lum_ph2)
                print('init: kTbb     :', self.kTbb2)
                print('init: enebb    :', self.enebb2)
                print('init: zeta_xinj:', self.zeta_xinj2, ' N:', self.Nph_inj2)
                print('init: wph_inj  :', self.wph_inj2, np.log10(self.wph_inj2 + 1e-100))
                print('init: Lp2      :', np.log10(Lp2_ergs + 1e-100), ' erg/s')



            # NOTE we inject electrons AND positrons w/ Nep_inj prtcls so 2x here
            Le_ergs = 2.0*self.Nep_inj*(self.wep_inj/self.wsum0)*gplaw/self.dt 
            Le_ergs *= 1/self.N_qdt # TODO normalize by weight since wep has N_qdt
            Le_ergs *= vol*self.nz0*me*c**2

            if do_print:
                print("init:-------------------- ep inj ----------------------")
                #print("init: pslope:   ", self.pslope)
                #print("init: pmin:     ", self.pmin)
                #print("init: pmax:     ", self.pmax)
                #print("init: pmean:    ", gplaw)
                #print("init: zeta_pinj:", self.zeta_pinj, ' N:', self.Nep_inj)
                #print("init: wep_inj:  ", self.wep_inj, np.log10(self.wep_inj + 1e-100))
                print("init: ell:      ", self.lum_ep)
                print("init: Le:       ", np.log10(Le_ergs + 1e-100), ' erg/s')
            
            # injected compactness from B field is sigma * \tau/t_A 
            # tau = self.nz0*sigT*H
            #self.lum_ant  = self.sigma*self.nz0*sigT*H*self.drive_ampl*self.drive_freq
            #self.lum_ant *= 4.0
            #Le_ant_ergs = self.lum_ant*me*c**3*H/sigT

            #if do_print:
            #    print("init:-------------------- ep ant ----------------------")
            #    print("init: ell:      ", self.lum_ant)
            #    print("init: Le:       ", np.log10(Le_ant_ergs + 1e-100), ' erg/s')

        else:

            # manually set some quantities to avoid problems later on
            self.N_wgt = 1
            self.N_box = 1
            self.N_time = 1
            self.npc_ref = 1
            self.N_qdt = self.qed_step


        # DONE
        if do_print:
            print("init:--------------------------------------------------")

        return

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

    # TODO debug
    return 0 # NOTE no injection in the beginning of the simulation
    
    #if xloc[0] > conf.height_atms + 1:
    #    return 0
    #if xloc[0] < conf.height_atms:
    #    return 0

    #if ispcs in [0,1]:
    #    return conf.ppc
    #if ispcs == 2:
    #    return conf.xpc

# Particle weight that is added to cell at location 'xloc'
def weigth_profile(xloc, ispcs, conf):
    return 1.0


