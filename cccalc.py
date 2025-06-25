#! /usr/bin/env python
from ordered import OrderedDict
from hfcalc import HFCalculation

class CCCalculation(HFCalculation):

    def __init__(this, config, interaction, space, core, nolabels, radius=False, com=False):
        OrderedDict.__init__(this)

        this.interaction = interaction
        this.modelspace = space
        this.core = core
        this.radius = radius
        this.com = com
        this.nolabels = nolabels

        this.config=config

        input_format = config["hf_format"]
        
        this['# Specify single-particle data and model-space parameters'] = None
        this['twobody_nmax'] = space.get_twobody_nmax()
        this['hbar_omega'] = space.get_hw()
        this['occ_protons'] = core.get_nprotons()
        this['occ_neutrons'] = core.get_nneutrons()
        this['mass_nucleus'] = core.get_mass()
        this['q2_constrained'] =  "no"
        this['lambda_q2'] = "0.00" + "\n"

        this['# Input files:'] = None
        this["hf_format"] = "standard"
        #if input_format == "hdf5":
        #    this["hf_input_file_hdf5"] = this.get_hdf5_outfile()
        #    this["hf_input_group_hdf5"] = this.get_hdf5_outgroup()
        #else:
        this['hf_input_file_initial_orbs'] = this.get_mscheme_init_outfile()
        this['hf_input_file_orbits'] = this.get_sp_outfile()
        this['hf_input_file_onebody'] = this.get_kinetic_outfile()
        this['hf_input_file_coefficients'] = this.get_umat_outfile()
        this['hf_input_file_fock_matrix'] = this.get_fock_outfile()
        # this['hf_input_file_channels'] = this.get_ch_outfile()
        # this['hf_input_file_twobody_pppp'] = this.get_pppp_outfile()
        # this['hf_input_file_twobody_rest'] = this.get_int_outfile()
        this["hf_input_file_twobody_h5"] = this.get_hdf5_outfile() + "\n"
        
        this["# read initial basis"] = None
        this['hf_initial_orbs'] = "yes" + "\n"

        this['# hdf5 input: yes/no'] = None
        this['hdf5_input'] = "yes" + "\n"
        
        this['# Use pre calculated fock-matrix = yes/no'] = None
        this['pre_calculated_fock_matrix'] = "yes" + "\n"
        
        this['# number of iterations and tolerance for cc ground-state'] = None
        this['cc_iterations'] = config["cc_iter"]
        this['cc_tolerance'] = str(config["cc_tolerance"]) + "\n"

        this['# which type of g.s. calculation: ccsd/ccsdt1/ccsd(t)'] = None
        this['cc_approximation'] = str(config["cc_approx"]) + "\n"

        this['# e3max and occupation number cuts for triples:'] = None
        this['cc_e3max_cut'] = "1000"
        this['cc_sp_occ_cut'] = "0.00" + "\n"

        this['# input parameters for diis algorithm to improve ccsd convergence'] = None
        this['diis_subspace']  = "10"
        this['diis_step'] = "10"
        this['diis_mixing'] = "0.7"
        
        this['# output file for t- and lambda amplitudes'] = None
        this['groundstate_file'] = this.get_cc_gs_outfile(config) + "\n"
        
        this['# precalculated CC ground-state yes/no'] = None
        this['pre_calculated_groundstate'] = str(config["precalc_cc_gs"]) 
        this['pre_calculated_left_groundstate'] = str(config["precalc_cc_left_gs"]) + "\n"
        
        this['# type of cc-eom calculations: none/lit/0vbb/eomccsd(t)/transitions/pa_eomccsd(t)/groundstate'] = None
        this['eom_calc_type'] = str(config["cc_calculation"])
        this['eom_approximation'] = "ccsd"
        this['eom_e3max_cut'] = "0"
        this['eom_l3max_cut'] = "100"
        this['eom_sp_nocc_cut'] = "0.000"
        this['eom_vec_file'] = "eom_vec_ccsdt1_e3max200_inmed450_N08E16_hw12_K22A22.h5"
        this['eom_vector_stored'] = "no"
        this['eom_expectation_value'] = str(config["expectation_value"]) + "\n"
        
        this['number_of_states'] = "1"
        this['arnoldi_iter'] = "150"
        this['arnoldi_tolerance'] = "1.D-6" + "\n"
        
        this['# speficy quantum numbers for equation-of-motion calculation and parameters for arnoldi algorithm'] = None
        this['j_eom']    =  "0"
        this['ipar_eom'] =  "0"
        this['tz_eom']   =  "0" + "\n"
        
        this['# Specify symmetry restoration options'] = None
        this['symmetry_restoration'] = str(config["symmetry_restoration"]) 
        this['linear_response'] = "yes"
        this['pre_calculated_projection_vectors']="no"

        this['symproj_approx'] = "ode_lr"
        this['all_angles'] = "yes"
        this['number_of_angles'] = "15"
        this['which_angle'] = "1"
        this['cc_symproj_file'] = this.get_cc_symproj_outfile(config) + "\n"
        this['symproj_w0w1w2'] = this.get_cc_w0w1w2_outfile(config) + "\n"
        
        #this['projected_angmom'] = "0"
        #this['t2_vap_file'] = "t2_vap_test_be8.h5"
        #this['vap'] = "no"
             

#        this['file_ccsd_data'] = this.get_ccsd_file()
    @staticmethod
    def factory(config, interactions, modelspaces, core, radius = False, com = False, nolabels=False):
        res = []
        for interaction in interactions:
            for space in modelspaces:
                calc = CCCalculation(config, interaction, space, core, nolabels)
                res.append(calc)
                if radius: res.append(CCCalculation(config, interaction, space, core, nolabels, radius=radius))
                if com: res.append(CCCalculation(config, interaction, space, core, nolabels, com=com))

        return res
    def get_ccsd_file(this):
        return "ccsd_%s.h5" % this.get_outfile_tail()

class LTCalculation(CCCalculation):

    def __init__(this, config, interaction, space, core, nolabels, radius=False, com=False):
        OrderedDict.__init__(this)

        this.interaction = interaction
        this.modelspace = space
        this.core = core
        this.radius = radius
        this.com = com
        this.nolabels = nolabels

        input_format = config["hf_format"]
        this["hf_input_io_format"] = input_format
        if input_format == "hdf5":
            this["hf_input_file_hdf5"] = this.get_hdf5_outfile()
            this["hf_input_group_hdf5"] = this.get_hdf5_outgroup()
        else:
            this['hf_input_file_orbits'] = this.get_sp_outfile()
            this['hf_input_file_onebody'] = this.get_kinetic_outfile()
            this['hf_input_file_coefficients'] = this.get_umat_outfile()
            this['hf_input_file_fock_matrix'] = this.get_fock_outfile()
            this['hf_input_file_channels'] = this.get_ch_outfile()
            this['hf_input_file_twobody_pppp'] = this.get_pppp_outfile()
            this['hf_input_file_twobody_rest'] = this.get_int_outfile()
        this['normalordered_input'] = config["normalordered_input"]

        this['file_ccsd_data'] = this.get_ccsd_file()
        this['ccsd_input'] = "yes"

        this['twobody_nmax'] = space.get_twobody_nmax()
        this['hbar_omega'] = space.get_hw()

        this['occ_protons'] = core.get_nprotons()
        this['occ_neutrons'] = core.get_nneutrons()
        this['mass_nucleus'] = core.get_mass()

        this['arnoldi_iter'] = config["eom_iter"]
        this['arnoldi_tolerance'] = config["eom_tolerance"]

        this['eom_calc_type'] = "lambda_triples"
        this['eom_truncation'] = "none"

    @staticmethod
    def factory(config, interactions, modelspaces, core, radius = False, com = False, nolabels=False):
        res = []
        for interaction in interactions:
            for space in modelspaces:
                calc = LTCalculation(config, interaction, space, core, nolabels)
                res.append(calc)
                if radius: res.append(LTCalculation(config, interaction, space, core, nolabels, radius=radius))
                if com: res.append(LTCalculation(config, interaction, space, core, nolabels, com=com))

        return res
class DMCalculation(CCCalculation):

    def __init__(this, config, interaction, space, core, nolabels, radius=False, com=False):
        OrderedDict.__init__(this)

        this.interaction = interaction
        this.modelspace = space
        this.core = core
        this.radius = radius
        this.com = com
        this.nolabels = nolabels

        input_format = config["hf_format"]
        this["hf_input_io_format"] = input_format
        if input_format == "hdf5":
            this["hf_input_file_hdf5"] = this.get_hdf5_outfile()
            this["hf_input_group_hdf5"] = this.get_hdf5_outgroup()
        else:
            this['hf_input_file_orbits'] = this.get_sp_outfile()
            this['hf_input_file_onebody'] = this.get_kinetic_outfile()
            this['hf_input_file_coefficients'] = this.get_umat_outfile()
            this['hf_input_file_fock_matrix'] = this.get_fock_outfile()
            this['hf_input_file_channels'] = this.get_ch_outfile()
            this['hf_input_file_twobody_pppp'] = this.get_pppp_outfile()
            this['hf_input_file_twobody_rest'] = this.get_int_outfile()
        if config['continuum']:
            this["include_continuum"] = "yes"
            this["number_continuum"] = config["continuum_meshsize"]
            this['hf_input_file_continuum_mesh'] = "cmesh_%s.dat" % this.get_outfile_tail()
            this['hf_input_file_continuum_basis'] = "cbasis_%s.dat" % this.get_outfile_tail()

        this['normalordered_input'] = config["normalordered_input"]

        this['file_ccsd_data'] = this.get_ccsd_file()
        this['ccsd_input'] = "yes"

        this['twobody_nmax'] = space.get_twobody_nmax()
        this['hbar_omega'] = space.get_hw()

        this['occ_protons'] = core.get_nprotons()
        this['occ_neutrons'] = core.get_nneutrons()
        this['mass_nucleus'] = core.get_mass()

        this['arnoldi_iter'] = config["eom_iter"]
        this['arnoldi_tolerance'] = config["eom_tolerance"]

        this['eom_calc_type'] = "density_matrix"
        this['density_matrix_file'] = this.get_density_matrix_file()
        this['eom_truncation'] = "none"

    def get_density_matrix_file(this):
        return "density_matrix_%s.dat" % this.get_outfile_tail()
    @staticmethod
    def factory(config, interactions, modelspaces, core, radius = False, com = False, nolabels = False):
        res = []
        for interaction in interactions:
            for space in modelspaces:
                calc = DMCalculation(config, interaction, space, core, nolabels)
                res.append(calc)
                if radius: res.append(DMCalculation(config, interaction, space, core, nolabels, radius=radius))
                if com: res.append(DMCalculation(config, interaction, space, core, nolabels, com=com))

        return res
