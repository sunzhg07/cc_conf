#! /usr/bin/env python
from ordered import OrderedDict
from cccalc import CCCalculation

class EOMCalculation(CCCalculation):
    eom_calcs = {
        "3p1h" : "2paeomccsd",
        "3p1h_left" : "2paeomccsd_left",
        "3p1h_both" : "2paeomccsd_both",
        "2p2h" : "eomccsd_ex",
        "2p2h_left" : "eomccsd_ex_left",
        "2p2h_both" : "eomccsd_ex_both",
        "2p1h" : "paeomccsd_new",
        "3p2h" : "paeomccsd_new",
        "1p2h" : "preomccsd" }

    def __init__(this, config, interaction, space, run, nolabels, radius=False, com=False):
        OrderedDict.__init__(this)

        this.interaction = interaction
        this.modelspace = space
        this.run = run
        this.state = run.state
        this.core = run.core
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

        this['occ_protons'] = this.core.get_nprotons()
        this['occ_neutrons'] = this.core.get_nneutrons()
        this['mass_nucleus'] = this.core.get_mass()

        this['arnoldi_iter'] = config["eom_iter"]
        this['arnoldi_tolerance'] = config["eom_tolerance"]

        this['result_prefix'] = this.get_prefix()
        this['eom_calc_type'] = this.eom_calcs[run.config['eom_type']]
        this['eom_truncation'] = run.config['eom_type'].split("_")[0]
        this['compare_strategy'] = config['compare_strategy']
        this['overwrite_result'] = config['overwrite_result']
        this['number_eom_states'] = 1
        j = this.state.angular_momentum
        tz = this.get_isospin()
        pi = this.get_parity()
        this["state_01"] = "%d, %d, %d, %d" % (j, pi, tz, this.state.number_of_states)

    def get_isospin(this):
        isomarker = this.state.isospin
        isospin = -99
        if this.eom_calcs[this.run.config['eom_type']] == "2paeomccsd":
            if isomarker == "nn" : isospin = 2
            if isomarker == "pn" : isospin = 0
            if isomarker == "pp" : isospin = -2
        elif this.eom_calcs[this.run.config['eom_type']] == "2paeomccsd_left":
            if isomarker == "nn" : isospin = 2
            if isomarker == "pn" : isospin = 0
            if isomarker == "pp" : isospin = -2
        elif this.eom_calcs[this.run.config['eom_type']] == "2paeomccsd_both":
            if isomarker == "nn" : isospin = 2
            if isomarker == "pn" : isospin = 0
            if isomarker == "pp" : isospin = -2
        elif this.eom_calcs[this.run.config['eom_type']] == "paeomccsd_new":
            if isomarker == "p" : isospin = -1
            if isomarker == "n" : isospin = 1
        elif this.eom_calcs[this.run.config['eom_type']] == "preomccsd":
            if isomarker == "p" : isospin = -1
            if isomarker == "n" : isospin = 1
        elif this.eom_calcs[this.run.config['eom_type']] == "eomccsd_ex":
            if isomarker == "pp" : isospin = -2
            if isomarker == "p" : isospin = -1
            if isomarker == "0" : isospin = 0
            if isomarker == "n" : isospin = 1
            if isomarker == "nn" : isospin = 2
        elif this.eom_calcs[this.run.config['eom_type']] == "eomccsd_ex_left":
            if isomarker == "pp" : isospin = -2
            if isomarker == "p" : isospin = -1
            if isomarker == "0" : isospin = 0
            if isomarker == "n" : isospin = 1
            if isomarker == "nn" : isospin = 2
        elif this.eom_calcs[this.run.config['eom_type']] == "eomccsd_ex_both":
            if isomarker == "pp" : isospin = -2
            if isomarker == "p" : isospin = -1
            if isomarker == "0" : isospin = 0
            if isomarker == "n" : isospin = 1
            if isomarker == "nn" : isospin = 2

        return isospin

    def get_prefix(this):
        return 'eomresult_%s' % this.get_outfile_tail()
    def get_parity(this):
        pmark = this.state.parity
        if pmark == "+": return 0
        elif pmark == "-" : return 1
        raise Exception("Unknown parity: %s" % pmark)

    @staticmethod
    def factory(config, interactions, modelspaces, run, nolabels):
        res = []
        radius = run.calc_radius()
        com = run.calc_com()
        for interaction in interactions:
            for space in modelspaces:
                calc = EOMCalculation(config, interaction, space, run, nolabels)
                res.append(calc)
                if radius: res.append(EOMCalculation(config, interaction, space, run, nolabels, radius=radius))
                if com: res.append(EOMCalculation(config, interaction, space, run, nolabels, com=com))

        return res
