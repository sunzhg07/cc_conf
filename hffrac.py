#! /usr/bin/env python
from ordered import OrderedDict

class HFfracCalculation(OrderedDict):

    def __init__(this, config, interaction, space, core, nolabels, radius=False, com=False):
        OrderedDict.__init__(this)
        
        this.interaction = interaction
        this.modelspace = space
        this.core = core
        this.radius = radius
        this.com = com
        this.nolabels = nolabels

        this.config=config
        
        # Interation parameters
        format = interaction.get_twobody_format()
        this['nn_input_io_format'] = format
        if format == "ascii" or format == "binary":
            this['nn_input_file_orbits'] = this.get_sp_file()
            this['nn_input_file_onebody'] = this.get_kinetic_file()
            this['nn_input_file_twobody'] = this.get_interaction_file()
        elif format == "hdf5":
            this['nn_input_file_hdf5'] = this.get_twobody_hdf5_file()
            this['nn_input_group_hdf5'] = this.get_twobody_hdf5_group()
        else:
            raise Exception("Unknown twobody format: %s" % format)

        
        this['twobody_nmax'] = space.get_twobody_nmax()
        this['hbar_omega'] = space.get_hw()
        
        include_threebody = interaction.use_threebody()
        this['threebody'] = "no"
        
        if include_threebody:
            this['threebody'] = "yes"
            this['max_in_memory'] = config["max_in_memory"]
            format = interaction.get_threebody_format()
            this['threebody_interaction_format'] = format
            if format == "hdf5":
                this['file_threebody_hdf5'] = this.get_threebody_hdf5_file()
            elif format == "hdf5_split":
                this['file_threebody_hdf5'] = this.get_threebody_hdf5_split_file()
            this['threebody_topology'] = interaction.get_threebody_topology()
            this['threebody_nmax'] = this.modelspace.get_threebody_nmax()
            lecs = this.interaction.get_lecs()
            this['cD'] = lecs.cD
            this['cE'] = lecs.cE
            this['c1'] = lecs.c1
            this['c3'] = lecs.c3
            this['c4'] = lecs.c4



        this[''] = None
        this['occ_protons'] = core.get_nprotons_jscheme()
        this['occ_neutrons'] = core.get_nneutrons_jscheme()
        this['fracfill_file'] = core.get_fracfill_file()
        this['generate_mscheme_input'] = "yes"

        this['mass_nucleus'] = core.get_mass()

        this['com_switch'] = "0.00"
        if config['com_correction']: this['com_switch'] = "1.00"

        if radius: this['radius_beta'] = "0.001"
        else: this['radius_beta'] = "0.00"

        if com: this['com_beta'] = "0.001"
        else: this['com_beta'] = "0.00"
        
        this['j2_beta'] = "0.00"

        this['include_hf'] = "no"
        if config['hartree_fock']: this['include_hf'] = "yes"


        this['include_nat'] = "no"
        if config['natural_orbitals']: this['include_nat'] = "yes"
        if this['include_nat'] == "yes":        
            this['nat_modelspace_in_nmax'] = space.get_nat_twobody_nmax()

        
      
        this['precalc_hf'] = "no"
        if config['precalc_hf']: this['precalc_hf'] = "yes"

        this['precalc_density'] = "no"
        if config['precalc_density']: this['precalc_density'] = "yes"
        
        this['use_broyden'] = "yes"+"\n"
        
        
        mscheme = False
        if this['occ_protons'] + this['occ_neutrons'] == this['mass_nucleus']: mscheme = True
        
        
        this['fractional_filling'] = core.get_fractional_filling()
        if this['fractional_filling'] == "yes":
            this['occ_protons_fermi'] = core.get_nprotons_fermi()
            this['occ_neutrons_fermi'] = core.get_nneutrons_fermi()

        
        # check if we are using m-scheme code:
        #if (this['occ_protons'] + this['occ_neutrons'] == this['mass_nucleus']):
        
        this[''] = None
        this['prolate_filling'] = core.get_prolate_filling()
#        print this['prolate_filling']
        this['constrained'] = "no"
        if config['constrained_hf']: this['constrained'] = "yes"
        this['restricted_sort'] = "no"
        if config['restricted_sort']: this['restricted_sort'] = "yes"
        this['perform_vap'] = "no"
        this['perform_pav'] = "no"
        
        if this['constrained'] == "yes":
            this['q2_value'] = "0.00"
            this['alm_parameter'] = "0.00"

#        if this['fractional_filling'] == "yes":
#            this['no2b_output_file_zerobody'] = this.get_no0b_file()
#            if not mscheme:
#                this['no2b_output_file_onebody'] = this.get_no1b_file()
#                this['no2b_output_file_twobody'] = this.get_no2b_file()
#            else:
#                this['nn_input_file_onebody'] = this.get_no1b_file()
#                this['nn_input_file_twobody'] = this.get_no2b_file()

        if this['fractional_filling'] == "yes":
            this['no2b_output_file_zerobody'] = this.get_no0b_file()
            this['no2b_output_file_onebody'] = this.get_no1b_file()
            this['no2b_output_file_twobody'] = this.get_no2b_file()
                    
            
        if config["continuum"]:
            this["include_continuum"] = "yes"
            this["number_continuum"] = config["continuum_meshsize"]

        output_format = config["hf_format"]
        
        this["hf_output_io_format"] = output_format
#        if output_format == "hdf5":
#            this["hf_output_file_hdf5"] = this.get_hdf5_outfile()
#            this["hf_output_group_hdf5"] = this.get_hdf5_outgroup()
#        else:

        if this['fractional_filling'] == "no" or (this['fractional_filling'] == "yes" and mscheme and core.get_prolate_filling() == "yes"):
            this['hf_output_file_orbits_initial'] = this.get_mscheme_init_outfile()
            this['hf_output_file_orbits'] = this.get_sp_outfile() 
            this['hf_output_file_onebody'] = this.get_kinetic_outfile()
            this['hf_output_file_coefficients'] = this.get_umat_outfile()
            this['hf_output_file_fock_matrix'] = this.get_fock_outfile()
#            this['hf_output_file_channels'] = this.get_ch_outfile()
#            this['hf_output_file_twobody_pppp'] = this.get_pppp_outfile()
#            this['hf_output_file_twobody_rest'] = this.get_int_outfile()
            this["hf_output_file_twobody_h5"] = this.get_hdf5_outfile()
        if (this['fractional_filling'] == "no" and core.get_prolate_filling() == "no") or (this['fractional_filling'] == "yes" and mscheme and core.get_prolate_filling() == "no"):
            this['hf_output_file_orbits_initial'] = this.get_mscheme_init_outfile_oblate()
            this['hf_output_file_orbits'] = this.get_sp_outfile_oblate() 
            this['hf_output_file_onebody'] = this.get_kinetic_outfile_oblate()
            this['hf_output_file_coefficients'] = this.get_umat_outfile_oblate()
            this['hf_output_file_fock_matrix'] = this.get_fock_outfile_oblate()
#            this['hf_output_file_channels'] = this.get_ch_outfile()
#            this['hf_output_file_twobody_pppp'] = this.get_pppp_outfile()
#            this['hf_output_file_twobody_rest'] = this.get_int_outfile()
            this["hf_output_file_twobody_h5"] = this.get_hdf5_outfile_oblate()

        if config["continuum"]:
            this['hf_output_file_continuum_mesh'] = "cmesh_%s.dat" % this.get_outfile_tail()
            this['hf_output_file_continuum_basis'] = "cbasis_%s.dat" % this.get_outfile_tail()

    @staticmethod
    def factory(config, interactions, modelspaces, core, radius = False, com = False, nolabels=False):
        res = []
        for interaction in interactions:
            for space in modelspaces:
                calc = HFfracCalculation(config, interaction, space, core, nolabels)
                res.append(calc)
                if radius: res.append(HFfracCalculation(config, interaction, space, core, nolabels, radius=radius))
                if com: res.append(HFfracCalculation(config, interaction, space, core, nolabels, com=com))

        return res
    def get_twobody_hdf5_file(this):
        s = "%s_%s.h5" % (this.interaction.get_twobody_prefix(), this.get_file_tail())
        return s

    def get_twobody_hdf5_group(this): return this.get_file_tail()

    def get_threebody_hdf5_file(this):
        s = "%s_%s.h5" % (this.interaction.get_threebody_prefix(), this.get_threebody_file_tail())
        return s

    def get_threebody_hdf5_split_file(this):
        s = "%s_%s_split.h5" % (this.interaction.get_threebody_prefix(), this.get_threebody_file_tail())
        return s

    def get_interaction_file(this):
        s = "%s_%s.dat" % (this.interaction.get_twobody_prefix(), this.get_file_tail())
        return s

    def get_no0b_file(this):
        return "no0b_%s.dat" % (this.get_no_outfile_tail())

    def get_no1b_file(this):
        s = "no1b_%s.dat" % (this.get_no_outfile_tail())
        return s

    def get_no2b_file(this):
        s = "no2b_%s.dat" % (this.get_no_outfile_tail())
        return s

    def get_kinetic_file(this):
        s = "kinetic_%s.dat" % this.get_file_tail()
        return s

    def get_sp_file(this):
        s = "sp_energy_%s.dat" % this.get_file_tail()
        return s

    def get_sp_outfile(this):
        s = "sp_energy_%s.dat" % this.get_outfile_tail()
        return s

    def get_mscheme_init_outfile(this):
        s = "mscheme_init_%s.dat" % this.get_outfile_tail()
        return s

    def get_hdf5_outfile(this):
        return "vno2b_%s.h5" % (
            this.get_outfile_tail())

    
    def get_cc_gs_outfile(this,config):
        return "t2_l2_%s_%s.h5" % (config["cc_approx"],
            this.get_outfile_tail())

    def get_cc_symproj_outfile(this,config):
        return "%s_symproj_ode_%s.h5" % (config["cc_approx"],
            this.get_outfile_tail())

    def get_sp_outfile_oblate(this):
        s = "sp_energy_%s_oblate.dat" % this.get_outfile_tail()
        return s

    def get_mscheme_init_outfile_oblate(this):
        s = "mscheme_init_%s_oblate.dat" % this.get_outfile_tail()
        return s

    def get_hdf5_outfile_oblate(this):
        return "vno2b_%s_oblate.h5" % (
            this.get_outfile_tail())

    
    def get_cc_gs_outfile_oblate(this,config):
        return "t2_l2_%s_%s_oblate.h5" % (config["cc_approx"],
            this.get_outfile_tail())

    def get_cc_symproj_outfile_oblate(this,config):
        return "%s_symproj_ode_%s_oblate.h5" % (config["cc_approx"],
            this.get_outfile_tail())

    def get_kinetic_outfile_oblate(this):
        return "kinetic_%s_oblate.dat" % this.get_outfile_tail()

    def get_umat_outfile_oblate(this):
        return "umat_%s_oblate.dat" % this.get_outfile_tail()

    def get_fock_outfile_oblate(this):
        return "fock_%s_oblate.dat" % this.get_outfile_tail()

#    def get_hdf5_outfile(this):
#        s = "hf_%s.h5" % this.get_outfile_tail()
#        return s

    def get_hdf5_outgroup(this): return this.get_outfile_tail()

    def get_kinetic_outfile(this):
        return "kinetic_%s.dat" % this.get_outfile_tail()

    def get_umat_outfile(this):
        return "umat_%s.dat" % this.get_outfile_tail()

    def get_fock_outfile(this):
        return "fock_%s.dat" % this.get_outfile_tail()

    def get_ch_outfile(this):
        return "channels_%s.dat" % this.get_outfile_tail()

    def get_pppp_outfile(this):
        return "%s_pppp_%s.dat" % (
            this.interaction.get_twobody_prefix(),
            this.get_outfile_tail())

    def get_int_outfile(this):
        return "%s_rest_%s.dat" % (
            this.interaction.get_twobody_prefix(),
            this.get_outfile_tail())

    def get_outfile_tail(this):
        infix = ""
        if not this.nolabels: infix = "%s_" %  this.interaction.label

        
        if this.config['natural_orbitals']:
            calc = "NAT"
        else:
            calc = "HF"
            
        if this.config['constrained_hf']:
            calc = calc +"_Q"

#        print core.get_prolate_filling()
#        if not this.config['prolate_filling']: print calc
#            calc = calc +"_oblate"
#            print calc
#        this.config['nat_modelspace_in_nmax']
#        print this.config['prolate_filling']

        n3max = ""
        if this.interaction.use_threebody():
            lecs = this.interaction.get_lecs()
            #n3max = "E%02d_cD%g_cE%g"% (this.modelspace.get_threebody_nmax(), lecs.cD, lecs.cE)
            #n3max = "_N3N%02d"% (this.modelspace.get_threebody_nmax())
            n3max = "E%02d"% (this.modelspace.get_threebody_nmax())
        if this.com: infix = "cm1_"
        if this.radius: infix = "ra1_"
        
        if this.config['natural_orbitals']:
            s = "%sNnat%02d_N%02d%s_hw%g_%s_%s" % (
                infix,
                this.modelspace.nat_nmax,
                this.modelspace.nmax,
                n3max,
                this.modelspace.hw,
                this.core.label,calc)
        else:
             s = "%sN%02d%s_hw%g_%s_%s" % (
                infix,
                this.modelspace.nmax,
                n3max,
                this.modelspace.hw,
                this.core.label,calc)
       

            #,
#                    this.core.mass_correction)
        return s
    def get_file_tail(this):
        s = "N%02d_hw%g" % (this.modelspace.nmax, this.modelspace.hw)
        return s

    def get_threebody_file_tail(this):
        e3max = int(this.modelspace.get_threebody_nmax())
        nmax = int(this.modelspace.nmax)
        if nmax == e3max: return this.get_file_tail()
      #  s = "N%02dE%02d_hw%g" % (nmax, e3max, this.modelspace.hw)
        s = "N14E%02d_hw%g" % (e3max, this.modelspace.hw)
        return s

    def get_no_outfile_tail(this):
        infix = ""
        if not this.nolabels: infix = "%s_" %  this.interaction.label
        
        calc = "HF"
        n3max = ""
        if this.interaction.use_threebody():
            lecs = this.interaction.get_lecs()
            #n3max = "E%02d_cD%g_cE%g"% (this.modelspace.get_threebody_nmax(), lecs.cD, lecs.cE)
            #n3max = "_N3N%02d"% (this.modelspace.get_threebody_nmax())
            n3max = "E%02d"% (this.modelspace.get_threebody_nmax())
        if this.com: infix = "cm1_"
        if this.radius: infix = "ra1_"

        s = "%sN%02d%s_hw%g_%s_%s" % (
                    infix,
                    this.modelspace.nmax,
                    n3max,
                    this.modelspace.hw,
                    this.core.label,calc)#,
#                    this.core.mass_correction)
        return s
