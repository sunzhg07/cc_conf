#! /usr/bin/env python

import json
from hfcalc import HFCalculation
from cccalc import CCCalculation, LTCalculation, DMCalculation
from eomcalc import EOMCalculation

default_values =    {
                            "scheduler_type"        :   "serial",
                            "num_mpi_groups"        : 1,
                            "com_correction"        : True,
                            "hartree_fock"          : True,
                            "continuum"             : False,
                            "continuum_meshsize"    : 50,
                            "normalordered_input"  :   "yes",
                            "compare_strategy"      :   "lowest",
                            "cc_iter"               :   1000,
                            "cc_tolerance"          :   1e-8,
                            "eom_iter"              :   150,
                            "eom_tolerance"         :   1e-8,
                            "max_in_memory"         :   10,
                            "hf_format"             :   "hdf5",
                            "overwrite_result"      :   "yes"
                    }
class EomState(object):

    def __init__(this, c, j, tz, pi, n):
        this.radius = False
        if c.has_key("calculate_radius"): this.radius = c['calculate_radius']
        this.check_com_contamination = False
        if c.has_key('check_com_contamination') :
            this.check_com_contamination = c['check_com_contamination']

        this.angular_momentum = j
        this.isospin = tz
        this.parity = pi
        this.number_of_states = n

    def __str__(this):
        s = "%d, tz = %s, pi = %s, n = %d" % (this.angular_momentum, this.isospin, this.parity, this.number_of_states)

        return s

    @staticmethod
    def parse(configtree):
        res = []
        for s in configtree:
            nstates = s['number_of_states']
            for tz in s['isospin']:
                for pi in s['parity']:
                    for j in s['angular_momentum']:
                        res.append(EomState(s, j, tz, pi, nstates))
        return res

class Interaction(object):

    def __init__(this, interaction, lecs):
        this.type = interaction['type']
        this.twobody = interaction['twobody']
        this.label = interaction["label"]
        this.threebody = {'include' : False}
        this.lecs = lecs
        if interaction.has_key('threebody'):
            this.threebody = interaction['threebody']

    def __eq__(this, other):
        res = True
        if this.type != other.type: res = False
        if this.twobody != other.twobody: res = False
        if this.threebody != other.threebody: res = False
        return res

    def __str__(this):
        s = "\t\tType: %s" % this.type
        s += "\n\t\tTwobody_format: %s" % this.twobody['format']
        s += "\n\t\tTwobody_prefix: %s" % this.twobody['input_prefix']
        s += "\n\t\tThreebody_include: %s" % this.use_threebody()
        if this.use_threebody():
            s += "\n\t\tThreebody_format: %s" % this.threebody['format']
            s += "\n\t\tThreebody_topology: %s" % this.threebody['topology']
        s += "\n\t\tThreebody_prefix: %s" % this.threebody['input_prefix']

        return s
    def use_threebody(this): return this.threebody['include']
    def get_lecs(this) : return this.lecs

    def get_twobody_prefix(this): return this.twobody['input_prefix']
    def get_twobody_format(this): return this.twobody['format']
    def get_threebody_format(this): return this.threebody['format']
    def get_threebody_prefix(this): return this.threebody['input_prefix']
    def get_threebody_topology(this): return this.threebody['topology']

    @staticmethod
    def parse(configtree):

        res = {}
        for interaction in configtree:
            label = interaction['label']
            if interaction.has_key('threebody'):
                cD = interaction['threebody']['cD']
                cE = interaction['threebody']['cE']
                c1 = interaction['threebody']['c1']
                c3 = interaction['threebody']['c3']
                c4 = interaction['threebody']['c4']
                lecs = Lec.parse(cD, cE, c1, c3, c4)
                res[label] = []
                for lec in lecs:
                    res[label].append(Interaction(interaction, lec))
            else: res[label] = [Interaction(interaction, None)]
        return res

class Lec(object):
    def __init__(this, cD, cE, c1, c3, c4):
        this.cD = cD
        this.cE = cE
        this.c1 = c1
        this.c3 = c3
        this.c4 = c4

    @staticmethod
    def parse(cD, cE, c1, c3, c4):
        res = []

        if type(cE) is not list: cE = [cE]
        if type(cD) is not list: cD = [cD]
        for e in cE:
            for d in cD:
                res.append(Lec( d, e, c1, c3, c4))
        return res


    def __str__(this):
        s = "\t\tcD = %g, cD = %g, c1 = %g, c3 = %g, c4 = %g" % (this.cD, this.cE, this.c1, this.c3, this.c4)
        return s

    def __eq__(this, other):
        if \
            this.cD == other.cD \
            and this.cE == other.cE \
            and this.c1 == other.c1 \
            and this.c3 == other.c3 \
            and this.c4 == other.c4: return True

        return False

class Modelspace(object):

    def __init__(this, nmax, hw, n3max):
        this.nmax = nmax
        this.hw = hw
        this.n3max = n3max

    def __str__(this):
        s = "\t\tnmax = %g, hw = %g, n3max = %g" % (this.nmax, this.hw, this.n3max)
        return s

    def get_threebody_nmax(this): return this.n3max
    def get_twobody_nmax(this): return this.nmax
    def get_hw(this) : return this.hw
    @staticmethod
    def parse_outer(nmax, hws, n3max=None):
        use_n3max = True
        if not n3max: use_n3max = False

        res = []
        for n in nmax:
            for hw in hws:
                if not use_n3max: n3max = [n]
                for n3 in n3max:
                    res.append(Modelspace(n, hw, n3))

        return res

    def __eq__(this, other):
        res = True
        if this.nmax != other.nmax: res = False
        if this.hw != other.hw: res = False
        if this.n3max != other.n3max: res = False

        return res

    @staticmethod
    def parse(configtree):

        res = {}
        for space in configtree:
            label = space['label']
            nmax = space['nmax']
            hw = space['hw']
            if space.has_key("n3max"):
                n3max = space['n3max']
                res[label] = Modelspace.parse_outer(nmax, hw, n3max)
            else: res[label] = Modelspace.parse_outer(nmax, hw)

        return res

class Core(object):
    def __init__(this, label, np, nn, mass):
        this.label = label
        this.proton_holes = np
        this.neutron_holes = nn
        this.mass_correction = mass

    def get_nprotons(this) : return this.proton_holes
    def get_nneutrons(this) : return this.neutron_holes
    def get_mass(this) : return this.mass_correction
    def __eq__(this, other):
        res = True
        if this.label != other.label: res = False
        if this.proton_holes != other.proton_holes: res = False
        if this.neutron_holes != other.neutron_holes: res = False
        if this.mass_correction != other.mass_correction: res = False

        return res

    def __str__(this):
        s = "%s, %s, %s, %s" % \
            (this.label, this.proton_holes, this.neutron_holes, this.mass_correction)
        return s

class Run(object):
    keys = [ "eom_type" ]

    optional_keys = [
        "calculate_radius",
        "check_com_contamination" ]

    def is_empty(this):
        if this.state == None: return True
        return False
    def __init__(this, config, core, state):
        
        this.config = config.fromkeys(this.keys)

        for key in this.keys:
            if config.has_key(key): this.config[key] = config[key]

        for key in this.optional_keys:
            val = False
            if config.has_key(key) : val = config[key]
            this.config[key] = val

        # Override if state contains radius or com
        if state: this.config["calculate_radius"] = state.radius
        if state: this.config["check_com_contamination"] = state.check_com_contamination
        
        this.core = core
        this.state = state
    
    def calc_radius(this): return this.config["calculate_radius"]
    def calc_com(this): return this.config["check_com_contamination"]
    def __str__(this):
        s = ""
        for key, val in this.config.iteritems():
            s += "\t\t%s: %s\n" % (key, val)
        s += "%s\n" % this.core
        s += "\t\tState: %s\n" % this.state

        return s
    @staticmethod
    def parse(configtree):
        res = []

        for r in configtree:
            for mass in r['mass_correction']:
                core = Core(    r['core'],
                                r['proton_holes'],
                                r['neutron_holes'],
                                mass )
                if r.has_key('states'):
                    states = EomState.parse(r['states'])
                    for s in states:
                        res.append(Run( r, core, s))
                else:
                    res.append(Run( r, core, None))

        return res

class Calculation(object):
    keys = [
        "density_matrix",
        "lambda_triples" ]
    optional_keys = [
        "normalordered_input",
        "com_correction",
        "hartree_fock",
        "continuum",
        "continuum_meshsize",
        "cc_iter",
        "cc_tolerance",
        "eom_iter",
        "compare_strategy",
        "eom_tolerance",
        "overwrite_result",
        "max_in_memory",
        "hf_format",
        "scheduler_type" ]

    def generate_hf_calculations(this, nolabels):
        res = []

        for core in this.cores:
            radius = this.get_radius_flag(core)
            com = this.get_com_flag(core)
            res.extend(HFCalculation.factory(this.config, this.interactions,\
                this.modelspaces, core, radius, com, nolabels))
        return res

    def generate_cc_calculations(this, nolabels):
        res = []

        for core in this.cores:
            radius = this.get_radius_flag(core)
            com = this.get_com_flag(core)
            res.extend(CCCalculation.factory(this.config, this.interactions,\
                this.modelspaces, core, radius, com, nolabels))
        return res

    def generate_eom_calculations(this, nolabels):
        res = []

        if this.config['density_matrix']:
            for core in this.cores:
                radius = this.get_radius_flag(core)
                com = this.get_com_flag(core)
                res.extend(DMCalculation.factory(this.config, this.interactions,\
                    this.modelspaces, core, radius, com, nolabels))

        if this.config['lambda_triples']:
            for core in this.cores:
                radius = this.get_radius_flag(core)
                com = this.get_com_flag(core)
                res.extend(LTCalculation.factory(this.config, this.interactions,\
                    this.modelspaces, core, radius, com, nolabels))

        for run in this.runs:
            if run.is_empty(): continue
            res.extend(EOMCalculation.factory(this.config, this.interactions,\
                this.modelspaces, run, nolabels))

        return res

    def get_radius_flag(this, core):
        res = False
        if core in this.radius_flags: res = True
        return res

    def get_com_flag(this, core):
        res = False
        if core in this.com_flags: res = True
        return res

    def __str__(this):
        s = "Calculation: \n"

        for key, val in this.config.iteritems():
            s += "\t%s: %s\n" % (key, val)

        s += "\tRuns: \n"
        for r in this.runs: s += "%s\n" % r

        s += "\tInteractions: \n"
        for r in this.interactions: s += "%s\n" % r

        s += "\tModelspaces: \n"
        for r in this.modelspaces: s += "%s\n" % r

        return s
    def __init__(this, config, runs, ints, spaces, params):
        this.config = config.fromkeys(this.keys)

        for key in this.config.keys():
            this.config[key] = config[key]

        for key in this.optional_keys:
            val = default_values[key]
            if params.has_key(key) :
                val = params[key]
            else:
                params[key] = val
            if config.has_key(key) : val = config[key]
            this.config[key] = val

        this.cores = []
        this.radius_flags = []
        this.com_flags = []

        this.runs = runs
        this.cores = this.get_unique_cores()
        this.update_core_flags()

        this.interactions = []
        for label in config['interaction']:
            this.interactions.extend(ints[label])

        this.modelspaces = []
        for label in config['modelspace']:
            this.modelspaces.extend(spaces[label])

    def get_unique_cores(this):

        cores = []
        for r in this.runs:
            core = r.core
            if core not in cores: cores.append(core)

        return cores

    def update_core_flags(this):
        for r in this.runs:
            core = r.core
            if r.calc_radius() and core not in this.radius_flags:
                this.radius_flags.append(core)
            if r.calc_com() and core not in this.com_flags:
                this.com_flags.append(core)

    @staticmethod
    def parse(configtree, interactions, spaces, params):
        res = []
            
        for c in configtree:
            res.append(Calculation(c, Run.parse(c['run']), interactions, spaces, params))
        return res


class Parameters(dict):
    optional_keys = [
        "scheduler_type",
        "num_mpi_groups"
    ]
    def __init__(this, configtree):
        dict.__init__(this, configtree.copy())
        for key in this.optional_keys:
            val = default_values[key]
            if not this.has_key(key) :
                this[key] = val

class Configuration(object):
    def __init__(this, calcs, ints, spaces, params):
        this.calculations = calcs
        this.interactions = ints
        this.modelspaces = spaces
        this.parameters = params

    def __str__(this):
        s = ""
        for c in this.calculations: s += "%s\n" % c
        return s

    def get_num_groups(this) : return this.parameters['num_mpi_groups']
    def get_scheduler_type(this) : return this.parameters['scheduler_type']

    def generate_hf_calculations(this, nolabels):
        res = []

        calcs = []
        for c in this.calculations:
            calcs.extend(c.generate_hf_calculations(nolabels))
        for c in calcs:
            if c not in res: res.append(c)
        return res

    def generate_cc_calculations(this, nolabels):
        res = []

        calcs = []
        for c in this.calculations:
            calcs.extend(c.generate_cc_calculations(nolabels))
        for c in calcs:
            if c not in res: res.append(c)
        return res

    def generate_eom_calculations(this, nolabels):
        res = []

        calcs = []
        for c in this.calculations:
            calcs.extend(c.generate_eom_calculations(nolabels))
        for c in calcs:
            if c not in res: res.append(c)
        return res

    def write_config_files(this, filename, nolabels=False):
        hf = this.generate_hf_calculations(nolabels)
        cc = this.generate_cc_calculations(nolabels)
        eom = this.generate_eom_calculations(nolabels)
        print "Number of HF configs: %d" % len(hf)
        print "Number of CC configs: %d" % len(cc)
        print "Number of eom configs: %d" % len(eom)

        fh_main = open(filename, 'w')
        fh_main.write("scheduler_type = %s\n" % this.get_scheduler_type())
        fh_main.write("num_groups = %d\n" % this.get_num_groups())
        fh_main.write("num_hf_runs = %d\n" % len(hf))
        fh_main.write("num_cc_runs = %d\n" % len(cc))
        fh_main.write("num_eom_runs = %d\n" % len(eom))

        for i, c in zip(range(1, len(hf)+1), hf):
            filename = "hf_conf_%04d.conf" % i
            fh_main.write("hf_config_filename_%04d = %s\n" % (i, filename))
            fh_conf = open(filename, 'w')
            for k,v in c.items():
                fh_conf.write("%s = %s\n" % (k,v))
            fh_conf.close()

        for i, c in zip(range(1, len(cc)+1), cc):
            filename = "ccsd_conf_%04d.conf" % i
            fh_main.write("cc_config_filename_%04d = %s\n" % (i, filename))
            fh_conf = open(filename, 'w')
            for k,v in c.items():
                fh_conf.write("%s = %s\n" % (k,v))
            fh_conf.close()

        for i, c in zip(range(1, len(eom)+1), eom):
            filename = "eom_conf_%04d.conf" % i
            fh_main.write("eom_config_filename_%04d = %s\n" % (i, filename))
            fh_conf = open(filename, 'w')
            for k,v in c.items():
                fh_conf.write("%s = %s\n" % (k,v))
            fh_conf.close()

        fh_main.close()


class ConfigFactory(object):

    @staticmethod
    def factory(filename):

        fh = open(filename, 'r')
        lines = fh.readlines(); fh.close()
        res = []
        for line in lines:
            if not line.strip().startswith("#"): res.append(line)

        configuration = json.loads("\n".join(res))

        content = configuration['configuration']
        # Parse parameters
        if content.has_key("parameters"): parameters = Parameters(content["parameters"])
        else: parameters = Parameters({})

        # Parse modelspace
        modelspaces = Modelspace.parse(content["modelspace"])

        # Parse interactions
        interactions = Interaction.parse(content["interaction"])

        # Parse calculations
        calculations = Calculation.parse(content["calculation"], interactions,
            modelspaces, parameters)

        # Create configuration object
        config = Configuration(calculations, interactions, modelspaces, parameters)

        return config



if __name__ == "__main__":
    import sys

    usage = """Usage: %s configfile""" % sys.argv[0]

    try:
        configfile = sys.argv[1].strip()
    except:
        print(usage)

    config = ConfigFactory.factory(configfile)
    config.write_config_files("main_conf.ini")
