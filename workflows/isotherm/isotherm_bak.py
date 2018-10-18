from aiida.orm import CalculationFactory, DataFactory
from aiida.orm import load_node
from aiida.orm.code import Code
from aiida.orm.data.base import Bool, Str
from aiida.work import workfunction as wf
from aiida.work.run import run, submit
from aiida.work.workchain import WorkChain, ToContext, if_, while_, Outputs

# subworkflows
from charges import DdecChargesWorkChain

import matplotlib.pyplot as plt

# calculations
Cp2kCalculation = CalculationFactory('cp2k')
DdecCalculation = CalculationFactory('ddec')
RaspaCalculation = CalculationFactory('raspa')
ZeoppCalculation = CalculationFactory('zeopp.network')

# parameters
ArrayData = DataFactory('array')
CifData = DataFactory('cif')
NetworkParameters = DataFactory('zeopp.parameters')
ParameterData = DataFactory('parameter')
StructureData = DataFactory('structure')

@wf
def from_cif_to_structuredata(cif):
    """Helper function that converts CifData object into StructureData"""
    a = cif.get_ase()
    return StructureData(ase=a)
    
    
default_options = {
        "resources": {
            "num_machines": 1,
            "tot_num_mpiprocs": 1,
            "num_mpiprocs_per_machine": 1,
            },
        "max_wallclock_seconds": 30 * 60,
        "withmpi": False,
        }


class Isotherm(WorkChain):
    """
    Workchain that for a given matherial will compute an isotherm of a certain gaz adsorption
    """

    @classmethod
    def define(cls, spec):
        """
        This is the most important method of a Workchain, that defines the
        inputs that it takes, the logic of the execution and the outputs
        that are generated in the process 
        """
        super(Isotherm, cls).define(spec)
        
        # structure, adsorbant, pressures
        spec.input('structure', valid_type=CifData)
        spec.input("probe_molecule", valid_type=ParameterData)
        spec.input("pressures", valid_type=ArrayData)

        # cp2k
        spec.input('cp2k_code', valid_type=Code)
        spec.input('cp2k_parameters', valid_type=ParameterData, required=False,
                default=ParameterData(dict={}))
        spec.input("cp2k_options", valid_type=ParameterData,
                default=None, required=False)
        spec.input('cp2k_parent_folder', valid_type=RemoteData,
                default=None, required=False)

        # ddec
        spec.input('ddec_code', valid_type=Code)
        spec.input('ddec_parameters', valid_type=ParameterData, required=False,
                default=ParameterData(dict=default_ddec_params))
        spec.input("ddec_options", valid_type=ParameterData,
                default=None, required=False)

        # zeo++ 
        spec.input("zeopp_code", valid_type=Code)
        spec.input("zeopp_options", valid_type=ParameterData,
                default=ParameterData(dict=default_options))

        # raspa
        spec.input("raspa_code", valid_type=Code)
        spec.input("raspa_parameters", valid_type=ParameterData, required=True)
        spec.input("raspa_options", valid_type=ParameterData,
                default=ParameterData(dict=default_options))

        spec.input("_interactive", valid_type=bool, required=False, default=False)
        spec.input("_usecharges", valid_type=bool, required=False, default=True)
       
        # workflow
        spec.outline(
            cls.init,
            if_(cls.should_use_charges)(
                cls.run_point_charges,
            ),
            cls.run_geom_zeopp,
            cls.parse_geom_zeopp,
            if_(cls.should_run_block_zeopp)(
                cls.run_block_zeopp,
                cls.parse_block_zeopp,
            ),
            cls.run_henry_raspa,
            while_(cls.should_run_loading_raspa)(
                cls.run_loading_raspa,
                cls.parse_loading_raspa,
            ),
            cls.return_result,
        )
        
        # Here we define the output the Workchain will generate and
        # return. Dynamic output allows a variety of AiiDA data nodes
        # to be returned
        spec.dynamic_output()

    def init(self):
        """
        Initialize variables and the pressures we want to compute
        """
        self.ctx.current_p_indx = 0
        self.ctx.pressures = self.inputs.pressures.get_array("pressures")
        self.ctx.result = []
        
        self.ctx.processed_structure = self.inputs.structure

        self.ctx.raspa_parameters = self.inputs.raspa_parameters.get_dict()
        if self.inputs._usecharges:
            self.ctx.raspa_parameters['GeneralSettings']['UseChargesFromCIFFile'] = "yes"

        if self.inputs._interactive == True:
            self.fig, self.ax = plt.subplots(1,1)
            self.ax.set_xlabel(u"Pressure [bar]")
            self.ax.set_ylabel(u"Loading average [molecules/unit cell]")

    def should_use_charges(self):
        """
        Whether it is needed to employ charges
        """
        return self.inputs._usecharges

    def run_cp2k_charge_density(self):
        """Compute the charge-density of a structure that can be later
        used for extracting ddec point charges."""

        inputs = {
            'structure'       : from_cif_to_structuredata(self.ctx.structure),
            'cp2k_code'       : self.inputs.cp2k_code,
            'cp2k_parameters' : self.inputs.cp2k_parameters,
            'cp2k_options'    : self.inputs.cp2k_options,
            'ddec_code'       : self.inputs.ddec_code,
            'ddec_parameters' : self.inputs.ddec_parameters,
            'ddec_options'    : self.inputs.ddec_options,
            '_label'          : "run_point_charges",
        }

        # Create the calculation process and launch it
        future  = submit(DdecChargesWorkChain, **inputs)
        self.report("pk: {} | Running cDdecChargesWorkChain to compute the "
                "point charges")
        return ToContext(point_charges_calc=Outputs(future))

    def run_geom_zeopp(self):
        """This is the main function that will perform a raspa
        calculation for the current pressure"""
        
        
        # network parameters
        sigma = self.inputs.probe_molecule.get_dict()['sigma']
        NetworkParameters = DataFactory('zeopp.parameters')
        params = NetworkParameters({
            'ha': True,
            'res': True,
            'sa': [sigma, sigma, 100000],
            'volpo': [sigma, sigma, 100000],
        }

        # Create the input dictionary
        inputs = {
            'code'       : self.inputs.zeopp_code,
            'structure'  : self.ctx.structure,
            'parameters' : params,
            '_options'   : self.inputs.zeopp_options,
            '_label'     : "run_geom_zeopp",
        }

        # Create the calculation process and launch it
        future  = submit(ZeoppCaclculation.process(), **inputs)
        self.report("pk: {} | Running geometry analysis with zeo++".format(future.pid))

        return ToContext(zeopp=Outputs(future))

    def parse_geom_zeopp(self):
        """
        Extract the pressure and loading average of the last completed raspa calculation
        """
        self.ctx.raspa_parameters['GeneralSettings']['HeliumVoidFraction'] = \
        self.ctx.zeopp["pore_volume_volpo"].dict.POAV_Volume_fraction

    def should_run_block_zeopp(self):
        """If the pore non-accessible volume is 0 - there is no need to run"""
        return self.ctx.zeopp["pore_volume_volpo"].dict.PONAV_Volume_fraction <= 0.001
        
    def run_block_zeopp(self):
        """This is the main function that will perform a raspa
        calculation for the current pressure."""

        # Create the input dictionary
        sigma = self.inputs.probe_molecule.get_dict()['sigma']
        params = NetworkParameters({
            'ha':True,
            'block': [sigma, 200],
        })

        inputs = {
            'code'       : self.inputs.zeopp_code,
            'structure'  : self.ctx.structure,
            'parameters' : params,
            '_options'   : self.inputs._options,
            '_label'     : "run_block_zeopp",

        }

        # Create the calculation process and launch it
        future  = submit(ZeoppCaclulation.process(), **inputs)
        self.report("pk: {} | Running zeo++ block volume calculation".format(future.pid))
        self.ctx.block_pk = future.pid
        return ToContext(zeopp=Outputs(future))
    
    def parse_block_zeopp(self):
        """Extract the pressure and loading average of the last completed raspa calculation"""

        self.ctx.raspa_parameters['Component'][0]['BlockPockets'] = True
        self.ctx.raspa_parameters['Component'][0]['BlockPocketsPk'] = self.ctx.block_pk

    def should_run_loading_raspa(self):
        """
        This is the main condition of the while loop, as defined
        in the outline of the Workchain. We only run another
        raspa calculation if the current iteration is smaller than
        the total number of pressures we want to compute
        """
        return self.ctx.cuttent_p_index < len(self.ctx.pressures)

    def run_loading_raspa(self):
        """This is the main function that will perform a raspa
        calculation for the current pressure"""

        pressure = self.ctx.pressures[self.ctx.p]
        self.ctx.raspa_parameters['GeneralSettings']['ExternalPressure'] = pressure
        if self.ctx.prev_pk is not None:
            self.ctx.raspa_parameters['GeneralSettings']['RestartFile'] = True
            self.ctx.raspa_parameters['GeneralSettings']['RestartFilePk'] = self.ctx.prev_pk

        # Create the input dictionary
        inputs = {
            'code'       : self.inputs.raspa_code,
            'structure'  : self.ctx.structure,
            'parameters' : ParameterData(dict=self.ctx.raspa_parameters),
            '_options'   : self.inputs._options,
            '_label'     : "run_loading_raspa",
        }

        # Create the calculation process and launch it
        process = RaspaCalculation.process()
        future  = submit(process, **inputs)
        self.report("pk: {} | Running raspa for the pressure {} [bar]".format(future.pid, pressure/1e5))

        self.ctx.p += 1
        self.ctx.prev_pk = future.pid

        return ToContext(raspa=Outputs(future))

    def parse_loading_raspa(self):
        """
        Extract the pressure and loading average of the last completed raspa calculation
        """
        pressure = self.ctx.raspa_parameters['GeneralSettings']['ExternalPressure']/1e5
        loading_average = self.ctx.raspa["component_0"].dict.loading_absolute_average
        self.ctx.result.append((pressure, loading_average))        
        self.plot_data()

    def return_result(self):
        """
        Attach the results of the raspa calculation and the initial structure to the outputs
        """
        result = {
            "initial_structure": self.inputs.structure,
            "result": ParameterData(dict={"isotherm": self.ctx.result}),
        }

        for link_name, node in result.iteritems():
            self.out(link_name, node)

        self.report("Workchain <{}> completed successfully".format(self.calc.pk))
        return
    def run_henry_raspa(self):
        """
        This is the main function that will perform a raspa
        calculation for the current pressure
        """
        raspa_parameters = self.inputs.raspa_parameters.get_dict()
        raspa_parameters['GeneralSettings'].pop('ExternalPressure')
        for i, comp in enumerate(raspa_parameters['Component']):
            name = comp['MoleculeName']
            raspa_parameters['Component'][0] = {
                "MoleculeName"                     : name,
                "MoleculeDefinition"               : "TraPPE",
                "WidomProbability"                 : 1.0,
                "CreateNumberOfMolecules"          : 0,
            }
        # Create the input dictionary
        inputs = {
            'code'       : Code.get_from_string(self.inputs.raspa_codename.value),
            'structure'  : self.ctx.processed_structure,
            'parameters' : ParameterData(dict=raspa_parameters),
            '_options'   : self.inputs._options,
            '_label'     : "run_henry_raspa",
        }

        # Create the calculation process and launch it
        process = RaspaCalculation.process()
        future  = submit(process, **inputs)
        self.report("pk: {} | Running raspa for the Henry coefficients".format(future.pid))

        return 

    def plot_data(self, init=False):
        if self.inputs._interactive == False:
            return
        self.ax.plot(*zip(*self.ctx.result), marker='o', linestyle='--', color='r')
        self.fig.canvas.draw()



import ipywidgets as ipw

        
class IsothermSettings():
    """
    A class that contain jupyter widgets simplifying to setup Isotherm workflow
    """
    def __init__(self):
        self.layout = ipw.Layout(width="400px")
        self.style = {"description_width":"180px"}

    def settings_panel(self):
        self.cutoff = ipw.FloatText(
            value=14.0,
            step=0.2,
            description='Cutoff [A]:',
            disabled=False,
            layout=self.layout,
            style=self.style,
        )
        
        self.ff = ipw.Dropdown(
            options=(
                'GenericMOFs',
                'LSMO_DREIDING-TraPPE',
                'LSMO_UFF-TraPPE'
            ),
            value='GenericMOFs',
            description='Select forcefield:',
            layout=self.layout,
            style=self.style,
            )
        mol_opt = [
            ("Carbon dioxide", {'molecule':'CO2','sigma': 1.525}),
            ('methane', {'molecule':'methane','sigma': 1.86}),
            ("Nitrogen", {'molecule':'N2', 'sigma': 1.86})
        ]
        self.probe_molecule = ipw.Dropdown(
            options=mol_opt,
            #value=0,
            description='Select molecule:',
            layout=self.layout,
            style=self.style,
            )

        self.init_cycles = ipw.IntText(
                value=5000,
                step=1000,
                description = "Number of initialization cycles",
                disabled=False,
                layout=self.layout,
                style=self.style,
            )
        self.prod_cycles = ipw.IntText(
                value=5000,
                step=1000,
                description = "Number of production cycles",
                disabled=False,
                layout=self.layout,
                style=self.style,
            )

        return ipw.VBox([
            ipw.HBox([self.init_cycles, self.prod_cycles]),
            self.ff,
            self.probe_molecule,
            self.cutoff,
        ])

    def return_raspa_parameters (self, cif):
        from math import cos, sin, sqrt, pi
        import numpy as np
        deg2rad=pi/180.

        struct=cif.values.dictionary.itervalues().next()

        a = float(struct['_cell_length_a'])
        b = float(struct['_cell_length_b'])
        c = float(struct['_cell_length_c'])

        alpha = float(struct['_cell_angle_alpha'])*deg2rad
        beta  = float(struct['_cell_angle_beta'])*deg2rad
        gamma = float(struct['_cell_angle_gamma'])*deg2rad

        # this code makes sure that the structure is replicated enough times in x, y, z direction
        # in order to be compatible with the cutoff value
        
        # first step is computing cell parameters according to  https://en.wikipedia.org/wiki/Fractional_coordinates
        v = sqrt(1-cos(alpha)**2-cos(beta)**2-
                 cos(gamma)**2+2*cos(alpha)*cos(beta)*cos(gamma))
        cell=np.zeros((3,3))
        cell[0,:] = [a, 0, 0]
        cell[1,:] = [b*cos(gamma), b*sin(gamma),0]
        cell[2,:] = [c*cos(beta), c*(cos(alpha)-cos(beta)*cos(gamma))/(sin(gamma)),c*v
                   /sin(gamma)]
        cell=np.array(cell)

        # diagonalizing the cell matrix
        diag = np.diag(cell)
        # and computing nx, ny and nz
        nx, ny, nz = tuple(int(i) for i in np.ceil(self.cutoff.value/diag*2.))

        params ={
                "GeneralSettings":
                {
                "SimulationType"                   : "MonteCarlo",
                "NumberOfCycles"                   : self.prod_cycles.value,  
                "NumberOfInitializationCycles"     : self.init_cycles.value,   # 20000

                "PrintEvery"                       : 2000,

                "ChargeMethod"                     : "Ewald",
#                "UseChargesFromCIFFile"            : "yes",
                "CutOff"                           : 14.0,
                "Forcefield"                       : "{}".format(self.ff.value),
                "EwaldPrecision"                   : 1e-6,

                "Framework"                        : 0,
                "UnitCells"                        : "{} {} {}".format(nx, ny, nz),
                "HeliumVoidFraction"               : 0.0,

                "ExternalTemperature"              : 298.0,
                "ExternalPressure"                 : 58e4,
                },
                "Component":
                [{
                "MoleculeName"                     : "{}".format(self.probe_molecule.value['molecule']),
                "MoleculeDefinition"               : "TraPPE",
                "TranslationProbability"           : 0.5,
                "ReinsertionProbability"           : 0.5,
                "SwapProbability"                  : 1.0,
                "CreateNumberOfMolecules"          : 0,
                }],
                }
        return params
    
    def return_cp2k_parameters (self, structure, smearing=False):
        multiplicity = 1 + structure.get_ase().get_chemical_symbols().count('Cu')
        params = {
            'FORCE_EVAL': {
                'METHOD': 'Quickstep',
                'DFT': {
                    'MULTIPLICITY': multiplicity,
                    'UKS': True if multiplicity != 1 else False,
                    'CHARGE': 0,
                    'BASIS_SET_FILE_NAME': 'BASIS_MOLOPT',
                    'POTENTIAL_FILE_NAME': 'GTH_POTENTIALS',
                    'QS': {
                        'METHOD':'GPW',
                    },
                    'POISSON': {
                        'PERIODIC': 'XYZ',
                    },
                    'MGRID': {
                        'CUTOFF':     400,
                        'NGRIDS':       4,
                        'REL_CUTOFF':  50,
                    },
                    'SCF':{
                        'SCF_GUESS': 'ATOMIC',
                        'EPS_SCF': 1.0e-6,
                        'MAX_SCF': 50,
                        'MAX_ITER_LUMO': 10000,
                    },
                    'XC': {
                        'XC_FUNCTIONAL': {
                            '_': 'PBE',
                        },
                    },
                    'PRINT': {
                        'E_DENSITY_CUBE': {
                            'STRIDE': '1 1 1',
                            },
                        'MO_CUBES': {
                            'STRIDE': '1 1 1',
                            'WRITE_CUBE': 'F',
                            'NLUMO': 1,
                            'NHOMO': 1,
                        },
                    },
                },
                'SUBSYS': {
                    'KIND': [
                        {'_': 'H', 'BASIS_SET': 'DZVP-MOLOPT-SR-GTH',
                            'POTENTIAL': 'GTH-PBE'},
                        {'_': 'Li', 'BASIS_SET': 'DZVP-MOLOPT-SR-GTH',
                            'POTENTIAL': 'GTH-PBE'},
                        {'_': 'B', 'BASIS_SET': 'DZVP-MOLOPT-SR-GTH',
                            'POTENTIAL': 'GTH-PBE'},
                        {'_': 'N', 'BASIS_SET': 'DZVP-MOLOPT-SR-GTH',
                            'POTENTIAL': 'GTH-PBE'},
                        {'_': 'C', 'BASIS_SET': 'DZVP-MOLOPT-SR-GTH',
                            'POTENTIAL': 'GTH-PBE'},
                        {'_': 'O', 'BASIS_SET': 'DZVP-MOLOPT-SR-GTH',
                            'POTENTIAL': 'GTH-PBE'},
                        {'_': 'F', 'BASIS_SET': 'DZVP-MOLOPT-SR-GTH',
                            'POTENTIAL': 'GTH-PBE'},
                        {'_': 'Si', 'BASIS_SET': 'DZVP-MOLOPT-SR-GTH',
                            'POTENTIAL': 'GTH-PBE'},
                        {'_': 'P', 'BASIS_SET': 'DZVP-MOLOPT-SR-GTH',
                            'POTENTIAL': 'GTH-PBE'},
                        {'_': 'S', 'BASIS_SET': 'DZVP-MOLOPT-SR-GTH',
                            'POTENTIAL': 'GTH-PBE'},
                        {'_': 'Cl', 'BASIS_SET': 'DZVP-MOLOPT-SR-GTH',
                            'POTENTIAL': 'GTH-PBE'},
                        {'_': 'Ni', 'BASIS_SET': 'DZVP-MOLOPT-SR-GTH',
                            'POTENTIAL': 'GTH-PBE'},
                        {'_': 'Cu', 'BASIS_SET': 'DZVP-MOLOPT-SR-GTH',
                            'POTENTIAL': 'GTH-PBE'},
                        {'_': 'Zn', 'BASIS_SET': 'DZVP-MOLOPT-SR-GTH',
                            'POTENTIAL': 'GTH-PBE'},
                        {'_': 'Br', 'BASIS_SET': 'DZVP-MOLOPT-SR-GTH',
                            'POTENTIAL': 'GTH-PBE'},
                        {'_': 'I', 'BASIS_SET': 'DZVP-MOLOPT-SR-GTH',
                            'POTENTIAL': 'GTH-PBE'},
                    ],
                },
            },
        }
        
        scf = params['FORCE_EVAL']['DFT']['SCF']

        if smearing:
            scf['MAX_SCF'] = 500
            scf['ADDED_MOS'] = 1000

            scf['SMEAR'] = {
                '_'                        : 'ON',
                'METHOD'                   : 'FERMI_DIRAC',
                'ELECTRONIC_TEMPERATURE'   : 300,
            }

            scf['DIAGONALIZATION'] = {'ALGORITHM': 'STANDARD'}
            
            scf['MIXING'] = {
                'METHOD'    : 'BROYDEN_MIXING',
                'ALPHA'     : 0.1,
                'NBROYDEN'  : 8,
            }
            
        else:
            scf['OUTER_SCF'] = {
                        'MAX_SCF': 10,
                        'EPS_SCF': 1.0e-6,
            }

            scf['OT'] = {
                    'MINIMIZER': 'CG',
                    'PRECONDITIONER': 'FULL_ALL',
            },

        return params
    
    def return_ddec_parameters (self, atomic_densities="/scratch/snx3000/ongari/atomic_densities/"):
        params = {
            "net charge": 0.0,
            "charge type": "DDEC6",
            "periodicity along A, B, and C vectors" : [True, True, True,],
            "compute BOs" : False,
            "atomic densities directory complete path" : atomic_densities,
            "input filename" : "valence_density",
            "number of core electrons":[
                "1  0",
                "2  0",
                "3  0",
                "4  0",
                "5  2",
                "6  2",
                "7  2",
                "8  2",
                "9  2",
                "10 2",
                "11 2",
                "12 2",
                "13 10",
                "14 10",
                "15 10",
                "16 10",
                "17 10",
                "18 10",
                "19 10",
                "20 10",
                "21 10",
                "22 10",
                "23 10",
                "24 10",
                "25 10",
                "26 10",
                "27 10",
                "28 10",
                "29 18",
                "30 18",
                "35 28",
                "53 46",
                ]
        }
        return params

# EOF
