from aiida.orm import CalculationFactory, DataFactory


from aiida.work.workchain import WorkChain, ToContext, if_, while_, Outputs
from aiida.work.run import run, submit


from aiida.orm.data.parameter import ParameterData
from aiida.orm.data.array import ArrayData
from aiida.orm.data.cif import CifData
from aiida.orm.data.base import Bool, Str

from aiida.orm.code import Code
#from aiida.orm
RaspaCalculation = CalculationFactory('raspa')
ZeoppCalculation = CalculationFactory('zeopp.network')


import matplotlib.pyplot as plt


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
        
        # First we define the inputs, specifying the type we expect
        options = {
            "resources": {
                "num_machines": 1,
                "tot_num_mpiprocs": 1,
                "num_mpiprocs_per_machine": 1,
            },
            "max_wallclock_seconds": 30 * 60,
            "max_memory_kb": 2e6,
#            "queue_name":"serial",
        }
        spec.input("probe_molecule", valid_type=ParameterData, required=True)
        spec.input("parameters", valid_type=ParameterData, required=True)
        spec.input("pressures", valid_type=ArrayData, required=True)
        spec.input("structure", valid_type=CifData, required=True)
        spec.input("zeopp_codename", valid_type=Str, required=True)
        spec.input("raspa_codename", valid_type=Str, required=True)
        spec.input("_options", valid_type=dict, required=False, default=options)
        spec.input("_interactive", valid_type=bool, required=False, default=False)
        
        # The outline describes the business logic that defines
        # which steps are executed in what order and based on
        # what conditions. Each `cls.method` is implemented below
        spec.outline(
            cls.init,
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
        self.ctx.p = 0
        self.ctx.prev_pk = None
        self.ctx.pressures = self.inputs.pressures.get_array("pressures")
        self.ctx.result = []

        self.ctx.parameters = self.inputs.parameters.get_dict()
        if self.inputs._interactive == True:
            self.fig, self.ax = plt.subplots(1,1)
            self.ax.set_xlabel(u"Pressure [bar]")
            self.ax.set_ylabel(u"Loading average [molecules/unit cell]")
       
    def run_geom_zeopp(self):
        """
        This is the main function that will perform a raspa
        calculation for the current pressure
        """

        NetworkParameters = DataFactory('zeopp.parameters')
        # Create the input dictionary
        sigma = self.inputs.probe_molecule.get_dict()['sigma']
        params = {
            'ha': True,
            'res': True,
            'sa': [sigma, sigma, 100000],
            'volpo': [sigma, sigma, 100000],
        }
        inputs = {
            'code'       : Code.get_from_string(self.inputs.zeopp_codename.value),
            'structure'  : self.inputs.structure,
            'parameters' : NetworkParameters(dict=params),
            '_options'   : self.inputs._options,
            '_label'     : "run_geom_zeopp",
        }

        # Create the calculation process and launch it
        process = ZeoppCalculation.process()
        future  = submit(process, **inputs)
        self.report("pk: {} | Running geometry analysis with zeo++".format(future.pid))

        return ToContext(zeopp=Outputs(future))

    def parse_geom_zeopp(self):
        """
        Extract the pressure and loading average of the last completed raspa calculation
        """
        self.ctx.parameters['GeneralSettings']['HeliumVoidFraction'] = self.ctx.zeopp["pore_volume_volpo"].dict.POAV_Volume_fraction

    def should_run_block_zeopp(self):
        """
        If the pore non-accessible volume is 0 - there is no need to run
        """
        return self.ctx.zeopp["pore_volume_volpo"].dict.PONAV_Volume_fraction != 0
        
    def run_block_zeopp(self):
        """
        This is the main function that will perform a raspa
        calculation for the current pressure
        """

        NetworkParameters = DataFactory('zeopp.parameters')
        # Create the input dictionary
        sigma = self.inputs.probe_molecule.get_dict()['sigma']
        params = {
            'ha':True,
            'block': [sigma, 200],
        }
        inputs = {
            'code'       : Code.get_from_string(self.inputs.zeopp_codename.value),
            'structure'  : self.inputs.structure,
            'parameters' : NetworkParameters(dict=params),
            '_options'   : self.inputs._options,
            '_label'     : "run_block_zeopp",

        }

        # Create the calculation process and launch it
        process = ZeoppCalculation.process()
        future  = submit(process, **inputs)
        self.report("pk: {} | Running zeo++ block volume calculation".format(future.pid))
        
        self.ctx.block_pk = future.pid

        return ToContext(zeopp=Outputs(future))
    
    def parse_block_zeopp(self):
        """
        Extract the pressure and loading average of the last completed raspa calculation
        """
        self.ctx.parameters['Component'][0]['BlockPockets'] = True
        self.ctx.parameters['Component'][0]['BlockPocketsPk'] = self.ctx.block_pk

    def should_run_loading_raspa(self):
        """
        This is the main condition of the while loop, as defined
        in the outline of the Workchain. We only run another
        raspa calculation if the current iteration is smaller than
        the total number of pressures we want to compute
        """
        return self.ctx.p < len(self.ctx.pressures)

    def run_loading_raspa(self):
        """
        This is the main function that will perform a raspa
        calculation for the current pressure
        """
        pressure = self.ctx.pressures[self.ctx.p]
        self.ctx.parameters['GeneralSettings']['ExternalPressure'] = pressure
        if self.ctx.prev_pk is not None:
            self.ctx.parameters['GeneralSettings']['RestartFile'] = True
            self.ctx.parameters['GeneralSettings']['RestartFilePk'] = self.ctx.prev_pk

        # Create the input dictionary
        inputs = {
            'code'       : Code.get_from_string(self.inputs.raspa_codename.value),
            'structure'  : self.inputs.structure,
            'parameters' : ParameterData(dict=self.ctx.parameters),
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
        pressure = self.ctx.parameters['GeneralSettings']['ExternalPressure']/1e5
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
        parameters = self.inputs.parameters.get_dict()
        for i, comp in enumerate(parameters['Component']):
            name = comp['MoleculeName']
            parameters['Component'][0] = {
                "MoleculeName"                     : name,
                "MoleculeDefinition"               : "TraPPE",
                "WidomProbability"                 : 1.0,
                "CreateNumberOfMolecules"          : 0,
            }
        # Create the input dictionary
        inputs = {
            'code'       : Code.get_from_string(self.inputs.raspa_codename.value),
            'structure'  : self.inputs.structure,
            'parameters' : ParameterData(dict=parameters),
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
            value=12.0,
            step=0.2,
            description='Cutoff [A]:',
            disabled=False,
            layout=self.layout,
            style=self.style,
        )
        
        self.ff = ipw.Dropdown(
            options=('LSMO_DREIDING-TraPPE', 'LSMO_UFF-TraPPE'),
            value='LSMO_DREIDING-TraPPE',
            description='Select forcefield:',
            layout=self.layout,
            style=self.style,
            )
        mol_opt = [('methane', {'molecule':'methane','sigma': 1.86}), ("Carbon dioxide", {'molecule':'CO2','sigma': 1.525}), ("Nitrogen", {'molecule':'N2', 'sigma': 1.86})]
        self.probe_molecule = ipw.Dropdown(
            options=mol_opt,
            #value=0,
            description='Select molecule:',
            layout=self.layout,
            style=self.style,
            )

        self.init_cycles = ipw.IntText(
                value=2000,
                step=1000,
                description = "Number of initialization cycles",
                disabled=False,
                layout=self.layout,
                style=self.style,
            )
        self.prod_cycles = ipw.IntText(
                value=2000,
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


    def return_parameters (self, cif):
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
                "CutOff"                           : 12.0,
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

# EOF
