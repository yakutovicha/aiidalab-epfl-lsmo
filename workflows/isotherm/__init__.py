from aiida.orm import CalculationFactory, DataFactory


from aiida.work.workchain import WorkChain, ToContext, while_, Outputs
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
        spec.input("parameters", valid_type=ParameterData, required=True)
        spec.input("pressures", valid_type=ArrayData, required=True)
        spec.input("structure", valid_type=CifData, required=True)
        spec.input("zeopp_codename", valid_type=Str, required=True)
        spec.input("raspa_codename", valid_type=Str, required=True)
        spec.input("_interactive", valid_type=bool, required=False, default=False)
        
        # The outline describes the business logic that defines
        # which steps are executed in what order and based on
        # what conditions. Each `cls.method` is implemented below
        spec.outline(
            cls.init,
            cls.run_zeopp,
            cls.parse_zeopp,
            while_(cls.should_run_raspa)(
                cls.run_raspa,
                cls.parse_raspa,
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
        self.ctx.options = {
            "resources": {
                "num_machines": 1,
                "tot_num_mpiprocs": 1,
                "num_mpiprocs_per_machine": 1,
            },
            "max_wallclock_seconds": 30 * 60,
#            "queue_name":"serial",
        }
        self.ctx.parameters = self.inputs.parameters.get_dict()
        if self.inputs._interactive == True:
            self.fig, self.ax = plt.subplots(1,1)
            self.ax.set_xlabel(u"Pressure [bar]")
            self.ax.set_ylabel(u"Loading average [molecules/unit cell]")


    def should_run_raspa(self):
        """
        This is the main condition of the while loop, as defined
        in the outline of the Workchain. We only run another
        raspa calculation if the current iteration is smaller than
        the total number of pressures we want to compute
        """
        return self.ctx.p < len(self.ctx.pressures)

    def run_zeopp(self):
        """
        This is the main function that will perform a raspa
        calculation for the current pressure
        """

        NetworkParameters = DataFactory('zeopp.parameters')
        # Create the input dictionary
        inputs = {
            'code'       : Code.get_from_string(self.inputs.zeopp_codename.value),
            'structure'  : self.inputs.structure,
            'parameters' : NetworkParameters(dict={'volpo': [1.4, 1.4, 50000],}),
            '_options'   : self.ctx.options,
        }

        # Create the calculation process and launch it
        self.report("Running zeo++ calculation")
        process = ZeoppCalculation.process()
        future  = submit(process, **inputs)

        return ToContext(zeopp=Outputs(future))

    def parse_zeopp(self):
        """
        Extract the pressure and loading average of the last completed raspa calculation
        """
        self.ctx.parameters['GeneralSettings']['HeliumVoidFraction'] = self.ctx.zeopp["pore_volume_volpo"].dict.POAV_Volume_fraction

    def run_raspa(self):
        """
        This is the main function that will perform a raspa
        calculation for the current pressure
        """
        pressure = self.ctx.pressures[self.ctx.p]
        self.ctx.parameters['GeneralSettings']['ExternalPressure'] = pressure
        if self.ctx.prev_pk is not None:
            self.ctx.parameters['GeneralSettings']['RestartFile'] = True
            self.ctx.parameters['restart_pk'] = self.ctx.prev_pk

        # Create the input dictionary
        inputs = {
            'code'       : Code.get_from_string(self.inputs.raspa_codename.value),
            'structure'  : self.inputs.structure,
            'parameters' : ParameterData(dict=self.ctx.parameters),
            '_options'   : self.ctx.options,
        }

        # Create the calculation process and launch it
        self.report("Running raspa for the pressure {} [bar]".format(pressure/1e5))
        process = RaspaCalculation.process()
        future  = submit(process, **inputs)

        self.ctx.p += 1
        self.ctx.prev_pk = future.pid

        return ToContext(raspa=Outputs(future))

    def parse_raspa(self):
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
        """
        self.He=ipw.FloatText(
            value=0.347,
            step=0.01,
            description='Helium Void Fraction:',
            disabled=False,
            layout=ipw.Layout(width="280px"),
            style = self.style
        )
        
        self.cellnx = ipw.IntText(
            value=1,
            description='Number of unit cells: X',
            disabled=False,
            layout=ipw.Layout(width="240px"),
            style = self.style
        )

        self.cellny = ipw.IntText(
            value=1,
            description='Y',
            disabled=False,
            layout=ipw.Layout(width="100px"),
            style = {"description_width":"40px"}
        )

        self.cellnz = ipw.IntText(
            value=1,
            description='Z',
            disabled=False,
            layout=ipw.Layout(width="100px"),
            style = {"description_width":"40px"}
        )
        """
        self.cutoff = ipw.FloatText(
            value=12.0,
            step=0.2,
            description='Cutoff [A]:',
            disabled=False,
            layout=self.layout,
            style=self.style,
        )
        
        self.ff = ipw.Dropdown(
            options=('GenericZeolites', 'tcc', "TraPPE", "GenericMOFs", "GarciaPerez2016" ),
            value='tcc',
            description='Select forcefield:',
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
            #self.He,
            self.cutoff,
        ])


    def return_parameters (self, cif):
        from math import cos, sin, sqrt, pi
        import numpy as np
        deg2rad=pi/180.

        
        a = float(cif.values['image0']['_cell_length_a'])
        b = float(cif.values['image0']['_cell_length_b'])
        c = float(cif.values['image0']['_cell_length_c'])
        
        alpha = float(cif.values['image0']['_cell_angle_alpha'])*deg2rad
        beta  = float(cif.values['image0']['_cell_angle_beta'])*deg2rad
        gamma = float(cif.values['image0']['_cell_angle_gamma'])*deg2rad
        
        
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
                "MoleculeName"                     : "methane",
                "MoleculeDefinition"               : "TraPPE",
                "TranslationProbability"           : 0.5,
                "ReinsertionProbability"           : 0.5,
                "SwapProbability"                  : 1.0,
                "CreateNumberOfMolecules"          : 0,
                }],
                }
        return params
