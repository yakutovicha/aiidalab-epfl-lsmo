from aiida.orm import CalculationFactory, DataFactory
from aiida.orm import load_node
from aiida.orm.code import Code
from aiida.orm.data.base import Bool, Str
from aiida.work import workfunction as wf
from aiida.work.run import run, submit
from aiida.work.workchain import WorkChain, ToContext, if_, while_, Outputs

# subworkflows
from workflows.charges import DdecChargesWorkChain
from aiida_raspa.workflows import RaspaConvergeWorkChain

# calculations
ZeoppCalculation = CalculationFactory('zeopp.network')

# data objects
ArrayData = DataFactory('array')
CifData = DataFactory('cif')
NetworkParameters = DataFactory('zeopp.parameters')
ParameterData = DataFactory('parameter')
RemoteData = DataFactory('remote')
StructureData = DataFactory('structure')

@wf
def from_cif_to_structuredata(cif):
    """Helper function that converts CifData object into StructureData"""
    a = cif.get_ase()
    return StructureData(ase=a)


class Isotherm(WorkChain):
    """Workchain that for a given matherial will compute an isotherm of a
    certain gaz adsorption."""

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
        spec.input("cp2k_options", valid_type=ParameterData, default=None, required=False)
        spec.input('cp2k_parent_folder', valid_type=RemoteData, default=None, required=False)

        # ddec
        spec.input('ddec_code', valid_type=Code)
        spec.input("ddec_options", valid_type=ParameterData, default=None, required=False)

        # zeopp
        spec.input('zeopp_code', valid_type=Code)
        spec.input("zeopp_options", valid_type=ParameterData, default=None, required=False)

        # raspa
        spec.input("raspa_code", valid_type=Code)
        spec.input("raspa_parameters", valid_type=ParameterData)
        spec.input("raspa_options", valid_type=ParameterData, default=None, required=False)

        # settings
        spec.input("_interactive", valid_type=bool, default=False, required=False)
        spec.input("_usecharges", valid_type=bool, default=False, required=False)

        # workflow
        spec.outline(
            cls.init,
            if_(cls.should_use_charges)(
                cls.run_point_charges,
                cls.parse_point_charges,
            ),
            cls.run_geom_zeopp,
            cls.parse_geom_zeopp,
            cls.run_henry_raspa,
            while_(cls.should_run_loading_raspa)(
                cls.run_loading_raspa,
                cls.parse_loading_raspa,
            ),
            cls.return_results,
        )

        # TODO: once the workflow is ready, explicitely specify the outputs
        spec.dynamic_output()

    def init(self):
        """Initialize variables and the pressures we want to compute"""
        self.ctx.structure = self.inputs.structure
        self.ctx.pressures = self.inputs.pressures.get_array("pressures")
        self.ctx.current_p_index = 0
        self.ctx.result = []

        self.ctx.raspa_parameters = self.inputs.raspa_parameters.get_dict()

        if self.inputs._usecharges:
            self.ctx.raspa_parameters['GeneralSettings']['UseChargesFromCIFFile'] = "yes"

        self.ctx.restart_raspa_calc = None

    def should_use_charges(self):
        """Whether it is needed to employ charges"""
        return self.inputs._usecharges

    def run_point_charges(self):
        """Compute the charge-density of a structure that can be later
        used for extracting ddec point charges."""

        inputs = {
            'structure'       : from_cif_to_structuredata(self.ctx.structure),
            'cp2k_code'       : self.inputs.cp2k_code,
            'cp2k_options'    : self.inputs.cp2k_options,
            'ddec_code'       : self.inputs.ddec_code,
            '_label'          : "run_point_charges",
        }

        # Create the calculation process and launch it
        running = submit(DdecChargesWorkChain, **inputs)
        self.report("pk: {} | Running cDdecChargesWorkChain to compute the point charges")
        return ToContext(point_charges_calc=Outputs(running))

    def parse_point_charges(self):
        """Extract structure with charges and put it into self.ctx.structure"""
        self.ctx.structure = self.ctx.point_charges_calc['output_structure']

    def run_geom_zeopp(self):
        """This is the main function that will perform a raspa calculation for the current pressure"""
        # network parameters
        sigma = self.inputs.probe_molecule.dict.sigma
        NetworkParameters = DataFactory('zeopp.parameters')
        params = NetworkParameters({
            'ha': True,
            'res': True,
            'sa': [sigma, sigma, 100000],
            'volpo': [sigma, sigma, 100000],
        })

        # Create the input dictionary
        inputs = {
            'code'       : self.inputs.zeopp_code,
            'structure'  : self.ctx.structure,
            'parameters' : params,
            '_options'   : self.inputs.zeopp_options.get_dict(),
            '_label'     : "run_geom_zeopp",
        }

        # Create the calculation process and launch it
        running = submit(ZeoppCalculation.process(), **inputs)
        self.report("pk: {} | Running geometry analysis with zeo++".format(running.pid))

        return ToContext(zeopp=Outputs(running))

    def parse_geom_zeopp(self):
        """Extract the pressure and loading average of the last completed raspa calculation"""
        self.ctx.raspa_parameters['GeneralSettings']['HeliumVoidFraction'] = \
        self.ctx.zeopp["pore_volume_volpo"].dict.POAV_Volume_fraction

    def should_run_loading_raspa(self):
        """We run another raspa calculation only if the current iteration is smaller than
        the total number of pressures we want to compute"""
        return self.ctx.current_p_index < len(self.ctx.pressures)

    def run_loading_raspa(self):
        """This function will run RaspaConvergeWorkChain for the current pressure"""
        pressure = self.ctx.pressures[self.ctx.current_p_index]
        self.ctx.raspa_parameters['GeneralSettings']['ExternalPressure'] = pressure

        # Create the input dictionary
        inputs = {
            'code'       : self.inputs.raspa_code,
            'structure'  : self.ctx.structure,
            'parameters' : ParameterData(dict=self.ctx.raspa_parameters),
            'options'   : self.inputs.raspa_options,
            '_label'     : "run_loading_raspa",
        }

        if self.ctx.restart_raspa_calc is not None:
            inputs['retrieved_parent_folder'] = self.ctx.restart_raspa_calc

        # Create the calculation process and launch it
        running = submit(RaspaConvergeWorkChain, **inputs)
        self.report("pk: {} | Running raspa for the pressure {} [bar]".format(running.pid, pressure/1e5))
        self.ctx.current_p_index += 1
        return ToContext(raspa_loading=Outputs(running))

    def parse_loading_raspa(self):
        """Extract the pressure and loading average of the last completed raspa calculation"""
        self.ctx.restart_raspa_calc = self.ctx.raspa_loading['retrieved_parent_folder']
        pressure = self.ctx.raspa_parameters['GeneralSettings']['ExternalPressure']/1e5
        loading_average = self.ctx.raspa_loading["component_0"].dict.loading_absolute_average
        self.ctx.result.append((pressure, loading_average))

    def run_henry_raspa(self):
        """This is the main function that will perform a raspaa calculation for the current pressure"""
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
            'code'       : self.inputs.raspa_code,
            'structure'  : self.ctx.structure,
            'parameters' : ParameterData(dict=raspa_parameters),
            'options'    : self.inputs.raspa_options,
            '_label'     : "run_henry_raspa",
        }

        # Create the calculation process and launch it
        running = submit(RaspaConvergeWorkChain, **inputs)
        self.report("pk: {} | Running raspa for the Henry coefficients".format(running.pid))

        return ToContext(raspa_henry=Outputs(running))

    def return_results(self):
        """Attach the results of the raspa calculation and the initial structure to the outputs."""
        self.out("result", ParameterData(dict={"isotherm": self.ctx.result}))
        self.report("Workchain <{}> completed successfully".format(self.calc.pk))
        return

# EOF
