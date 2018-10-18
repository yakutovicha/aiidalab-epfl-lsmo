from __future__ import print_function

from aiida.orm import CalculationFactory, DataFactory
from aiida.orm.code import Code
from aiida.work.run import submit
from aiida.work.workchain import WorkChain, ToContext, Outputs
from aiida_cp2k.workflows import Cp2kDftBaseWorkChain


# calculations 
DdecCalculation = CalculationFactory('ddec')

# data objects
CifData = DataFactory('cif')
ParameterData = DataFactory('parameter')
RemoteData = DataFactory('remote')
StructureData = DataFactory('structure')


default_options = {
    "resources": {
        "num_machines": 1,
        "num_mpiprocs_per_machine": 1,
        },
    "max_wallclock_seconds": 3 * 60 * 60,
    "withmpi": False,
    }

default_ddec_params = {
    "net charge"                               : 0.0,
    "charge type"                              : "DDEC6",
    "periodicity along A, B, and C vectors"    : [True, True, True,],
    "compute BOs"                              : False,
    "atomic densities directory complete path" : "/home/yakutovi/chargemol_09_26_2017/atomic_densities/",
    "input filename"                           : "valence_density",
    "number of core electrons"                 : [
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

class DdecChargesWorkChain(WorkChain):
    """A workchain that computes DDEC charges"""
    @classmethod
    def define(cls, spec):
        super(DdecChargesWorkChain, cls).define(spec)
        #TODO: Change to this when aiida 1.0.0 will be released
        #spec.expose_inputs(Cp2kDftBaseWorkChain, namespace='cp2k', exclude=('structure'))
        spec.input('structure', valid_type=StructureData)
        spec.input('cp2k_code', valid_type=Code)
        spec.input("cp2k_options", valid_type=ParameterData, default=None, required=False)
        spec.input('cp2k_parent_folder', valid_type=RemoteData, default=None, required=False)
        spec.input('ddec_code', valid_type=Code)
        spec.input('ddec_parameters', valid_type=ParameterData, default=ParameterData(dict=default_ddec_params),
                required=False)
        spec.input("ddec_options", valid_type=ParameterData, default=ParameterData(dict=default_options), required=False)

        spec.outline(
                cls.setup,
                cls.run_cp2k_charge_density,
                cls.run_ddec_point_charges,
                cls.return_results,
                )

        spec.output('output_structure', valid_type=CifData, required=False)

    def setup(self):
        """Perform initial setup"""
        pass

    def run_cp2k_charge_density(self):
        """Compute charge-density with CP2K"""
        #TODO Change to this when aiida 1.0.0 will be released
        # inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, namespace='base'))
        # inputs.structure = self.input.structure
        # inputs = prepare_process_inputs(Cp2kDftBaseWorkChain, inputs)
        parameters = ParameterData(dict={
                'FORCE_EVAL':{
                    'DFT':{
                        'PRINT':{
                            'E_DENSITY_CUBE':{
                                'STRIDE': '1 1 1',
                                }
                            },
                        },
                    },
                })
        inputs = {
            'code'       : self.inputs.cp2k_code,
            'structure'  : self.inputs.structure,
            'parameters' : parameters,
            'options'    : self.inputs.cp2k_options,
            '_guess_multiplisity' : True,
            }
        running = submit(Cp2kDftBaseWorkChain, **inputs)
        self.report("pk: {} | Running Cp2kDftBaseWorkChain to compute the charge-density")
        return ToContext(charge_density_calc=Outputs(running))


    def run_ddec_point_charges(self):
        """Compute ddec point charges from precomputed charge-density."""
        charge_density = self.ctx.charge_density_calc['remote_folder']
        #options['prepend_text'] = "export OMP_NUM_THREADS=12"
        inputs = {
            'code'                   : self.inputs.ddec_code,
            'parameters'             : self.inputs.ddec_parameters,
            'charge_density_folder'  : charge_density,
            '_options'               : self.inputs.ddec_options.get_dict(),
            '_label'                 : "run_pointcharges_ddec",
        }

        # Create the calculation process and launch it
        running = submit(DdecCalculation.process(), **inputs)
        self.report("pk: {} | Running ddec to compute point charges based on the charge-density")
        return ToContext(ddec_calc=Outputs(running))

    def return_results(self):
        self.out('output_structure', self.ctx.ddec_calc['structure'])

# EOF
