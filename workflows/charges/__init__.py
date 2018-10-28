from __future__ import print_function

from aiida.orm import CalculationFactory, DataFactory
from aiida.orm.code import Code
from aiida.work import workfunction as wf
from aiida.work.run import submit
from aiida.work.workchain import WorkChain, ToContext, Outputs
from aiida_cp2k.workflows import Cp2kDftBaseWorkChain
from copy import deepcopy

# calculations 
DdecCalculation = CalculationFactory('ddec')

# data objects
CifData = DataFactory('cif')
ParameterData = DataFactory('parameter')
RemoteData = DataFactory('remote')
StructureData = DataFactory('structure')


def dict_merge(dct, merge_dct):
    """ Taken from https://gist.github.com/angstwad/bf22d1822c38a92ec0a9
    Recursive dict merge. Inspired by :meth:``dict.update()``, instead of
    updating only top-level keys, dict_merge recurses down into dicts nested
    to an arbitrary depth, updating keys. The ``merge_dct`` is merged into
    ``dct``.
    :param dct: dict onto which the merge is executed
    :param merge_dct: dct merged into dct
    :return: None
    """
    import collections
    for k, v in merge_dct.iteritems():
        if (k in dct and isinstance(dct[k], dict)
                and isinstance(merge_dct[k], collections.Mapping)):
            dict_merge(dct[k], merge_dct[k])
        else:
            dct[k] = merge_dct[k]

@wf
def merge_ParameterData(p1, p2):
    p1_dict = p1.get_dict()
    p2_dict = p2.get_dict()
    dict_merge(p1_dict, p2_dict)
    return ParameterData(dict=p1_dict).store()

empty_pd = ParameterData(dict={}).store()

default_ddec_options = {
    "resources": {
        "num_machines": 1,
        "num_mpiprocs_per_machine": 1,
        },
    "max_wallclock_seconds": 3 * 60 * 60,
    "withmpi": False,
    }

default_ddec_params = ParameterData(dict={
    "net charge"                               : 0.0,
    "charge type"                              : "DDEC6",
    "periodicity along A, B, and C vectors"    : [True, True, True,],
    "compute BOs"                              : False,
    "atomic densities directory complete path" : "/home/yakutovi/chargemol_09_26_2017/atomic_densities/",
    "input filename"                           : "valence_density",
}).store()

@wf
def extract_core_electrons(charge_density_folder):
    fn = charge_density_folder.get_abs_path() + '/path/aiida.out'
    with open(fn) as f:
        content = f.readlines()
    for n_line, line in enumerate(content):
        if "- Atoms:" in line:
            n_atoms = int(line.split()[2])
        if "Atom  Kind  Element       X           Y           Z" in line:
            break
    res = { e.split()[3] : round(float(e.split()[7])) for e in content[n_line+2:n_line+n_atoms+2]}
    res = [ k+' '+str(int(k)-int(v)) for k,v in res.items() ]
    return ParameterData(dict={'number of core electrons': res}).store()


class DdecCp2kChargesWorkChain(WorkChain):
    """A workchain that computes DDEC charges"""
    @classmethod
    def define(cls, spec):
        super(DdecCp2kChargesWorkChain, cls).define(spec)
        #TODO: Change to this when aiida 1.0.0 will be released
        #spec.expose_inputs(Cp2kDftBaseWorkChain, namespace='cp2k', exclude=('structure'))

        # specify the inputs of the workchain
        spec.input('structure', valid_type=StructureData)
        spec.input('cp2k_code', valid_type=Code)
        spec.input("cp2k_parameters", valid_type=ParameterData, default=empty_pd)
        spec.input("_cp2k_options", valid_type=dict, default=None, required=False)
        spec.input('cp2k_parent_folder', valid_type=RemoteData, default=None, required=False)
        spec.input('ddec_code', valid_type=Code)
        spec.input('ddec_parameters', valid_type=ParameterData, default=default_ddec_params, required=False)
        spec.input("_ddec_options", valid_type=dict, default=deepcopy(default_ddec_options), required=False)

        # specify the chain of calculations to be performed
        spec.outline(
                cls.setup,
                cls.run_cp2k,
                cls.prepare_ddec,
                cls.run_ddec,
                cls.return_results,
                )

        # specify the outputs of the workchain
        spec.output('output_structure', valid_type=CifData, required=False)

    def setup(self):
        """Perform initial setup"""
        pass

    def run_cp2k(self):
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
                                '_'     : 'ON',
                                'STRIDE': '1 1 1',
                                }
                            },
                        },
                    },
                }).store()
        parameters = merge_ParameterData(parameters, self.inputs.cp2k_parameters)
        inputs = {
            'code'                : self.inputs.cp2k_code,
            'structure'           : self.inputs.structure,
            'parameters'          : parameters,
            '_label'              : 'Cp2kDftBaseWorkChain',
            }
        try:
            inputs['_options'] = self.inputs._cp2k_options
        except:
            pass
        try:
            inputs['parent_folder'] = self.inputs.cp2k_parent_folder
        except:
            self.report("Folder with converged wavefunction was not provided, CP2K will compute it from scratch")

        running = submit(Cp2kDftBaseWorkChain, **inputs)
        self.report("pk: {} | Running Cp2kDftBaseWorkChain to compute the charge-density".format(running.pid))
        return ToContext(charge_density_calc=Outputs(running))


    def prepare_ddec(self):
        """Prepare inputs for ddec point charges calculation."""
        # extract number of core electrons from the cp2k output
        self.report("self.ctx.charge_density_calc: {}".format(str(self.ctx.charge_density_calc)))
        core_e = extract_core_electrons(self.ctx.charge_density_calc['retrieved'])
        # prepare input dictionary
        self.ctx.ddec_inputs = {
            'code'                   : self.inputs.ddec_code,
            'charge_density_folder'  : self.ctx.charge_density_calc['remote_folder'],
            'parameters'             : merge_ParameterData(self.inputs.ddec_parameters, core_e),
            '_options'               : self.inputs._ddec_options,
            '_label'                 : "DdecCalculation",
        }

    def run_ddec(self):
        """Compute ddec point charges from precomputed charge-density."""
        # Create the calculation process and launch it
        running = submit(DdecCalculation.process(), **self.ctx.ddec_inputs)
        self.report("pk: {} | Running ddec to compute point charges based on the charge-density".format(running.pid))
        return ToContext(ddec_calc=Outputs(running))

    def return_results(self):
        self.out('output_structure', self.ctx.ddec_calc['structure'])

# EOF
