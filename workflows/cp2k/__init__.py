from aiida.common.extendeddicts import AttributeDict
from aiida.work.run import submit
from aiida.work.workchain import WorkChain, Outputs
from aiida.orm import Code
from aiida.orm.data.base import Str
from aiida.orm.data.parameter import ParameterData
from aiida.orm.data.remote import RemoteData
from aiida.orm.data.structure import StructureData
from aiida.orm.utils import CalculationFactory
from aiida.work.workchain import ToContext, if_, while_

Cp2kCalculation = CalculationFactory('cp2k')

spin = {
        "H"  : 0.0,
        "Li" : 0.0,
        "Be" : 0.0,
        "B"  : 0.0,
        "C"  : 0.0,
        "N"  : 0.0,
        "O"  : 0.0,
        "F"  : 0.0,
        "Na" : 0.0,
        "Mg" : 0.0,
        "Al" : 0.0,
        "Si" : 0.0,
        "P"  : 0.0,
        "S"  : 0.0,
        "Cl" : 0.0,
        "K"  : 0.0,
        "Ca" : 0.0,
        "Sc" : 1.0/2.0,
        "Ti" : 2.0/2.0,
        "V"  : 3.0/2.0,
        "Cr" : 4.0/2.0,
        "Mn" : 5.0/2.0,
        "Fe" : 4.0/2.0,
        "Co" : 3.0/2.0,
        "Ni" : 2.0/2.0,
        "Cu" : 1.0/2.0,
        "Zn" : 0.0,
        "Zr" : 2.0/2.0,
        }

params = {
    'FORCE_EVAL': {
        'METHOD': 'Quickstep',
        'DFT': {
            'CHARGE': 0,
            'BASIS_SET_FILE_NAME': 'BASIS_MOLOPT',
            'POTENTIAL_FILE_NAME': 'GTH_POTENTIALS',
            'RESTART_FILE_NAME'  : './parent_calc/aiida-RESTART.wfn',
            'QS': {
                'METHOD':'GPW',
            },
            'POISSON': {
                'PERIODIC': 'XYZ',
            },
            'MGRID': {
                'CUTOFF':     600,
                'NGRIDS':       4,
                'REL_CUTOFF':  50,
            },
            'SCF':{
                'SCF_GUESS': 'ATOMIC',
                'EPS_SCF': 1.0e-6,
                'MAX_SCF': 50,
                'MAX_ITER_LUMO': 10000,
                'OT':{
                    'MINIMIZER': 'DIIS',
                    'PRECONDITIONER': 'FULL_ALL',
                    },
                'OUTER_SCF':{
                    'EPS_SCF': 1.0e-6,
                    'MAX_SCF': 10,
                    },
            },
            'XC': {
                'XC_FUNCTIONAL': {
                    '_': 'PBE',
                },
            },
            'PRINT': {
                'MO_CUBES': {   # this is to print the band gap
                    'STRIDE': '1 1 1',
                    'WRITE_CUBE': 'F',
                    'NLUMO': 1,
                    'NHOMO': 1,
                },
            },
        },
        'SUBSYS': {
        },
    },
}


def last_scf_loop(fpath):
    """
    Simple function that extracts all the output starting from the last SCF
    loop.
    """
    with open(fpath) as f:
        content = f.readlines()
    # find the last scf loop in the cp2k output file
    for n, line in enumerate(reversed(content)):
        if "SCF WAVEFUNCTION OPTIMIZATION" in line:
            break
    return content[-n-1:]

def scf_converged(fpath):
    content = last_scf_loop(fpath)
    for line in content:
        if "SCF run converged in" in line:
            return True
    return False

def scf_was_diverging(fpath):
    content = last_scf_loop(fpath)
    for line in content:
        if "Minimizer" in line and "CG" in line:
            grep_string = "OT CG"
            break

        elif "Minimizer" in line and "DIIS" in line:
            grep_string = "OT DIIS"
            break
    
    n_change = 7
    difference = []
    n_positive = 0
    for line in content:
        if grep_string in line:
            difference.append(line.split()[n_change])
    for number in difference[-12:]:
        if float(number) > 0:
            n_positive +=1

    if n_positive>5:
        return True
    return False


def get_multiplicity(structure):
    multiplicity = 1
    all_atoms = structure.get_ase().get_chemical_symbols()
    for key, value in spin.iteritems():
        multiplicity += all_atoms.count(key) * value * 2.0
    return int(round(multiplicity))

def get_atom_kinds(structure):
    kinds = []
    all_atoms = set(structure.get_ase().get_chemical_symbols())
    for a in all_atoms:
        kinds.append({
            '_': a,
            'BASIS_SET': 'DZVP-MOLOPT-SR-GTH',
            'POTENTIAL': 'GTH-PBE',
            'MAGNETIZATION': spin[a]*2.0,
            })
    return kinds

options = {
    "resources": {
        "num_machines": 1,
        "num_mpiprocs_per_machine": 12,
    },
    "max_wallclock_seconds": 3 * 60 * 60,
    }
class Cp2kDftBaseWorkChain(WorkChain):
    """
    A base workchain to be used for DFT calculations with CP2K
    """
    @classmethod
    def define(cls, spec):
        super(Cp2kDftBaseWorkChain, cls).define(spec)
        spec.input('code', valid_type=Code)
        spec.input('structure', valid_type=StructureData)
        spec.input("parameters", valid_type=ParameterData,
                default=ParameterData(dict=params))
        spec.input("options", valid_type=ParameterData,
                default=ParameterData(dict=options))
        spec.input('parent_folder', valid_type=RemoteData, required=False) 
        
        spec.outline(
            cls.setup,
            while_(cls.should_run_calculation)(
                cls.prepare_calculation,
                cls.run_calculation,
                cls.inspect_calculation,
            ),
            cls.return_results,
        )
        spec.output('output_structure', valid_type=StructureData, required=False)
#        spec.output('output_parameters', valid_type=ParameterData)
#        spec.output('remote_folder', valid_type=RemoteData)

    def setup(self):
        self.ctx.done = False
        self.ctx.restart_calc = None
        # TODO: use restart from self.inputs.parent_folder
        self.ctx.nruns = 0

#       This (below) doesn't work :(
#       How one can assign labels to the workchain?
#        if self.inputs.label:
#            self.label = self.inputs.label

        self.ctx.structure = self.inputs.structure
        self.ctx.parameters = self.inputs.parameters.get_dict()
        self.ctx.options = self.inputs.options.get_dict()

        multiplicity = get_multiplicity(self.inputs.structure)
        kinds = get_atom_kinds(self.inputs.structure)

        self.ctx.parameters['FORCE_EVAL']['DFT']['MULTIPLICITY'] = multiplicity
        self.ctx.parameters['FORCE_EVAL']['SUBSYS']['KIND'] = kinds

        # TODO: implement merging of the input cp2k-parameters with the default
        # ones


    def should_run_calculation(self):
        return not self.ctx.done

    def prepare_calculation(self):
        """Prepare all the neccessary input links to run the calculation"""
        p = ParameterData(dict=self.ctx.parameters)
        p.store()

        self.ctx.inputs = AttributeDict({
            'code': self.inputs.code,
            'structure'  : self.ctx.structure,
            'parameters' : p,
            '_options'    : self.ctx.options,
            '_label'      : "DFTwithCP2k",
            })

        # restart from the previous calculation only if the necessary data are
        # provided
        if self.ctx.restart_calc:
            self.ctx.inputs['parent_folder'] = self.ctx.calculation['remote_folder']

    def run_calculation(self): 
        """Run cp2k calculation."""
        
#        with open("/home/epfl/work/mc-epfl-lsmo/workflows/cp2k/workflow.output", "w") as f:
#            f.write(str(self.ctx.inputs))
        
        # Create the calculation process and launch it
        process = Cp2kCalculation.process()
        future  = submit(process, **self.ctx.inputs)
        self.report("pk: {} | Running DFT calculation with"
                " cp2k".format(future.pid))
        self.ctx.nruns += 1
        return ToContext(calculation=Outputs(future))
    
    def inspect_calculation(self):
        """
        Analyse the reults of CP2K calculation and decide weahter there is a
        need to restart it. If yes, then decide exactly how to restart the
        calculation.
        """
        # I will try to disprove those statements. I will not succeed in doing
        # so - the calculation will be considered as completed
        converged_geometry = True
        converged_scf = True
        exceeded_time = False

        # File to analyze
        outfile = self.ctx.calculation['retrieved'].get_abs_path() + '/path/aiida.out'
        self.ctx.restart_calc = self.ctx.calculation['remote_folder']

        # First (and the simplest) check is whether the runtime was exceeded
        exceeded_time = self.ctx.calculation['output_parameters'].dict['exceeded_walltime']
        if exceeded_time:
            self.report("The time of the cp2k calculation has been exceeded")
        else:
            self.report("The time of the cp2k calculation has NOT been exceeded")

        # Second check is whether the last SCF did converge
        converged_scf = scf_converged(outfile)
        self.ctx.parameters['FORCE_EVAL']['DFT']['SCF']['SCF_GUESS'] = 'RESTART'
        if not converged_scf and scf_was_diverging(outfile):
            # If, however, scf was even diverging I should go for more robust
            # minimizer.
            # Aslo, to avoid being trapped in the wrong minimum I restart
            # from atomic wavefunctions.
            self.report("Going for more robust (but slow) SCF minimizer")
            self.ctx.parameters['FORCE_EVAL']['DFT']['SCF']['SCF_GUESS'] = 'ATOMIC'
            self.ctx.parameters['FORCE_EVAL']['DFT']['SCF']['MAX_SCF'] = 2000
            self.ctx.parameters['FORCE_EVAL']['DFT']['SCF']['OT']['MINIMIZER'] = 'CG' 
            self.ctx.parameters['FORCE_EVAL']['DFT']['SCF']['OUTER_SCF']['MAX_SCF'] = 0

            # TODO: I may also look for the forces here. For example a very
            # strong force may cause convergence problems, needs to be
            # implemented
                

       # Third check:
       # TODO: check for the geometry convergence/divergence problems
       # useful for geo/cell-opt restart
       # if aiida-1.restart in retrieved (folder):
       #    self.ctx.parameters['EXT_RESTART'] = {'RESTART_FILE_NAME': './parent_calc/aiida-1.restart'}

        if converged_geometry and converged_scf and not exceeded_time:
            self.report("Calculation converged, terminating the workflow")
            self.ctx.done = True

        

    def return_results(self):
        self.out('output_structure', self.ctx.structure)
