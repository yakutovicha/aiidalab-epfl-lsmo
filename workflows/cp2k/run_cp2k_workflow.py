from aiida.common.example_helpers import test_and_get_code  # noqa
from aiida.orm.data.structure import StructureData  # noqa
from aiida.orm.data.parameter import ParameterData  # noqa
from aiida.orm.data.base import Str
from aiida.work.run import submit

import ase.build
from cp2k import Cp2kDftBaseWorkChain

atoms = ase.build.molecule('H2O')
atoms.center(vacuum=2.0)
structure = StructureData(ase=atoms)
structure.store()
options_dict = {
    "resources": {
        "num_machines": 1,
        "num_mpiprocs_per_machine": 2,
    },
    "max_wallclock_seconds": 3 * 60 * 60,
    }
options = ParameterData(dict=options_dict)
code = test_and_get_code('cp2k@localhost', expected_code_type='cp2k')
submit(Cp2kDftBaseWorkChain,
        code=code,
        structure=structure,
        options=options,
        ) 
