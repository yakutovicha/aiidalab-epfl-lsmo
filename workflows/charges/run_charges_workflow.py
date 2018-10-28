from aiida.common.example_helpers import test_and_get_code  # noqa
from aiida.orm.data.structure import StructureData  # noqa
from aiida.orm.data.parameter import ParameterData  # noqa
from aiida.orm.data.base import Str
from aiida.work.run import submit

from ase.io import read
from workflows.charges import DdecCp2kChargesWorkChain

atoms = read('Cu-MOF-74.cif')

structure = StructureData(ase=atoms)
structure.store()

cp2k_options = {
    "resources": {
        "num_machines": 1,
    },
    "max_wallclock_seconds": 1 * 60 * 60,
    }

ddec_options = {
    "resources": {
        "num_machines": 1,
    },
    "max_wallclock_seconds": 1 * 60 * 60,
    "withmpi": False,
    }
cp2k_code = test_and_get_code('cp2k@fidis-debug', expected_code_type='cp2k')
ddec_code = test_and_get_code('ddec@fidis-debug', expected_code_type='ddec')
submit(DdecCp2kChargesWorkChain,
        structure=structure,
        cp2k_code=cp2k_code,
        _cp2k_options=cp2k_options,
#        cp2k_parent_folder=load_node(5337),
        ddec_code=ddec_code,
        _ddec_options=ddec_options,
        ) 
