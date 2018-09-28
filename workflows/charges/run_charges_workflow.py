from aiida.common.example_helpers import test_and_get_code  # noqa
from aiida.orm.data.structure import StructureData  # noqa
from aiida.orm.data.parameter import ParameterData  # noqa
from aiida.orm.data.base import Str
from aiida.work.run import submit

from ase.io import read
from charges import DdecChargesWorkChain

atoms = read('Fe-MOF-74_h111.xyz')
atoms.cell = [[6.96775, 0.00000, 0.00000],
        [-2.33067, 15.22261, 0.00000],
        [ -2.32566, -7.57517, 13.22945]]

structure = StructureData(ase=atoms)
structure.store()

cp2k_options_dict = {
    "resources": {
        "num_machines": 2,
    },
    "max_wallclock_seconds": 8 * 60 * 60,
    }

ddec_options_dict = {
    "resources": {
        "num_machines": 1,
    },
    "max_wallclock_seconds": 8 * 60 * 60,
    "withmpi": False,
    }
cp2k_options = ParameterData(dict=cp2k_options_dict)
ddec_options = ParameterData(dict=ddec_options_dict)
cp2k_code = test_and_get_code('cp2k@fidis', expected_code_type='cp2k')
ddec_code = test_and_get_code('ddec@fidis', expected_code_type='ddec')
submit(DdecChargesWorkChain,
        structure=structure,
        cp2k_code=cp2k_code,
        cp2k_options=cp2k_options,
        ddec_code=ddec_code,
        ddec_options=ddec_options,
        ) 
