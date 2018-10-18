import os
import numpy as np
from aiida.common.example_helpers import test_and_get_code  # noqa
from aiida.orm import CalculationFactory, DataFactory
from aiida.orm.data.base import Str
from aiida.work.run import submit

from ase.io import read
from workflows.isotherm import Isotherm

# data objects
ArrayData = DataFactory('array')
ParameterData = DataFactory('parameter')
CifData = DataFactory('cif')

structure = CifData(file=os.getcwd()+'/BONWAD.cif')

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

zr_options_dict = {
    "resources": {
        "num_machines": 1,
        "tot_num_mpiprocs": 1,
    },
    "max_wallclock_seconds": 8 * 60 * 60,
    "withmpi": False,
    }

sigma = 1.6

zeopp_parameters = ParameterData(dict={
    'ha'    : True,
    'res'   : True,
    'sa'    : [sigma, sigma, 100000],
    'volpo' : [sigma, sigma, 100000],
    })

raspa_parameters = ParameterData(dict={
        "GeneralSettings":
        {
        "SimulationType"                   : "MonteCarlo",
        "NumberOfCycles"                   : 2000,
        "NumberOfInitializationCycles"     : 2000,   # 20000

        "PrintEvery"                       : 500,

        "ChargeMethod"                     : "Ewald",
        "CutOff"                           : 14.0,
        "Forcefield"                       : "GenericMOFs",
        "EwaldPrecision"                   : 1e-6,

        "Framework"                        : 0,
        "UnitCells"                        : "2 2 2",
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
        })

cp2k_options = ParameterData(dict=cp2k_options_dict)
ddec_options = ParameterData(dict=ddec_options_dict)
zeopp_options = ParameterData(dict=zr_options_dict)
raspa_options = ParameterData(dict=zr_options_dict)
cp2k_code = test_and_get_code('cp2k@fidis', expected_code_type='cp2k')
ddec_code = test_and_get_code('ddec@fidis', expected_code_type='ddec')
zeopp_code = test_and_get_code('zeopp@deneb', expected_code_type='zeopp.network')
raspa_code = test_and_get_code('raspa@deneb', expected_code_type='raspa')


pressures = ArrayData()
pressures.set_array("pressures", np.array([1e4, 1e5, 1e6, 1e7]))
submit(Isotherm,
        structure=structure,
        probe_molecule=ParameterData(dict={"sigma":1.5}),
        pressures=pressures,
        cp2k_code=cp2k_code,
        cp2k_options=cp2k_options,
        ddec_code=ddec_code,
        ddec_options=ddec_options,
        zeopp_code=zeopp_code,
        zeopp_options=zeopp_options,
        raspa_code=raspa_code,
        raspa_parameters=raspa_parameters,
        raspa_options=raspa_options,
        _usecharges=True,
        )
