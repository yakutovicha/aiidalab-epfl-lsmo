import sys
from ase.io import read,write

from aiida.orm.calculation.work import WorkCalculation
from aiida.orm.data.structure import StructureData
from aiida.orm.utils import CalculationFactory

cp2k_calc = CalculationFactory('cp2k')



def parse_out(cp2k_calc, grep_string, n_word_print):
    path = cp2k_calc.out.retrieved.get_abs_path() + '/path/aiida.out'
    with open(path) as f:
        for line in f.readlines():
            if grep_string in line:
                #print line
                result = line.split()[n_word_print]
    return result


def print_field(to_parse, cp2k_calc=None, title=False):
    for  element in to_parse:
        length = element["print_field"]
        if title:
            print element["title"].rjust(length)[:length],
        else:
            print parse_out(cp2k_calc, element["grep_string"], element["n_word_print"]).rjust(length)[:length],
        print " |",

def get_workcalcs_from_label(label="", type_workcalc=""):
    q = QueryBuilder()
    q.append(StructureData, filters={'label' : {'==': label}}, tag='structure')
    q.append(type_workcalc, descendant_of='structure')
    return list(zip(*q.all())[0])

to_parse = [
        {
            "title"        : "Total energy",
            "grep_string"  : "ENERGY|",
            "n_word_print" : 8, 
            "print_field"  : 14, 
            },
        {
            "title"        : "HOMO - LUMO beta",
            "grep_string"  : "HOMO - LUMO gap [eV]",
            "n_word_print" : 6, 
            "print_field"  : 18, 
            },
        ]


threshold = 1

print "calc/pk |",
print_field(to_parse, title=True)
print " Minimizer |",
print "                                       path                               |",
print ""

list_workcalcn = []

for inp in sys.argv[1:]:
    try:
        pkn = [ load_node(int(inp)) ]
    except:
        pkn = get_workcalcs_from_label(label=inp, type_workcalc=WorkCalculation)

    list_workcalcn += pkn

for workcalc in list_workcalcn:
    inp_struct = workcalc.inp.structure.get_ase().get_chemical_formula(), workcalc.inp.structure.label
    if not workcalc.has_finished():
        print "Skipping underway workchain: {}".format(workcalc.pk)
        continue

    if workcalc.has_failed():
        print "Skipping failed workchain: {}".format(workcalc.pk)
        continue

    list_of_calcs = [ calc for calc in workcalc.get_outputs() if type(calc) == cp2k_calc ]
    
    if len(list_of_calcs) >= threshold:
        print "------------------------------------------------------------------------------------------------------------------"
        print " workcalculation: ", workcalc, inp_struct 
        for number, calculation in enumerate(list_of_calcs):
            print str(number).rjust(7)[:7] + " |",
            print_field(to_parse=to_parse, cp2k_calc=calculation )
            minimizer = calculation.inp.parameters.get_dict()['FORCE_EVAL']['DFT']['SCF']['OT']['MINIMIZER']
            print minimizer.rjust(10)[:10] + " |",
            energy = calculation.out.output_parameters.get_attr('energy')
            path = calculation.out.remote_folder.get_remote_path()
            print path.rjust(73)[:73] + " |",
    
            print ""
