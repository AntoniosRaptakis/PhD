import numpy as np
from ase.io import read,write
from ase.visualize import view
import os,sys

def cell_formation(structure,factor):
    len_a3 = structure.get_cell_lengths_and_angles()[2]
    len_a1 = structure.get_cell_lengths_and_angles()[0]*(1+factor*0.001)
    len_a2 = structure.get_cell_lengths_and_angles()[1]*(1+factor*0.001)
    angle_a1_a2 = structure.get_cell_lengths_and_angles()[3]
    angle_a1_a3 = structure.get_cell_lengths_and_angles()[4]
    angle_a2_a3 = structure.get_cell_lengths_and_angles()[5]

    structure.set_cell([len_a1,len_a2,len_a3,angle_a1_a2,angle_a1_a3,angle_a2_a3],scale_atoms=True)
    len_a3 = 30
    structure.set_cell([len_a1,len_a2,len_a3,angle_a1_a2,angle_a1_a3,angle_a2_a3])
    return structure

elements = ['Br','C','B','Cl','F','H','I','K','Mg','N','Na','O','P','S','Zn']

elements_orbitals = {'H':'H="s"',
                     'C':'C="p"',
              	     'N':'N="p"',
                     'O':'O="p"',
                     'B':'B="p"',
                     'S':'S="d"'}

initial_structure = read('anthracene.gen')
write('anthracene.gen',initial_structure)

#searches for the chemical elements of the structure
symbols=initial_structure.get_chemical_symbols()
structure_elements=[symbols[0]]
for symbol in symbols:
    if symbol not in structure_elements:
        structure_elements.append(symbol)

#clarifies if there is paramentrization for the specific elements included
if all(x in elements for x in structure_elements):
    exist='True'
else:
    exist='False'
    sys.exit('It does not exist an element!')

# define the name of the directory to be created
path = "/scratch/ws/anra181a-SPECint/Bulk_modulus/porphyrin/anthracene/"

for i in range(-8,10,2):

    s = read('anthracene.gen')
    directory_name = "optimization_compress_stretch_"+str(i*0.001)
    os.mkdir(os.path.join(path,directory_name))
    os.chdir(directory_name)
    fin = open("../dftb1_in.hsd")
    fout = open("dftb_in.hsd", "wt")
    new_s = cell_formation(s,i)
    write('structure_'+str(i*0.001)+'.gen',new_s)

    for line in fin:
        if line=='    <<<\n':
            fout.write(line.replace('    <<<','    <<<"structure_%s.gen"'%(str(i*0.001))))
        elif line=='      atoms=orbital\n':
            for element in structure_elements:
                fout.write(line.replace('      atoms=orbital','                %s' %(elements_orbitals[element])))
        else:
            fout.write(line)
    fin.close()
    fout.close()

    fin = open("../run1.sh")
    fout = open("run.sh", "wt")
    for line in fin:
        if line=='#SBATCH -J\n':
            fout.write(line.replace('#SBATCH -J','#SBATCH -J por_anth_%s'%(str(i*0.001))))
        else:
            fout.write(line)
    fin.close()
    fout.close()
    cmd = 'sbatch run.sh'
    os.system(cmd)
    os.chdir(path)
