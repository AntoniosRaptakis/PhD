from ase.io import read,write
from ase.calculators.dftb import Dftb
import sys,os
from ase.eos import EquationOfState
from ase.visualize import view
from math import sqrt,cos,radians
import numpy as np

def structure_learning(structure):
    symbols=structure.get_chemical_symbols()
    structure_elements=[symbols[0]]
    for symbol in symbols:
        if symbol not in structure_elements:
            structure_elements.append(symbol)

    return structure_elements

def optimization(structure):
    calc = Dftb(atoms=structure,
                run_manyDftb_steps=True,
                Driver_='ConjugateGradient',
                Driver_MovedAtoms='1:-1',
                Driver_MaxForceComponent='1E-6',
                Driver_MaxSteps=50000,
                Driver_OutputPrefix = "opttot",
                Driver_AppendGeometries = 'Yes',
                Hamiltonian_SCC='Yes',
                Hamiltonian_SCCTolerance=1.0e-6,
                Hamiltonian_MaxSCCIterations=1000,
                Hamiltonian_SlaterKosterFiles_Prefix="/home/anra181a/matsci-0-3/",
                Hamiltonian_MaxAngularMomentum_='',
                Hamiltonian_MaxAngularMomentum_C='"p"',
                Hamiltonian_MaxAngularMomentum_H='"s"',
                Hamiltonian_MaxAngularMomentum_N='"p"',
                Hamiltonian_KPointsAndWeights_='',
                Hamiltonian_KPointsAndWeights_SupercellFolding='{\n    3 0 0\n    0 3 0\n    0 0 1\n    0.5 0.5 0.0}',
                Hamiltonian_Charge=0)

    structure.set_calculator(calc)
    calc.calculate(structure)
    return calc.get_potential_energy()

def Total_energy():
    with open('detailed.out') as f1:
        lines = f1.readlines()

    term = "Total energy"
    Energy= []
    for line in lines:
        if term in line: Energy.append(line)

    Total_energy = float(Energy[0][59:69])
    return Total_energy


def find_lattice_constant(v_0,structure):

    a1_init = structure.get_cell_lengths_and_angles()[0]
    a2_init = structure.get_cell_lengths_and_angles()[1]
    a3_init = structure.get_cell_lengths_and_angles()[2]
    structure.set_cell([a1_init,a2_init,30,s.get_cell_lengths_and_angles()[3],s.get_cell_lengths_and_angles()[4],s.get_cell_lengths_and_angles()[5]])
    V_init = structure.get_volume()
    factor = sqrt(v_0/V_init) - 1
    new_a1 = structure.get_cell_lengths_and_angles()[0]*(1+factor)
    new_a2 = structure.get_cell_lengths_and_angles()[1]*(1+factor)

    return new_a1, new_a2

s = read('anthracene.gen')
atoms = structure_learning(s)
path = "/scratch/ws/anra181a-SPECint/Bulk_modulus/porphyrin/anthracene/"
energy=[]
volume=[]
f=open('bulk_modulus.dat','w')
f.write("# The length of a1 vector before lattice optimization is %s" %s.get_cell_lengths_and_angles()[0]+"\n")
f.write("# The length of a2 vector before lattice optimization is %s" %s.get_cell_lengths_and_angles()[1]+"\n")
f.write("#   Stretch         |a1|       |a2|       Volume        Energy "+"\n")

for i in range(-8,10,2):
    path2 = "optimization_compress_stretch_%s/" %str(i*0.001)
    os.chdir(path2)
    s = read('opttot.gen')
    E = Total_energy()
    energy.append(E)
    volume.append(s.get_volume())
    a1 = s.get_cell_lengths_and_angles()[0]
    a2 = s.get_cell_lengths_and_angles()[1]
    os.chdir("../")
    f.write("{:10.3f}    {:10.3f}   {:10.3f}    {:10.3f}    {:10.3f}".format(i*0.001*100,a1,a2,s.get_volume(),E)+"\n")

eos = EquationOfState(volume, energy)
v0, e0, B = eos.fit()
s = read('anthracene.gen')
l1, l2 = find_lattice_constant(v0,s)
s.set_cell([l1,l2,s.get_cell_lengths_and_angles()[2],
           s.get_cell_lengths_and_angles()[3],
           s.get_cell_lengths_and_angles()[4],
           s.get_cell_lengths_and_angles()[5]],scale_atoms=True)
s.set_cell([l1,l2,30,s.get_cell_lengths_and_angles()[3],s.get_cell_lengths_and_angles()[4],s.get_cell_lengths_and_angles()[5]])
write('optimized_structure.gen',s)

f.write("# The optimized structure has volume %s and energy %s eV." %(v0,e0)+"\n")
f.write("# The length of the a1 vector is:%s." %s.get_cell_lengths_and_angles()[0]+"\n")
f.write("# The length of the a2 vector is:%s." %s.get_cell_lengths_and_angles()[1]+"\n")
f.write("# The length of the a3 vector is:%s." %s.get_cell_lengths_and_angles()[2]+"\n")
f.write("# The angle between a1 and a2 is:%s." %s.get_cell_lengths_and_angles()[3]+"\n")
f.write("# The angle between a1 and a3 is:%s." %s.get_cell_lengths_and_angles()[4]+"\n")
f.write("# The angle between a2 and a3 is:%s." %s.get_cell_lengths_and_angles()[5]+"\n")
f.close()

eos.plot('modulus.png')

best_structure=read('optimized_structure.gen')
directory_name="optimized_structure"
os.mkdir(os.path.join(path,directory_name))
os.chdir(directory_name)
optimization(best_structure)
