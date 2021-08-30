import xlwt
from xlwt import Workbook
from ase.io import read,write
from ase.calculators.dftb import Dftb
import sys,os
from math import sqrt,cos,radians
import numpy as np

########################################################################################################################################################################################
###################################################################   Band  structure G-X-M-G   ########################################################################################
########################################################################################################################################################################################

def band_structure_GXMG(structure):
    calc = Dftb(atoms=structure,
                run_manyDftb_steps=True,
                Hamiltonian_SCC='Yes',
                Hamiltonian_SCCTolerance=1.0e-6,
                Hamiltonian_MaxSCCIterations=1,
                Hamiltonian_SlaterKosterFiles_Prefix="/home/antonis/Desktop/parametrization/matsci-0-3/",
                Hamiltonian_MaxAngularMomentum_='',
                Hamiltonian_MaxAngularMomentum_C='"p"',
                Hamiltonian_MaxAngularMomentum_N='"p"',
                Hamiltonian_MaxAngularMomentum_H='"s"',
                Hamiltonian_MaxAngularMomentum_B='"p"',
                Hamiltonian_MaxAngularMomentum_O='"p"',
                Hamiltonian_KPointsAndWeights_='',
                Hamiltonian_KPointsAndWeights_KLines='{\n    1   0.0  0.0  0.0\n    50  0.5  0.0  0.0\n    50  0.5  0.5  0.0\n    50  0.0  0.0  0.0}')

    structure.set_calculator(calc)
    calc.calculate(structure)
    calc.get_potential_energy()

########################################################################################################################################################################################
###################################################################   Band  structure rhombic   ########################################################################################
########################################################################################################################################################################################

def band_structure_rhombic(structure):
    calc = Dftb(atoms=structure,
                run_manyDftb_steps=True,
                Hamiltonian_SCC='Yes',
                Hamiltonian_SCCTolerance=1.0e-6,
                Hamiltonian_MaxSCCIterations=1,
                Hamiltonian_SlaterKosterFiles_Prefix="/home/antonis/Desktop/parametrization/matsci-0-3/",
                Hamiltonian_MaxAngularMomentum_='',
                Hamiltonian_MaxAngularMomentum_C='"p"',
                Hamiltonian_MaxAngularMomentum_N='"p"',
                Hamiltonian_MaxAngularMomentum_H='"s"',
                Hamiltonian_MaxAngularMomentum_B='"p"',
                Hamiltonian_MaxAngularMomentum_O='"p"',
                Hamiltonian_KPointsAndWeights_='',
                Hamiltonian_KPointsAndWeights_KLines='{\n    1   0.0  0.0  0.0\n    50  0.475, -0.475  0.0\n    50  0.525, 0.475  0.0\n    50  0.5  0.5  0.0\n    50  0.0  0.0  0.0}')

    structure.set_calculator(calc)
    calc.calculate(structure)
    calc.get_potential_energy()

########################################################################################################################################################################################
###################################################################   Fermi level calculation   ########################################################################################
########################################################################################################################################################################################


def Fermi_level(structure):
    calc = Dftb(atoms=structure,
                Hamiltonian_SCC='Yes',
                Hamiltonian_SCCTolerance=1.0e-6,
                Hamiltonian_MaxSCCIterations=1000,
                Hamiltonian_Filling_ = '',
                Hamiltonian_Filling_Fermi_ = '',
                Hamiltonian_Filling_Fermi_Temperature = 0.02585,
                Hamiltonian_SlaterKosterFiles_Prefix="/home/antonis/Desktop/parametrization/matsci-0-3/",
                Hamiltonian_MaxAngularMomentum_='',
                Hamiltonian_MaxAngularMomentum_C='"p"',
                Hamiltonian_MaxAngularMomentum_H='"s"',
                Hamiltonian_MaxAngularMomentum_N='"p"',
                Hamiltonian_MaxAngularMomentum_O='"p"',
                Hamiltonian_MaxAngularMomentum_B='"p"',
                Hamiltonian_KPointsAndWeights_='',
                Hamiltonian_KPointsAndWeights_SupercellFolding='{\n    3 0 0\n    0 3 0\n    0 0 1\n    0.5 0.5 0.0}',
                Hamiltonian_Charge=0)

    structure.set_calculator(calc)
    calc.calculate(structure)

    with open('detailed.out') as f:
        lines = f.readlines()

    term="Fermi level"
    Fermi_level=[]
    for line in lines:
        if term in line: Fermi_level.append(line)
    E_fermi=round(float(Fermi_level[0][62:68]),2)

    return E_fermi

########################################################################################################################################################
###############################################################   polymers's  x-axis conversion   ######################################################
########################################################################################################################################################

def convert_square(l):
    points = [[0.0,0.0,0.0],        # Gamma
              [0.5,0.0,0.0],        # X
              [0.5,0.5,0.0],        # M
              [0.0,0.0,0.0]]        # Gamma

    n_points = [50,50,50]
    d = []

    for i in range(len(points)-1):
        s = 0
        for j in range(len(points[0])):
            s = s + abs(points[i][j]**2-points[i+1][j]**2)

        d.append(sqrt(s))
    x_points = []
    for i in range(len(n_points)):
        x_points.append(d[i]/n_points[i])

    x = [0]
    s = 0
    for i in range(len(n_points)):
        for j in range(n_points[i]):
            s = s + x_points[i]
            x.append(s)

    x_ticks = [x[0]]
    s = 0
    for i in range(len(n_points)):
        s=s+n_points[i]
        x_ticks.append(x[s])

    return x


def convert_rhombic(l):
    points = [[0.0  ,0.0   ,0.0],        # Gamma
              [0.475,-0.475,0.0],        # X
              [0.525, 0.475,0.0],        # A1
              [0.5  ,0.5   ,0.0],        # Y
              [0.0  ,0.0   ,0.0]]        # Gamma

    n_points = [50,50,50,50]
    d = []

    for i in range(len(points)-1):
        s = 0
        for j in range(len(points[0])):
            s = s + abs(points[i][j]**2-points[i+1][j]**2)

        d.append(sqrt(s))
    x_points = []
    for i in range(len(n_points)):
        x_points.append(d[i]/n_points[i])

    x = [0]
    s = 0
    for i in range(len(n_points)):
        for j in range(n_points[i]):
            s = s + x_points[i]
            x.append(s)

    x_ticks = [x[0]]
    s = 0
    for i in range(len(n_points)):
        s=s+n_points[i]
        x_ticks.append(x[s])

    return x

##########################################################################################################################################################
##################################################################   COF's  gap   ########################################################################
##########################################################################################################################################################

def gap():
    with open('band.out') as f:
        lines = f.readlines()

    term="KPT"
    find=[]
    for line in lines:
        if term in line:
            find.append(lines.index(line))
    electron_filling="0.00000"
    for i in range(find[0]+1,find[1]):
        if electron_filling in lines[i]:
            lumo=i
            break
    homo=lumo-1

    directory = 'Fermi'
    os.mkdir(directory)
    os.chdir(directory)
    s = read('../opttot.gen')
    E_f = Fermi_level(s)
    os.chdir('../')
    os.system('rm -r Fermi')

    data=np.loadtxt("band_tot.dat")
    gap_min_value=[]
    gap_index_min_value=[]
    for i in range(len(data[:,homo])):
        gap_finder=[]
        for j in range(len(data[:,lumo])):
            gap_finder.append(abs(data[:,lumo][i]-data[:,homo][j]))

        gap_min_value.append(min(gap_finder))
        gap_index_min_value.append([i+1,gap_finder.index(min(gap_finder))+1])

    Gap = round(min(gap_min_value),2)   # this is the Gap

    k_point_valence = gap_index_min_value[gap_min_value.index(min(gap_min_value))][0]  # this is the k-point for the gap of
    k_point_conduction = gap_index_min_value[gap_min_value.index(min(gap_min_value))][1]

    E_valence = round(data[:,lumo][k_point_valence-1],2)   # valence band
    E_conduction = round(data[:,homo][k_point_conduction-1],2)   # conduction band

    valence_band = []
    conduction_band = []
    for i in range(len(data)):
        valence_band.append(data[:,homo][i])
        conduction_band.append(data[:,lumo][i])
    return [Gap,E_conduction,E_valence,k_point_valence,k_point_conduction, E_f], conduction_band, valence_band

################################################################################################################################################################
################################################################################################################################################################


central_path = "/home/antonis/Desktop/comparison_3vs4_arms/porphyrin/"

sub_paths = ["anthracene/","phenazine/","quinoline/","quinoline_reversed/","quinoline_reversed_without_rings/","quinoline_without_rings/", "imine_CC_+1/","imine_CC_+2/","imine_NC_+2/","imine_NC_+2/"]

sub_paths = ["imine_NC_+1/","imine_NC_+2/"]

datas = []
for i in range(len(sub_paths)):

    print sub_paths[i]

    os.chdir(sub_paths[i] + "COF/")
    
    s = read("opttot.gen")

    if i>=6:
        try:
            band_structure_rhombic(s)
        except (RuntimeError, TypeError, NameError):
            pass
    else:
        try:
            band_structure_GXMG(s)
        except (RuntimeError, TypeError, NameError):
            pass
    
    cmd = "dp_bands band.out band"
    os.system(cmd)

    basics = gap()
    datas.append(basics)

    os.chdir(central_path)


valence = []
conduction = []
k_min = []
gap_COF = []
Fermi = []
for i in range(len(datas)):
    valence.append(datas[i][2])
    conduction.append(datas[i][1])
    k_min.append(datas[i][0][3])
    gap_COF.append(datas[i][0][0])
    Fermi.append(datas[i][0][5])

t = []
give_values = []
for i in range(len(sub_paths)):

    t3 = (-(conduction[i][0] + valence[i][0]) + (conduction[i][k_min[i]] + valence[i][k_min[i]]))/8

    e1_plus_e2 = ((conduction[i][0] + valence[i][0]) + (conduction[i][k_min[i]] + valence[i][k_min[i]]))/2
    e1_minus_e2 = (conduction[i][k_min[i]] - valence[i][k_min[i]]) + 4*t3

    t1 = np.sqrt(((conduction[i][0] - valence[i][0])**2 - (e1_minus_e2 + 4*t3)**2)/64)

    t.append(t1)

    E_1 = (e1_plus_e2 + e1_minus_e2)/2
    E_2 = (e1_plus_e2 - e1_minus_e2)/2



    param = [E_1, E_2, round(t1,4), round(t3,4)]
    print param

    give_values.append(param)


wb = Workbook()
sheet1 = wb.add_sheet('Sheet 1')
sheet1.write(0,0,"Linker")
sheet1.write(0,1,"Fermi level")
sheet1.write(0,2,"Gap")
sheet1.write(0,3,"LUMO")
sheet1.write(0,4,"HOMO")
sheet1.write(0,5,"t")
sheet1.write(0,6,"t*")

for i in range(len(give_values)):
    sheet1.write(1+i,0,sub_paths[i][:-1])
    sheet1.write(1+i,1,Fermi[i])
    sheet1.write(1+i,2,gap_COF[i])
    sheet1.write(1+i,3,give_values[i][0])
    sheet1.write(1+i,4,give_values[i][1])
    sheet1.write(1+i,5,give_values[i][2])
    sheet1.write(1+i,6,give_values[i][3])


wb.save('datas_COFs_check.xls')
