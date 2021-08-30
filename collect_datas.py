import xlwt 
from xlwt import Workbook 
import os,sys
import numpy as np

 
wb = Workbook() 

sheet1 = wb.add_sheet('Sheet 1') 
sheet1.write(0, 2, '1D polymer')
sheet1.write(0, 6, 'COF')
sheet1.write(1, 0, 'linker')
sheet1.write(1, 1, 'number of rings')
sheet1.write(1, 2, 'gap')
sheet1.write(1, 3, 'E1')
sheet1.write(1, 4, 'E2')
sheet1.write(1, 5, 't1')
sheet1.write(1, 6, 't2')
sheet1.write(1, 7, 'gap')
sheet1.write(1, 8, 'E1')
sheet1.write(1, 9, 'E2')
sheet1.write(1, 10, 't')
sheet1.write(1, 11, 't*')

central_path = "/home/antonis/Desktop/porphyrin_1D_polymer/"
path_linkers = ["simple_phenyls/", "anthracene/","phenazine/","phthalimidophthalimide/","catecholato_diboron/"]

data_polymer, data_COF = [], []
for i in range(len(path_linkers)):
    
    os.chdir(path_linkers[i])

    data_polymer.append(np.loadtxt("tight_binding_imines_polymer.dat"))
    data_COF.append(np.loadtxt("tight_binding_imines_cof.dat"))
    
    os.chdir(central_path)

for i in range(len(data_polymer)):
    
    for j in range(len(data_polymer[i])):

        sheet1.write(2 + i + j, 1, data_polymer[i][j][0])
        sheet1.write(2 + i + j, 2, abs(data_polymer[i][j][1]-data_polymer[i][j][2]))
        sheet1.write(2 + i + j, 3, data_polymer[i][j][1])
        sheet1.write(2 + i + j, 4, data_polymer[i][j][2])
        sheet1.write(2 + i + j, 5, data_polymer[i][j][3])
        sheet1.write(2 + i + j, 6, data_polymer[i][j][4])

    for j in range(len(data_COF[i])):

        sheet1.write(2 + i + j, 7, abs(data_COF[i][j][1]-data_COF[i][j][2]))
        sheet1.write(2 + i + j, 8, data_COF[i][j][1])
        sheet1.write(2 + i + j, 9, data_COF[i][j][2])
        sheet1.write(2 + i + j, 10, data_COF[i][j][3])
        sheet1.write(2 + i + j, 11, data_COF[i][j][4])

wb.save('tight_binding_imines.xls') 
