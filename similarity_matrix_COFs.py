import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.transforms as mtransforms
from math import sqrt
import matplotlib as mpl
import os,sys
from matplotlib import rc
import latex
from scipy.ndimage.filters import gaussian_filter1d
import statistics as stat


plt.rc('text', usetex=True)
plt.rc('font', family='Serif')
##############################################################################################################################################
#########################################################   COF's  x-axis conversion   #######################################################
##############################################################################################################################################
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

    top_valence = round(data[:,homo][k_point_valence-1],2)   # valence band
    valence_band = data[:,homo]
    pre_valence_band = data[:,homo-1]
    bottom_conduction = round(data[:,lumo][k_point_conduction-1],2)   # conduction band
    conduction_band = data[:,lumo]
    after_conduction_band = data[:,lumo+1]
    
    return [Gap,top_valence,valence_band,pre_valence_band,bottom_conduction,conduction_band,after_conduction_band,k_point_conduction,k_point_valence]



def ticks(a):
    tick=np.linspace(0,max(a),5)
    tick=[round(a,1) for a in tick]
    return tick

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
fig, ax = plt.subplots(figsize=(16,8))


path = "/home/antonis/Desktop/paper_electronic_properties/"

#core = ["structure_5C/rings_6/","some_others/2_phenazines_6_rings/","some_others/1_pyrene_linker/","saturation/1_pyrene_linker/+C/"]
core = ["structure_5C/rings_3/","structure_5N/rings_3/","structure_1/","structure_6/","structure_3/","structure_4","structure_7/rings_1/","structure_9/rings_1/"]

porphyrin, datas = [], []
os.chdir(path)
for i in range(len(core)):

    os.chdir("porphyrazine/" + core[i])

    porphyrin.append(np.loadtxt("band_tot.dat"))
    datas.append(gap())

    os.chdir(path)


gap, top_HOMO_porphyrin, HOMO_porphyrin, pre_HOMO_porphyrin, bottom_LUMO_porphyrin, LUMO_porphyrin, after_LUMO_porphyrin, gap_point, gap_upper = [], [], [], [], [], [], [], [], []
for i in range(len(datas)):

    gap.append(datas[i][0])
    top_HOMO_porphyrin.append(datas[i][1])
    HOMO_porphyrin.append(datas[i][2])
    pre_HOMO_porphyrin.append(datas[i][3])
    bottom_LUMO_porphyrin.append(datas[i][4])
    LUMO_porphyrin.append(datas[i][5])
    after_LUMO_porphyrin.append(datas[i][6])
    gap_point.append(datas[i][7])
    gap_upper.append(datas[i][8])

dE = []
for i in range(len(HOMO_porphyrin)):
    value = []
    for j in range(len(HOMO_porphyrin)):
        check_HOMO, check_LUMO = [], []
        for k in range(len(HOMO_porphyrin[i])):
            check_HOMO.append((HOMO_porphyrin[i][k] - HOMO_porphyrin[j][k])**2)
            check_LUMO.append((LUMO_porphyrin[i][k] - LUMO_porphyrin[j][k])**2)
        value.append((np.sum(np.array(check_HOMO)) + np.sum(np.array(check_LUMO)))/(2*len(HOMO_porphyrin[i])))
    dE.append(-(np.array(value) - min(value))/(max(value) - min(value)))       # for max value per row
#    dE.append(-np.array(value)/0.0259**2)                                     # for room temperature

similarity = np.exp(dE)
check_similarity = (similarity + similarity.T)/2

map_heat = ax.imshow(check_similarity, cmap='YlGnBu', interpolation='nearest')

matrix_ticks = np.round(np.linspace(np.min(check_similarity),np.max(check_similarity),5),2)

ax.set_xticklabels([])
ax.set_yticklabels([])

cbar = fig.colorbar(map_heat, ax=ax, ticks=matrix_ticks)
cbar.ax.set_yticklabels(matrix_ticks) 
cbar.ax.tick_params(labelsize=35)


for i in range(len(check_similarity)):
    for j in range(len(check_similarity[i])):
        text = ax.text(j, i, round(check_similarity[i, j],2),
                       ha="center", va="center", fontsize=25, color="w")

plt.show()
