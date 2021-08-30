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


def ticks(a):
    tick=np.linspace(0,max(a),5)
    tick=[round(a,1) for a in tick]
    return tick

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
fig, ax = plt.subplots(figsize=(16,8))


path = "/home/antonis/Desktop/paper_electronic_properties/"

core = ["structure_5C/rings_3/","structure_5N/rings_3/","structure_1/","structure_6/","structure_3/","structure_4","structure_7/rings_1/","structure_9/rings_1/"]

porphyrin, datas = [], []
os.chdir(path)
for i in range(len(core)):

    os.chdir("porphyrin/" + core[i])

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


distance_HOMO, distance_LUMO, gap_difference = [], [], []
for k in range(len(HOMO_porphyrin)):
    for_HOMO, for_LUMO, for_gap = [], [], []
    for i in range(len(HOMO_porphyrin)):
        check_HOMO, check_LUMO, check_gap = [], [], []
        for j in range(len(HOMO_porphyrin[k])):
            d = ((HOMO_porphyrin[k][j] - top_HOMO_porphyrin[k]) - (HOMO_porphyrin[i][j] - top_HOMO_porphyrin[i]))**2
            check_HOMO.append(d)

            d = ((bottom_LUMO_porphyrin[k] - LUMO_porphyrin[k][j]) - (bottom_LUMO_porphyrin[i] - LUMO_porphyrin[i][j]))**2
            check_LUMO.append(d)

            d = ((HOMO_porphyrin[k][j] - LUMO_porphyrin[k][j]) - (HOMO_porphyrin[i][j] - LUMO_porphyrin[i][j]))**2
            check_gap.append(d)

        for_HOMO.append(stat.stdev(check_HOMO))
        for_LUMO.append(stat.stdev(check_LUMO))
        for_gap.append(stat.stdev(check_gap))

    distance_HOMO.append(for_HOMO)
    distance_LUMO.append(for_LUMO)
    gap_difference.append(for_gap)


how_similar = np.exp(-np.array(distance_HOMO)) + np.exp(-np.array(distance_LUMO)) + np.exp(-np.array(gap_difference))
for k in range(len(how_similar)):
    similar_min, similar_max = min(how_similar[k]), max(how_similar[k])
    for i,val in enumerate(how_similar[k]):
        how_similar[k][i] = (val-similar_min) / (similar_max-similar_min)


check_similarity = how_similar.copy()
for i in range(len(how_similar)):
    for j in range(len(how_similar[0])):
        if i>j:
            if how_similar[i][j]!=how_similar[j][i]:
                check_similarity[i][j]=(how_similar[j][i] + how_similar[i][j])/2

for i in range(len(how_similar)):
    for j in range(len(how_similar[0])):
        if j>i:
            if how_similar[j][i]!=how_similar[i][j]:
                check_similarity[j][i]=check_similarity[i][j]

map_heat = ax.imshow(check_similarity, cmap='YlGnBu', interpolation='nearest')

ax.set_xticklabels([])
ax.set_yticklabels([])


ticks=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
cbar = fig.colorbar(map_heat, ax=ax, ticks=ticks)
cbar.ax.set_yticklabels(ticks) 
cbar.ax.tick_params(labelsize=35)


for i in range(len(check_similarity)):
    for j in range(len(check_similarity[i])):
        text = ax.text(j, i, round(check_similarity[i, j],2),
                       ha="center", va="center", fontsize=25, color="w")

plt.show()
