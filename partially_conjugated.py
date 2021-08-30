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


plt.rc('text', usetex=True)
plt.rc('font', family='Serif')
##############################################################################################################################################
#########################################################   COF's  x-axis conversion   #######################################################
##############################################################################################################################################
def convert_DFTB(l):
    points = [[0.0,0.0,0.0],
              [0.5,0.0,0.0],
              [0.5,0.5,0.0],
              [0.0,0.0,0.0]]

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

    return x,x_ticks

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

    E_valence = round(data[:,lumo][k_point_valence-1],2)   # valence band
    E_conduction = round(data[:,homo][k_point_conduction-1],2)   # conduction band

    return [Gap,E_conduction,E_valence,k_point_conduction,k_point_valence]



def ticks(a):
    tick=np.linspace(0,max(a),5)
    tick=[round(a,1) for a in tick]
    return tick

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
fig = plt.figure(figsize=(8,6))
gs = gridspec.GridSpec(4, 4)
ax1 = plt.subplot(gs[0:3,0])
ax2 = plt.subplot(gs[0:3,1])
ax3 = plt.subplot(gs[0:3,2])
ax4 = plt.subplot(gs[0:3,3])
ax = [ax1,ax2,ax3,ax4] 
fig.subplots_adjust(wspace=0, hspace=0.5)

path = "/home/antonis/Desktop/paper_electronic_properties/"

core = ["structure_1/","structure_6/","structure_3/","structure_4/"]

porphyrin, datas = [], []
os.chdir(path)
for i in range(len(core)):

    os.chdir("porphyrin/" + core[i])

    porphyrin.append(np.loadtxt("band_tot.dat"))
    datas.append(gap())

    os.chdir(path)


HOMO_porphyrin, LUMO, gap_point, gap_upper = [], [], [], []
for i in range(len(datas)):
    HOMO_porphyrin.append(datas[i][1])
    LUMO.append(datas[i][2])
    gap_point.append(datas[i][3])
    gap_upper.append(datas[i][0])


porphyrazine, datas = [], []
os.chdir(path)
for i in range(len(core)):

    os.chdir("porphyrazine/" + core[i])

    porphyrazine.append(np.loadtxt("band_tot.dat"))
    datas.append(gap())

    os.chdir(path)


HOMO_porphyrazine, LUMO, gap_point, gap_upper = [], [], [], []
for i in range(len(datas)):
    HOMO_porphyrazine.append(datas[i][1])
    LUMO.append(datas[i][2])
    gap_point.append(datas[i][3])
    gap_upper.append(datas[i][0])


x_ticklabels = ["$\Gamma$","X","M","$\Gamma$"]
for i in range(len(ax)):
    x_axis ,x_ticks = convert_DFTB(len(porphyrazine[i][0]))
    for j in range(len(porphyrazine[i][0])):
        if j!=0:
            ax[i].plot(x_axis, porphyrin[i][:,j] - HOMO_porphyrin[i], color="tab:blue", label = r'$\mathrm{H2-TBPor}$')
            ax[i].plot(x_axis, porphyrazine[i][:,j] - HOMO_porphyrazine[i], color="tab:orange", label=r'$\mathrm{H2-Phthal}$')

    ax[i].set_ylim(-2,2)
    ax[i].set_xlim(min(x_axis),max(x_axis))
    ax[i].set_xticks(x_ticks)
    ax[i].set_xticklabels(x_ticklabels,fontsize=40)
    y = [-2,-1,0,1,2]
    if i==0:
        ax[i].set_yticks(y)
        ax[i].set_yticklabels(y,fontsize=40)
        ax[i].set_ylabel(r'$\mathrm{E\ -\ E_{VB}\ (eV)}$', labelpad=20, fontsize=45, weight='bold', style='oblique')
    else:
        ax[i].set_yticks([])

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
fig = plt.figure(figsize=(8,6))
gs = gridspec.GridSpec(2, 1)
ax1 = plt.subplot(gs[0])
fig.subplots_adjust(wspace=0, hspace=0.5)

df = pd.read_excel (r'datas_square.xls')

porphyrin = df["Gap"].tolist()[:-1]
porphyrazine = df["Gap.1"].tolist()[:-1]
phthalocyanine = df["Gap.2"].tolist()[:-1]

gap_porphyrin = porphyrin[:4]
gap_porphyrazine = porphyrazine[:4]


x = np.arange(len(gap_porphyrin))
width = 0.1

ax1.bar(x-width/2, gap_porphyrin, width,label=r'$\mathrm{TBPor}$')
ax1.bar(x+width/2, gap_porphyrazine, width,label=r'$\mathrm{Phthal}$')
ax1.legend(loc='upper center',frameon=False,prop={'size':25})
y_ticks = ticks(gap_porphyrin)
ax1.set_yticks(y_ticks)
ax1.set_yticklabels(y_ticks,fontsize=40)
ax1.set_xticks(x)
ax1.set_xticklabels([])
ax1.set_ylabel(r'$\mathrm{Gap\ (eV)}$',fontsize=45, labelpad=20, weight='bold', style='oblique')

plt.show()
