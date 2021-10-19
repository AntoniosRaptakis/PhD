import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Ridge
from sklearn.linear_model import Lasso
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import r2_score
from yellowbrick.regressor import ResidualsPlot
from yellowbrick.regressor import PredictionError
from sklearn.linear_model import RidgeCV
from sklearn.linear_model import LassoCV
from yellowbrick.regressor import AlphaSelection
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.transforms as mtransforms
import matplotlib as mpl
import latex
import math
from mpl_toolkits.mplot3d import Axes3D
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
from matplotlib.colors import ListedColormap

plt.rc('font', family='Times New Roman')
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
df_porphyrin = pd.read_csv(r'porphyrin.csv')
df_porphyrin.dropna(inplace=True)
porphyrin_COF = pd.concat([df_porphyrin['Fermi level'], df_porphyrin['Gap'], df_porphyrin['HOMO.2'], df_porphyrin['LUMO.2'], df_porphyrin['t'], df_porphyrin['t*']], axis=1)
porphyrin_1D = pd.concat([df_porphyrin['Fermi level.1'], df_porphyrin['Gap.1'], df_porphyrin['HOMO.3'], df_porphyrin['LUMO.3'], df_porphyrin['t1'], df_porphyrin['t2']], axis=1)
porphyrin_monomer = df_porphyrin.drop(['Bridge', 'Nr. of rings', 'Fermi level', 'Gap', 'LUMO.2', 'HOMO.2', 't', 't*', 'Fermi level.1', 'Gap.1', 'LUMO.3', 'HOMO.3', 't1', 't2'], axis=1)
#porphyrin_monomer = df_porphyrin.drop(['Bridge', 'Nr. of rings','HOMO-2', 'LUMO+2', 'HOMO-LUMO','Energy SCC', 'Repulsive energy', 'Total energy', 'Nr. of electrons', 'HOMO-2.1', 'LUMO+2.1', 'HOMO-LUMO.1','Energy SCC.1', 'Repulsive energy.1', 'Total energy.1', 'Nr. of electrons.1','Fermi level', 'Gap', 'LUMO.2', 'HOMO.2', 't', 't*', 'Fermi level.1', 'Gap.1', 'LUMO.3', 'HOMO.3', 't1', 't2'], axis=1)


df_pyrene = pd.read_csv(r'pyrene.csv')
df_pyrene.dropna(inplace=True)
pyrene_COF = pd.concat([df_pyrene['Fermi level'], df_pyrene['Gap'], df_pyrene['HOMO.2'], df_pyrene['LUMO.2'], df_pyrene['t'], df_pyrene['t*']], axis=1)
pyrene_1D = pd.concat([df_pyrene['Fermi level.1'], df_pyrene['Gap.1'], df_pyrene['HOMO.3'], df_pyrene['LUMO.3'], df_pyrene['t1'], df_pyrene['t2']], axis=1)
pyrene_monomer = df_pyrene.drop(['Bridge', 'Nr. of rings', 'Fermi level', 'Gap', 'LUMO.2', 'HOMO.2', 't', 't*', 'Fermi level.1', 'Gap.1', 'LUMO.3', 'HOMO.3', 't1', 't2'], axis=1)
#pyrene_monomer = df_pyrene.drop(['Bridge', 'Nr. of rings','HOMO-2', 'LUMO+2', 'HOMO-LUMO','Energy SCC', 'Repulsive energy', 'Total energy', 'Nr. of electrons', 'HOMO-2.1', 'LUMO+2.1', 'HOMO-LUMO.1','Energy SCC.1', 'Repulsive energy.1', 'Total energy.1', 'Nr. of electrons.1','Fermi level', 'Gap', 'LUMO.2', 'HOMO.2', 't', 't*', 'Fermi level.1', 'Gap.1', 'LUMO.3', 'HOMO.3', 't1', 't2'], axis=1)

X_porphyrin_train, X_porphyrin_test, y_porphyrin_train, y_porphyrin_test = train_test_split(porphyrin_monomer, porphyrin_COF['Gap'], test_size=0.25, random_state=40)
X_pyrene_train, X_pyrene_test, y_pyrene_train, y_pyrene_test = train_test_split(pyrene_monomer, pyrene_COF['Gap'], test_size=0.25, random_state=40)


linear_por, linear_pyr = LinearRegression(), LinearRegression()
#ridge = Ridge(alpha=0.22)
#lasso = Lasso(alpha=0.007)

linear_por.fit(X_porphyrin_train,y_porphyrin_train)
linear_pyr.fit(X_pyrene_train,y_pyrene_train)

columns = porphyrin_monomer.columns[0:12]
porphyrin_core_coef = linear_por.coef_[0:12]/linear_por.coef_.sum()
porphyrin_bridge_coef = linear_por.coef_[12:]/linear_por.coef_.sum()
pyrene_core_coef = linear_pyr.coef_[0:12]/linear_pyr.coef_.sum()
pyrene_bridge_coef = linear_pyr.coef_[12:]/linear_pyr.coef_.sum()
##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
fig = plt.figure(figsize=(8,6))
gs = gridspec.GridSpec(1,1)
ax1 = plt.subplot(gs[0])

ax1.scatter(porphyrin_COF['Gap'], linear_por.predict(porphyrin_monomer), color = "tab:blue", label='Porphyrin based COFs')
ax1.scatter(pyrene_COF['Gap'], linear_pyr.predict(pyrene_monomer), color = "tab:olive", label='Pyrene based COFs')
ax1.plot([0,2.5],[0,2.5], color='black')
ax1.grid(False)
ax1.legend(fontsize=18)
ax1.set_xlim(0,2.5)
ax1.set_ylim(0,2.5)

ticks = [0,0.5,1.0,1.5,2.0,2.5]
ax1.set_xticks(ticks)
ax1.set_xticklabels(ticks, fontsize=15)
ax1.set_yticks(ticks)
ax1.set_yticklabels(ticks, fontsize=15)

ax1.set_xlabel(r'$\mathrm{Calculated\ with\ DFTB\ gap\ (eV)}$', labelpad=10, fontsize=25,weight='bold', style='oblique')
ax1.set_ylabel(r'$\mathrm{Predicted\ gap\ (eV)}$', labelpad=10, fontsize=25,weight='bold', style='oblique')
#print("R^2 porphyrin: %s" %r2_score(porphyrin_COF['Gap'], linear_por.predict(porphyrin_monomer)))
#print("R^2 pyrene: %s" %r2_score(pyrene_COF['Gap'], linear_pyr.predict(pyrene_monomer)))
#print("MSE porphyrin: %s" %mean_squared_error(porphyrin_COF['Gap'], linear_por.predict(porphyrin_monomer)))
#print("MSE pyrene: %s" %mean_squared_error(pyrene_COF['Gap'], linear_pyr.predict(pyrene_monomer)))

print(mean_absolute_error(porphyrin_COF['Gap'], linear_por.predict(porphyrin_monomer)))

plt.show()
##################################################################################################################################################################
fig = plt.figure(figsize=(8,6))
gs = gridspec.GridSpec(1,2)
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])

y = np.arange(len(columns))
width = 0.1

ax1.barh(y - width/2, porphyrin_core_coef, width, color="tab:blue", label="Porphyrin based COFs")
ax1.barh(y + width/2, pyrene_core_coef, width, color="tab:olive", label='Pyrene based COFs')
ax1.legend(fontsize=14)
ax1.set_yticks(y)
ax1.set_yticklabels(columns, fontsize=15)
x_ticks = [-0.0041, 0.0, 0.0024, 0.0059]
ax1.set_xticks(x_ticks)
ax1.set_xticklabels(x_ticks, fontsize=15)
ax1.set_title(r'$\mathrm{Core}$', fontsize=20)

ax2.barh(y - width/2, porphyrin_bridge_coef, width, color="tab:blue")
ax2.barh(y + width/2, pyrene_bridge_coef, width, color="tab:olive")
ax2.yaxis.tick_right()
ax2.set_yticks(y)
ax2.set_yticklabels(columns, fontsize=15)
x_ticks = [-0.11, 0, 0.35, 0.70, 1.10]
ax2.set_xticks(x_ticks)
ax2.set_xticklabels(x_ticks, fontsize=15)
ax2.set_title(r'$\mathrm{Bridge}$', fontsize=20)


plt.show()
##################################################################################################################################################################
fig = plt.figure(figsize=(8,6))
gs = gridspec.GridSpec(3,2)
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])
ax3 = plt.subplot(gs[2])
ax4 = plt.subplot(gs[3])
ax5 = plt.subplot(gs[4])
ax6 = plt.subplot(gs[5])
ax = [ax1, ax2, ax3, ax4, ax5, ax6]

porphyrin = [porphyrin_COF['Gap'], linear_por.predict(porphyrin_monomer)]
pyrene = [pyrene_COF['Gap'], linear_pyr.predict(pyrene_monomer)]

ax1.hist(df_porphyrin['HOMO-LUMO'], color="tab:olive", alpha=0.4, label="Bridge HL gap")
ax1.hist(porphyrin[0], color="tab:blue", alpha=0.4, label="Calculated gap with DFTB")
ax1.hist(df_porphyrin['Gap.1'], color="tab:purple", alpha=0.4, label="1D-polymer")
ax1.axvline(x = df_porphyrin['HOMO-LUMO.1'][0], color="tab:red", label="Core HL gap")
ax1.legend(fontsize=15)
ax1.set_title("Porphyrin", fontsize=20)

ax2.hist(df_pyrene['HOMO-LUMO'], color="tab:olive", alpha=0.4)
ax2.hist(pyrene[0], color="tab:blue", alpha=0.4)
ax2.axvline(x = df_pyrene['HOMO-LUMO.1'][0], color="tab:red")
ax2.hist(df_pyrene['Gap.1'], color="tab:purple", alpha=0.4)
ax2.set_title("Pyrene", fontsize=20)

ax3.hist(porphyrin[0], color="tab:blue", alpha=0.4, label="Calculated gap with DFTB")
ax3.hist(porphyrin[1], color="tab:orange", alpha=0.4, label="Predicted gap")
ax3.legend(fontsize=15)

ax4.hist(pyrene[0], color="tab:blue", alpha=0.4)
ax4.hist(pyrene[1], color="tab:orange", alpha=0.4)

ax5.boxplot(porphyrin, vert=False)
ax5.set_yticks([1,2])
ax5.set_yticklabels(["Calculated gap", "Predicted gap"], rotation=30, fontsize=15)

ax6.boxplot(pyrene, vert=False)
ax6.set_yticks([1,2])
ax6.set_yticklabels(["Calculated gap", "Predicted gap"], rotation=30, fontsize=15)

x_ticks = [0,1,2,3,4]
for i in range(len(ax)):
    ax[i].set_xlim(0,4.5)
    ax[i].set_xticks(x_ticks)
    ax[i].set_xticklabels(x_ticks, fontsize=13)
    if i!=4 and i!=5:
        ax[i].set_ylim(0,14)

plt.show()
##################################################################################################################################################################
##################################################################################################################################################################
fig = plt.figure()
ax = Axes3D(fig)

ax.scatter(df_porphyrin['HOMO-LUMO'], df_porphyrin['HOMO-LUMO.1'], df_porphyrin['Gap'], marker = 'o', color = "tab:blue")
ax.scatter(df_pyrene['HOMO-LUMO'], df_pyrene['HOMO-LUMO.1'], df_pyrene['Gap'], marker = 'o', color = "tab:olive")
ax.set_xlabel("Bridge HL gap (eV)", fontsize=20)
ax.set_ylabel("Core HL gap (eV)", fontsize=20)
ax.set_zlabel("COF gap (eV)", fontsize=20)

plt.show()
