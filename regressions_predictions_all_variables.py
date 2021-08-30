import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression, Ridge, Lasso, RidgeCV, LassoCV
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
from yellowbrick.regressor import ResidualsPlot
from yellowbrick.regressor import PredictionError
from yellowbrick.regressor import AlphaSelection
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.transforms as mtransforms
import matplotlib as mpl
import latex

plt.rc('font', family='Times New Roman')
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
df = pd.read_excel(r'tight_binding_porphyrin.xlsx')

COF_LUMO = df["LUMO.1"]
COF_LUMO.dropna(inplace=True)
COF_HOMO = df["HOMO.1"]
COF_HOMO.dropna(inplace=True)
COF_gap = df["Gap.1"]
COF_gap.dropna(inplace=True)

polymer_HOMO = df["HOMO"]
polymer_LUMO = df["LUMO"]
polymer_coupling = df["t1"]
polymer_fermi = df["Fermi level"]

linker_HOMO = df["HOMO.2"]
linker_LUMO = df["LUMO.2"]
linker_fermi =  df["Fermi level.2"]

core_HOMO = df["HOMO.3"]
core_LUMO = df["LUMO.3"]
core_fermi = df["Fermi level.3"]

dimer_HOMO = df["HOMO.4"]
dimer_LUMO = df["LUMO.4"]
dimer_fermi = df["Fermi level.4"]

df_new = pd.concat([polymer_HOMO,polymer_LUMO,polymer_coupling,polymer_fermi,linker_HOMO,linker_LUMO,linker_fermi,core_HOMO,core_LUMO,core_fermi,dimer_HOMO,dimer_LUMO,dimer_fermi],axis=1)
df_new.dropna(inplace=True)

X_train, X_test, y_train, y_test = train_test_split(df_new,COF_HOMO,test_size=0.25,random_state=40)
#####################################################################################################################################################################
#####################################################################################################################################################################
############################################################################    HOMO     ############################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
############################################################################    Error    ############################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
fig = plt.figure(figsize=(8,6))
gs = gridspec.GridSpec(3,3)
ax = plt.subplot(gs[0])
fig.subplots_adjust(hspace=0.1,wspace=0.1)

alphas = np.logspace(-10,1,200)
visualizer=AlphaSelection(RidgeCV(alphas=alphas),ax=ax)
visualizer.fit(df_new,COF_HOMO)
visualizer.ax.grid(False)
visualizer.ax.legend(loc='lower right',frameon=False,prop={'size':20})
visualizer.ax.set_ylabel(r'$\mathrm{Error\ or\ score}$', labelpad=10, fontsize=25)
visualizer.ax.set_xlabel(r'$\mathrm{alpha}$', labelpad=10, fontsize=25)
visualizer.show()
#####################################################################################################################################################################
#####################################################################################################################################################################
######################################################################    Regressions     ###########################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#  Linear regression
linear_HOMO = LinearRegression()
visualizer1 = PredictionError(linear_HOMO,ax=ax)
visualizer1.fit(X_train,y_train)
visualizer1.score(X_test,y_test)

#  Ridge
linear_ridge_HOMO = Ridge(alpha=0.006)
visualizer1 = PredictionError(linear_ridge_HOMO,ax=ax)
visualizer1.fit(X_train,y_train)
visualizer1.score(X_test,y_test)

#  Lasso
linear_lasso_HOMO = Lasso(alpha=0.0015)
visualizer1 = PredictionError(linear_lasso_HOMO,ax=ax)
visualizer1.fit(X_train,y_train)
visualizer1.score(X_test,y_test)
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
############################################################################    LUMO     ############################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
############################################################################    Error    ############################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
X_train, X_test, y_train, y_test = train_test_split(df_new,COF_LUMO,test_size=0.25,random_state=40)

fig = plt.figure(figsize=(8,6))
gs = gridspec.GridSpec(3,3)
ax = plt.subplot(gs[0])
fig.subplots_adjust(hspace=0.1,wspace=0.1)

alphas = np.logspace(-10,1,200)
visualizer=AlphaSelection(RidgeCV(alphas=alphas),ax=ax)
visualizer.fit(df_new,COF_LUMO)
visualizer.ax.grid(False)
visualizer.ax.legend(loc='lower right',frameon=False,prop={'size':20})
visualizer.ax.set_ylabel(r'$\mathrm{Error\ or\ score}$', labelpad=10, fontsize=25)
visualizer.ax.set_xlabel(r'$\mathrm{alpha}$', labelpad=10, fontsize=25)
visualizer.show()
#####################################################################################################################################################################
#####################################################################################################################################################################
######################################################################    Regressions     ###########################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#  Linear regression
linear_LUMO = LinearRegression()
visualizer1 = PredictionError(linear_LUMO,ax=ax)
visualizer1.fit(X_train,y_train)
visualizer1.score(X_test,y_test)

#  Ridge
linear_ridge_LUMO = Ridge(alpha=0.001)
visualizer1 = PredictionError(linear_ridge_LUMO,ax=ax)
visualizer1.fit(X_train,y_train)
visualizer1.score(X_test,y_test)

#  Lasso
linear_lasso_LUMO = Lasso(alpha=0.0015)
visualizer1 = PredictionError(linear_lasso_LUMO,ax=ax)
visualizer1.fit(X_train,y_train)
visualizer1.score(X_test,y_test)
#####################################################################################################################################################################
#####################################################################################################################################################################
###################################################################    Coefficients    ##############################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
def ticks(a):
    tick=np.linspace(min(a),max(a),10)
    tick=[round(a,1) for a in tick]
    return tick

coefficients_linear_HOMO = linear_HOMO.coef_/linear_HOMO.coef_.sum()
coefficients_ridge_HOMO = linear_ridge_HOMO.coef_/linear_ridge_HOMO.coef_.sum()
coefficients_lasso_HOMO = linear_lasso_HOMO.coef_/linear_lasso_HOMO.coef_.sum()

coefficients_linear_LUMO = linear_LUMO.coef_/linear_LUMO.coef_.sum()
coefficients_ridge_LUMO = linear_ridge_LUMO.coef_/linear_ridge_LUMO.coef_.sum()
coefficients_lasso_LUMO = linear_lasso_LUMO.coef_/linear_lasso_LUMO.coef_.sum()

fig = plt.figure(figsize=(8,6))
gs = gridspec.GridSpec(1,2)
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])
fig.subplots_adjust(hspace=0.,wspace=0.)

y = np.arange(len(coefficients_linear_LUMO))
width = 0.1

ax1.barh(y-width,coefficients_linear_HOMO,width,color="tab:blue",label="Linear regression")
ax1.barh(y,coefficients_ridge_HOMO,width,color="tab:orange",label="Ridge")
ax1.barh(y+width,coefficients_lasso_HOMO,width,color="tab:green",label="Lasso")
ax1.legend(loc='upper right',frameon=False,prop={'size':20})
y_ticklabels = ["Polymer:HOMO","Polymer:LUMO","Polymer:t1","Polymer:Fermi","Linker:HOMO","Linker:LUMO","Linker:Fermi","Core:HOMO","Core:LUMO","Core:Fermi","Dimer:HOMO","Dimer:LUMO","Dimer:Fermi"]
ax1.set_yticks(y)
ax1.set_yticklabels(y_ticklabels,fontsize=20)
x_ticks = ticks(coefficients_linear_HOMO)
ax1.set_xticks(x_ticks)
ax1.set_xticklabels(x_ticks,fontsize=18)
ax1.set_title(r'$\mathrm{HOMO}$', y=1.07, fontsize=25)


ax2.barh(y-width,coefficients_linear_LUMO,width,color="tab:blue",label="Linear regression")
ax2.barh(y,coefficients_ridge_LUMO,width,color="tab:orange",label="Ridge")
ax2.barh(y+width,coefficients_lasso_LUMO,width,color="tab:green",label="Lasso")
ax2.set_yticks(y)
ax2.set_yticklabels([])
x_ticks = ticks(coefficients_linear_LUMO)
ax2.set_xticks(x_ticks)
ax2.set_xticklabels(x_ticks,fontsize=18)
ax2.set_title(r'$\mathrm{LUMO}$', y=1.07, fontsize=25)


fig.text(0.45, 0.03, r'$\mathrm{Coefficients}$', fontsize=25, va='center')


plt.show()
#####################################################################################################################################################################
#####################################################################################################################################################################
###################################################################    Predictions    ###############################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
def ticks(a):
    tick=np.linspace(min(a),max(a),10)
    tick=[round(a,1) for a in tick]
    return tick

fig = plt.figure(figsize=(8,6))
gs = gridspec.GridSpec(3,3)
ax = plt.subplot(gs[0])
fig.subplots_adjust(hspace=0.1,wspace=0.1)

ax.scatter(COF_gap,np.dot(df_new.values, coefficients_linear_LUMO) - np.dot(df_new.values, coefficients_linear_HOMO),color="tab:blue")
ax.scatter(COF_gap,np.dot(df_new.values, coefficients_ridge_LUMO) - np.dot(df_new.values, coefficients_ridge_HOMO),color="tab:orange")
ax.scatter(COF_gap,np.dot(df_new.values, coefficients_lasso_LUMO) - np.dot(df_new.values, coefficients_lasso_HOMO),color="tab:green")

ax.set_ylabel(r'$\mathrm{Predicted\ gap\ (eV)}$', labelpad=10, fontsize=25, weight='bold', style='oblique')
ax.set_xlabel(r'$\mathrm{Actual\ Gap\ (eV)}$', labelpad=10, fontsize=25, weight='bold', style='oblique')

plt.show()
