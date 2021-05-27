# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 15:51:35 2020

@author: linigodelacruz
"""
import numpy as np
from collections import defaultdict
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import norm
import os
os.chdir('D:/Users/Thomas/Studie/MEP/python-modules-for-bioinformatic-analyses/src')

path = 'D:/Users/Thomas/Studie/MEP/python-modules-for-bioinformatic-analyses/src/python_modules'
import sys
sys.path.insert(0, path) #hack to add module to path. otherwise it won't be found. Should find a better way to do this.

# set PYTHONPATH=%PYTHONPATH%;D:\Users\Thomas\Studie\MEP\python-modules-for-bioinformatic-analyses\src
from module_common_measures import common_go
from module_common_measures import common_partners
from module_common_measures import plotCommonHist
from module_common_measures import plotCommonScatter
# from python_modules.module_common_measures import common_partners_T

#%% getting the data

data_go=pd.read_excel('../datasets/slim-goterms-filtered-data.xlsx')

data_go.columns=['Gene','gene-id','go-aspect','go-term','go-id','feature-type']

data=pd.read_excel('../datasets/data-BioGrid-Yeast-doubled.xlsx')

#%% query
# query = pd.read_excel('../datasets/genesProteinFolding_SGD_amigo.xlsx')
# query = np.unique(query).tolist()
query = ['BEM1']

#%% Calling the function common_partners

common_partners_data=common_partners(query,data)

common_go_data=common_go(data_go=data_go,data_common_partners=common_partners_data)

#%% testing moduleT functions


# query = ['BEM1', 'BEM2']
# testingPartners = common_partners_T(query,data)


#%% Postprocessing the data

sl=common_partners_data[common_partners_data['Type']=='Synthetic Lethality']
pg=common_partners_data[common_partners_data['Type']=='Positive Genetic']
ng=common_partners_data[common_partners_data['Type']=='Negative Genetic']

common_go_data.loc[sl.index,'score']='SL'
common_go_data.loc[pg.index,'score']='PG'
common_go_data.loc[ng.index,'score']='NG'

common_go_data.loc[sl.index,'common_interactors']=common_partners_data.loc[sl.index,'fraction-of-common-partners']
common_go_data.loc[pg.index,'common_interactors']=common_partners_data.loc[pg.index,'fraction-of-common-partners']
common_go_data.loc[ng.index,'common_interactors']=common_partners_data.loc[ng.index,'fraction-of-common-partners']

#%% Messing around: Correlation coefficient etc.
# covar = np.cov(ng['fraction-of-common-partners'].astype(float),ng['fraction-of-common-partners'].astype(float))
# corcoeff = np.corrcoef(ng['fraction-of-common-partners'].astype(float),ng['fraction-of-common-partners'].astype(float))
# testdata = np.random.normal(loc=5.0, scale=2.0, size=1000)
# mean,std=norm.fit(data)
 
#%% Working pairplot ()
gene = 'BEM1'
# for gene in query:
tmp = common_go_data[common_go_data['query']==gene]
sns.set(style="ticks", color_codes=True)
# plot=sns.pairplot(tmp,vars=['fraction-of-common-go','common_interactors'],hue='score',hue_order=['NG','SL','PG'],
                  # diag_kind="hist", diag_kws = {'bins':int(np.ceil(np.sqrt(tmp.shape[0])))},corner=True)
plot=sns.pairplot(tmp,vars=['fraction-of-common-go','common_interactors'],hue='score',hue_order=['NG','SL','PG'],
                  corner=True, diag_kws = {'bw' : 10, 'kernel' : 'tri'})
# plt.title(query[0])
plot.fig.suptitle(gene)
# plot.savefig('../output_images/Testing/common-go-terms-of-'+ gene +'-based-on-their-type.png',dpi=300,format='png',transparent=True)


#%% Working Combined plot

# gene = 'BEM1'
# bins = int(np.ceil(np.sqrt(tmp.shape[0])))
sns.set(style="ticks", color_codes=True)
# plot=sns.pairplot(common_go_data,vars=['fraction-of-common-go','common_interactors'],hue='score',hue_order=['NG','PG','SL'],
                  # diag_kind="hist",diag_kws = {'bins':int(np.ceil(np.sqrt(common_go_data.shape[0])))},corner=True)
plot=sns.pairplot(common_go_data,vars=['fraction-of-common-go','common_interactors'],hue='score',hue_order=['NG','SL','PG'])
# plt.title(query[0])
plot.fig.suptitle('Protein Folding Genes',y=1.08)
# plot.savefig('../output_images/Testing/common-go-terms-of-'+ gene +'-based-on-their-type.png',dpi=300,format='png',transparent=True)

#%% testing pairplot


gene = 'BEM1'
tmp = common_go_data[common_go_data['query']==gene]
sns.set(style="ticks", color_codes=True)
# plot=sns.pairplot(tmp,vars=['fraction-of-common-go','common_interactors'],hue='score',hue_order=['NG','SL','PG'],
                  # diag_kind="hist", diag_kws = {'bins':int(np.ceil(np.sqrt(tmp.shape[0])))},corner=True)
plot=sns.pairplot(tmp,vars=['fraction-of-common-go','common_interactors'],hue='score',hue_order=['NG','SL','PG'],
                  diag_kind="hist", diag_kws = dict(bins=int(np.ceil(np.sqrt(tmp.shape[0]))), fill=False),corner=True)
# plt.title(query[0])
# sns.histplot
# plt.fig.suptitle(gene) #fixed title location

# plot.savefig('../output_images/Testing/EXAMPLE_common-go-terms-of-'+ gene +'-based-on-their-type.png',dpi=300,format='png',transparent=True)




#%% plot manually

gene = 'BEM1'
for gene in query:

    # pltData = common_go[common_go['query']==gene]
    geneCommonGoData = common_go_data[common_go_data['query']==gene]
    geneCommonInteractorData = common_partners_data[common_partners_data['query']==gene]
    pltData = common_go_data[common_go_data['query']==gene]

binList = int(np.ceil(np.sqrt(tmp.shape[0])))

ax1 = plt.subplot(2,2,1)
plotCommonHist(ax1,geneCommonGoData['fraction-of-common-go'],binList)
ax1.set_xlabel('% of GO-terms in common')
ax1.set_title('Common Go Terms')

ax2 = plt.subplot(2,2,2)
ax2.set_xlabel('% of interactors in common')
ax2.set_title('Common interactors') 
plotCommonHist(ax2,geneCommonInteractorData['fraction-of-common-partners'],binList)

ax3 = plt.subplot(2,1,2)
ax3.set_title('Common interactors') 
plotCommonScatter(ax3,geneCommonGoData['fraction-of-common-go'],geneCommonInteractorData['fraction-of-common-partners'])
# plt.xlim[0,100]
plt.subplots_adjust(wspace=0.3,hspace=0.7)

# plot.savefig('../output_images/Testing/EXAMPLE_alt_common-go-terms-of-'+ gene +'-based-on-their-type.png',dpi=300,format='png',transparent=True)

#%%

sns.jointplot(geneCommonGoData['fraction-of-common-go'],geneCommonInteractorData['fraction-of-common-partners'],height=5, ratio=2, marginal_kws=dict(bins=int(np.ceil(np.sqrt(tmp.shape[0])))));

#%% use jointplot for only 2 features

import seaborn as sns
# tmp = common_go_data[common_go_data['query']==gene]
tmp = common_go_data

# tmp["common_interactors"] = tmp["common_interactors"].astype(float)

tmp.astype({'fraction-of-common-go': 'float'}).dtypes
tmp.astype({'common_interactors': 'float'}).dtypes

# penguins = sns.load_dataset("penguins")
sns.jointplot(data=tmp, x="fraction-of-common-go", y="common_interactors", hue="score", height=5, ratio=2, marginal_ticks=(True))

    
#%% Saving the figure

#  plot.savefig('../output_images/common-go-terms-of-'+ genes +'-based-on-their-type.png',dpi=300,format='png',transparent=True)
