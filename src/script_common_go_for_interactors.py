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
from python_modules.module_common_measures import common_go
from python_modules.module_common_measures import common_partners
from python_modules.module_common_measures import plotCommonHist
from python_modules.module_common_measures import plotCommonScatter
#%% getting the data

data_go=pd.read_excel('../datasets/slim-goterms-filtered-data.xlsx')
# data_go=pd.read_excel('../datasets/testDataGo.xlsx')

data_go.columns=['Gene','gene-id','go-aspect','go-term','go-id','feature-type']

data=pd.read_excel('../datasets/data-BioGrid-Yeast-doubled.xlsx')
# data=pd.read_excel('../datasets/testDataInteractions.xlsx')
#%% query
query = pd.read_excel('../datasets/genesCellPolarity_SGD_amigo.xlsx')
query = np.unique(query).tolist()
# query=['BEM1','BEM2']
# query=['geneA','GeneB']

#%% Calling the function common_partners

common_partners_data=common_partners(query,data)

common_go_data=common_go(data_go=data_go,data_common_partners=common_partners_data)


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
 
#%% Working triangle pairplot (hist)
gene = 'BEM1'
# for gene in query:
tmp = common_go_data[common_go_data['query']==gene]
sns.set(style="ticks", color_codes=True)
# plot=sns.pairplot(tmp,vars=['fraction-of-common-go','common_interactors'],hue='score',hue_order=['NG','SL','PG'], kind = "hist", diag_kind="hist",corner=True)
plot=sns.pairplot(tmp,vars=['fraction-of-common-go','common_interactors'],hue='score',hue_order=['NG','SL','PG'], corner=True)
# plt.title(query[0])
plot.fig.suptitle(gene)
# plot.savefig('../output_images/Testing/common-go-terms-of-'+ gene +'-based-on-their-type.png',dpi=300,format='png',transparent=True)


#%% 

gene = 'BEM1'
# for gene in query:
tmp = common_go_data[common_go_data['query']==gene]
sns.set(style="ticks", color_codes=True)
plot=sns.pairplot(tmp,vars=['fraction-of-common-go','common_interactors'],hue='score',hue_order=['NG','SL','PG'], kind = "hist", diag_kind="hist",corner=True)
# plt.title(query[0])
plt.title(gene)
plot.savefig('../output_images/Testing/common-go-terms-of-'+ gene +'-based-on-their-type.png',dpi=300,format='png',transparent=True)

#%% viz
gene = 'BEM1'
tmp = common_go_data[common_go_data['query']==gene]
sns.set(style="ticks", color_codes=True)
plot=sns.pairplot(tmp,vars=['fraction-of-common-go','common_interactors'],hue='score',hue_order=['NG','SL','PG'],
                  diag_kind="hist", diag_kws = {'bins':int(np.ceil(np.sqrt(tmp.shape[0])))},corner=True)
# plt.title(query[0])
plot.fig.suptitle(gene) #fixed title location
plot.savefig('../output_images/Testing/EXAMPLE_common-go-terms-of-'+ gene +'-based-on-their-type.png',dpi=300,format='png',transparent=True)

#%% plot manually

for gene in query:
# gene = 'BEM1'
    # pltData = common_go[common_go['query']==gene]
    geneCommonGoData = common_go_data[common_go_data['query']==gene]
    geneCommonInteractorData = common_partners_data[common_partners_data['query']==gene]
    pltData = common_go_data[common_go_data['query']==gene]
# fig = plt.subplots(2,2)
# fig.suptitle(gene, fontsize=14)
# binList = np.linspace(0,100,11)
# binList = int(np.ceil(np.sqrt(tmp.shape[0])))
# diag_kws = {'bins':int(np.ceil(np.sqrt(tmp.shape[0])))}
# ax1 = plt.subplot(2,2,1)

# plotCommonHist(ax1,geneCommonGoData['fraction-of-common-go'],binList)
# ax1.set_xlabel('% of GO-terms in common')
# ax1.set_title('Common Go Terms')

# ax2 = plt.subplot(2,2,2)

# plotCommonHist(ax2,geneCommonInteractorData['fraction-of-common-partners'],binList)
# ax2.set_xlabel('% of interactors in common')
# ax2.set_title('Common interactors') 
# plotCommonHist(ax2,geneCommonInteractorData['fraction-of-common-partners'],binList)

# ax3 = plt.subplot(2,1,2)
# ax3.set_title('Common interactors') 
# plotCommonScatter(ax3,geneCommonGoData['fraction-of-common-go'],geneCommonInteractorData['fraction-of-common-partners'])
# # plt.xlim[0,100]
# plt.subplots_adjust(wspace=0.3,hspace=0.7)

# plot.savefig('../output_images/Testing/EXAMPLE_alt_common-go-terms-of-'+ gene +'-based-on-their-type.png',dpi=300,format='png',transparent=True)

#%%

sns.jointplot(geneCommonGoData['fraction-of-common-go'],geneCommonInteractorData['fraction-of-common-partners'],height=5, ratio=2, marginal_kws=dict(bins=int(np.ceil(np.sqrt(tmp.shape[0])))));


    
#%% Saving the figure

#  plot.savefig('../output_images/common-go-terms-of-'+ genes +'-based-on-their-type.png',dpi=300,format='png',transparent=True)
