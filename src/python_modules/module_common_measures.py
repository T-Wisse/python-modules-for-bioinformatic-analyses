# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 12:30:12 2020

@author: linigodelacruz
"""
import numpy as np
from collections import defaultdict 
import pandas as pd
import matplotlib.pyplot as plt


## Function to find common (how many) interaction partners of one gene with its interactors


def common_partners(query,data):
    
    """
    A function that finds the common interactors partners between the query genes provided from the data provided.
    query: dataframe with the target genes to look at
    data:dataframe with the info to take from
    """
        
    d2 = defaultdict(dict)
   
    # giant for loop
    
    for genes_names in query:
        #filtering the table just for the value of the query
        q1 = data[data['gene-query-name']==genes_names]
        q1_interact=q1['gene-target-name'].unique()
        
            
        for query2 in q1_interact:

        
            q2=data[data['gene-query-name']==query2] #these are get_query(q1[i])

            q2_interact=q2['gene-target-name'].unique()
        

            d = defaultdict(int)
            common = []
            for genes_names2  in q2_interact:
                if genes_names2 in q1_interact: # if a gene interactor of the query1 is in interactors of query 2 
                    common.append(genes_names2)
                    # tmp=np.array(q1[q1['gene-target-name']==genes_names2]['interaction-type'].tolist())
                    # type_interaction.append(tmp.ravel().ravel().ravel())
                    d[genes_names2] += 1
                   
            d2[query2]['query']=genes_names


            tmp=q1[q1['gene-target-name']==query2]['interaction-type'].tolist()
            #tmp1.append(tmp)
            tmp=np.unique(np.array(tmp))
           
            if "Synthetic Lethality" in tmp:
                d2[query2]["Type"] = "Synthetic Lethality"
                        
            else:
                d2[query2]["Type"] = tmp[0]

            d2[query2]["common"] = common
            d2[query2]["names of genes"]=query2
            d2[query2]["n_common"] = len(common)
            d2[query2]["number of partners of pairA"] = len(q1_interact)
            d2[query2]["number of partners of pairB"] = len(q2_interact)

            
            
            if len(q1)==0 :
                d2[query2]["fraction-of-common-partners"] = 0
            else:
                d2[query2]["fraction-of-common-partners"] = len(d)/len(q1_interact) *100     
       

        df=pd.DataFrame(d2).T
        #df.set_index=query2
        if len(df)==0:
            df_sorted=[]
        else:
            df_sorted=df.sort_values(by=["fraction-of-common-partners"])
            df_sorted=df_sorted[::-1]

        

    return df_sorted


def common_go_paralogs(data_go,data_paralogs):

 """
 a function that finds the common go terms for paralogs genes
 paralogs: paralogs are homologous genes that have evolved by duplication and code for protein with similar, but not identical functions.
 input: data_go= dataframe where all to go terms are 
        data_paralogs= dataframe where all the paralogs are
 output: dataframe with the fraction of common go for each paralog pair
 """
 
 query=np.unique(np.array(data_paralogs['query']))
 # big for loop for each gene analyzed in common partners
 d2=defaultdict(dict)

    
 for genes,i in zip(data_paralogs['paralogue-name'],query):
    d2[genes]['query']=i
    d2[genes]['names of paralogue']=genes

    tmp=data_go[data_go['Gene']==i]['go-term'].tolist()
    tmp=np.unique(tmp).tolist()

    

    tmp2=data_go[data_go['Gene']==genes]['go-term'].tolist()
    tmp2=np.unique(tmp2).tolist()

                


    d2[genes]['common-go-terms']=np.intersect1d(tmp,tmp2)
    if len(tmp)==0:
        d2[genes]['fraction-of-common-go']=0
    else:
        d2[genes]['fraction-of-common-go']=len(np.intersect1d(tmp,tmp2))/len(tmp) *100

    
    
 df=pd.DataFrame(d2).T
 df_sorted=df.sort_values(by='fraction-of-common-go',ascending=False)
      
    
    

 return df_sorted


def common_go(data_go,data_common_partners):
    """"
    function that computes the common go terms or interactors genes
    input: data_go= dataframe with all go terms per gene
    data_common_partners= dataframe output of the function common_partners
    
    """
    d3=defaultdict(dict)
    query=np.unique(np.array(data_common_partners['query'])) # Finds all query genes listed with a common partner
    # big for loop for each gene analyzed in common partners
    for i in np.arange(0,len(query)):
        partners=data_common_partners[data_common_partners['query']==query[i]]['names of genes']

        d2=defaultdict(dict)
        
        # print(partners)
        # print(query)
        for genes in partners:
            d2[genes]['query']=query[i]
            d2[genes]['names of genes']=genes

            tmp=data_go[data_go['Gene']==query[i]]['go-term'].tolist()
            tmp=np.unique(tmp).tolist()

            

            tmp2=data_go[data_go['Gene']==genes]['go-term'].tolist()
            tmp2=np.unique(tmp2).tolist()

                        

       
            d2[genes]['common-go-terms']=np.intersect1d(tmp,tmp2)
            if len(tmp)==0:
                d2[genes]['fraction-of-common-go']=0
            else:
                d2[genes]['fraction-of-common-go']=len(np.intersect1d(tmp,tmp2))/len(tmp) *100
                
      
        d3.update(d2)
        
    df=pd.DataFrame(d3).T
    df_sorted=df.sort_values(by='fraction-of-common-go',ascending=False)

    return df_sorted

def plotCommonHist(ax,data,binList):

    ax.set_ylabel('Count') #normalize?
   
    out =  ax.hist(data,binList)

    
    return out

def plotCommonScatter(ax,data1,data2):

    ax.set_xlabel('% of GO-terms in common')
    ax.set_ylabel('% of interactors in common')
       
    out =  ax.scatter(data1,data2)

    
    return out

def commonPhenotype(dataPhenotype,dataCommonPartners):
    
    d3=defaultdict(dict)
    query=np.unique(np.array(dataCommonPartners['query'])) # Finds all query genes listed with a common partner
    # big for loop for each gene analyzed in common partners
    for i in np.arange(0,len(query)):
        partners=dataCommonPartners[dataCommonPartners['query']==query[i]]['names of genes']

        d2=defaultdict(dict)

        for genes in partners:
            d2[genes]['query']=query[i]
            d2[genes]['names of genes']=genes

            tmp=dataPhenotype[dataPhenotype['gene-query-name']==query[i]]['interaction-type'].tolist()
            tmp=np.unique(tmp).tolist()

            tmp2=dataPhenotype[dataPhenotype['gene-query-name']==genes]['interaction-type'].tolist()
            tmp2=np.unique(tmp2).tolist()

            d2[genes]['common-phen-interactor']=np.intersect1d(tmp,tmp2)
            if len(tmp)==0:
                d2[genes]['fraction-of-common-phen']=0
            else:
                d2[genes]['fraction-of-common-phen']=len(np.intersect1d(tmp,tmp2))/len(tmp) *100
            
  
        d3.update(d2)
        
    df=pd.DataFrame(d3).T
    df_sorted=df.sort_values(by='fraction-of-common-phen',ascending=False)

    return df_sorted
    
#%% modulesT

def common_member(a, b): 
    a_set = set(a) 
    b_set = set(b) 
    if len(a_set.intersection(b_set)) > 0: 
        return(a_set.intersection(b_set))  
    return()

def common_partners_T(query,interactionData): #interaction data must have a 'gene query name' and a 'gene target name' column
    '''
    Finds common interactors for all genes in the query based on supplied interaction data
    query: list of genes
    ineractionData: dataframe of genetic interactions, contains gene query, gene target and type of interaction
    '''
    # d2 = defaultdict(dict)
    commonInteractors = defaultdict(dict)
    
    for queryGene in query: #For each query gene
    
        query1 = interactionData[interactionData['gene-query-name']==queryGene] #Reduce list to only query gene
        interactorList = query1['gene-target-name'].unique() #List genes interacting with query1 gene 
        
        for interactorGene in interactorList:
            
            query2 = interactionData[interactionData['gene-query-name']==interactorGene]
            interactors = query2['gene-target-name'].unique() #List genes interacting with query2 gene
            tempCommonInt = common_member(interactorList, interactors)
            commonInteractors[queryGene] = tempCommonInt
        
            # for interactorName  in interactors:
            #     if interactorName in interactorList: # if a gene interactor of the query1 is in interactors of query 2 
            #         commonInteractors.append(interactorName) #add gene to list
                    
    return commonInteractors #note, so far will only work properly for 1 gene in query...