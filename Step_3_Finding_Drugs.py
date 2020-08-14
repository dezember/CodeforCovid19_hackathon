"""Find drug candidates for the thresholded genes"""
# Author: Anurag Kanase
# Version: 1
# Github: dezember
# Date: August 12, 2020
# Project: Precision Drug Prediction

import pandas as pd
import numpy as np
from pyensembl import EnsemblRelease as pe
from urllib import request
from bs4 import BeautifulSoup
import json
from pyensembl import ensembl_grch38, cached_release

# find drugs 
def drugsfinder(geneName):
    """Find druggable candidate for the given gene name """
    drugs =[]
    max_drug = 4; # Maximum number of druggable candidates
    url = "http://dgidb.org/api/v2/interactions.json?genes="+geneName
    page = request.urlopen(url).read()
    soup = BeautifulSoup(page,'html.parser')
    drugpage=json.loads(soup.text)

    for druglist in drugpage['matchedTerms']:
        for drugnames in druglist['interactions'][:max_drug]: #remove max_drug to consider all the druggable candidadates
            drugs.append(drugnames['drugName'])
    return (drugs)


# Read the exported Differential Gene Expression find to find the Drug Candidates
df = pd.read_csv("diff_exp_results.csv")

#Stripping Ensembl ID name 
df['Gene']=df['Gene'].str[5:]

df = df[df.log2FoldChange >3 ]

#Store the Newly found Drugs for the thresholded Genes in df_DG
df_DG = pd.DataFrame(columns = ('Gene','Drug'))

#Import annotations. You will need to install: 
#pyensembl install --release 100 --species homo_sapiens
#if the file doesn't work
geneDB = cached_release(100, "human")

#Reading finding drugs available for the genes
print("finding drugs available for the genes...")
genes_list = []
drugs_list = []
for i in range(len(df['Gene'])):
    try:
        genename = geneDB.gene_name_of_gene_id(df['Gene'][i])
        drugs = drugsfinder(genename)
        genes_list = genes_list+[genename]*len(drugs)
        drugs_list = drugs_list+drugs
    except:
        pass

df_DG['Gene']=genes_list
df_DG['Drug']=drugs_list
df_DG.to_csv("Gene_Drug.csv")
print("Data saved to file: Gene_Drug.csv")