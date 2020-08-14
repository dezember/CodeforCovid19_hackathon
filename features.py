#Count features
#Author: Anurag Kanase
# Version: 1
# Github: dezember
# Date: August 12, 2020
# Project: Precision Drug Prediction

import pandas as pd
base = "/home/anurag/Seq/Counts/"
H1="SRR11804725"
H2="SRR11804726"
H3="SRR11804727"
H4="SRR11804728"

D1="SRR11804721"
D2="SRR11804722"
D3="SRR11804723"
D4="SRR11804724"

df1 = pd.read_csv(base+H1+".gene.txt", delimiter="\t", skiprows=1) 
df2 = pd.read_csv(base+H2+".gene.txt", delimiter="\t", skiprows=1)  
df3 = pd.read_csv(base+H3+".gene.txt", delimiter="\t", skiprows=1)
df4 = pd.read_csv(base+H4+".gene.txt", delimiter="\t", skiprows=1)
df5 = pd.read_csv(base+D1+".gene.txt", delimiter="\t", skiprows=1) 
df6 = pd.read_csv(base+D2+".gene.txt", delimiter="\t", skiprows=1)  
df7 = pd.read_csv(base+D3+".gene.txt", delimiter="\t", skiprows=1)
df8 = pd.read_csv(base+D4+".gene.txt", delimiter="\t", skiprows=1)

df = pd.concat([df1[df1.columns[0]],df1[df1.columns[6]],df2[df2.columns[6]],df3[df3.columns[6]],df4[df4.columns[6]],df5[df5.columns[6]],df6[df6.columns[6]],df7[df7.columns[6]],df8[df8.columns[6]]],axis=1,keys=['Gene','H1','H2','H3','H4','I1','I2','I3','I4'])

df.to_csv(base+"Features_1.csv",index=False)