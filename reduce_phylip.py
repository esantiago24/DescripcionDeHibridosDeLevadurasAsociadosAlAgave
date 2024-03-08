#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thursday March 16, 2022
@author: Erick Santiago

This script receives a Phylip file and reduces the multiple sequence alignment to only the sites for which the desired sample (-s) has available information.

"""
import re
import pandas as pd
import csv
import argparse
from pathlib import Path
import numpy

ag=argparse.ArgumentParser()
ag.add_argument("-p", "--phylip", required=True, help="Absolute path to the Phylip file in strict or relaxed format." )
ag.add_argument("-s", "--sample_ID", required=True, help="ID of the sample of interest.")

args=vars(ag.parse_args())
sample=args["sample_ID"]
sample=sample + " " #Add a whitespace at the end to keep the phylip format requirements
input_file=args["phylip"]

cwd=Path(input_file)
cwd=cwd.parent


ldf=[]
i=0
#Parse the phylip file to store information in a pandas dataframe.
with open(input_file) as file:
    for line in file.readlines():
        fmtdrow=[]
        if (i==0):
            i+=1
            continue
        else:
            row=line.split()
            SNPS=list(row[1])
            ID=row[0] + " "
            fmtdrow.append(ID)
            fmtdrow.extend(SNPS)
            ldf.append(fmtdrow)
    df=pd.DataFrame(ldf)

#Transpose de dataframe to keep samples as columns and SNPs as rows
df_t=df.transpose()
df_t.columns=df_t.iloc[0]
df_t=df_t[1:]
df_t.reset_index(inplace=True,drop=True)

#Fetch the indexes of the rows of the column "Sample" to only keep those which not contain a N or - value
pos2keep=df_t.index[(df_t[sample] != "N")].tolist()
pos2rm=(df_t.index[df_t[sample] == ("-")].tolist())
pos2keep=[x for x in pos2keep if x not in pos2rm]

df_tr=df_t.iloc[pos2keep]
df_tr.reset_index(inplace=True,drop=True)

df_tr=df_tr.transpose()
print("Number of sites removed: ", (df.shape[1] - df_tr.shape[1]))
#print("Size ",df_t.shape[1]) #Number of columns: df.shape[1]

#Concatenate all the SNPS into a single string to print them as a single row on a new file
phylip_content=[]
first_line= str(df_tr.shape[0]) + " " + str(df_tr.shape[1]) + "\n"
phylip_content.append(first_line)
IDs=df[0].values.tolist()
df_tr["Concat"]= IDs + df_tr.values.sum(axis=1) + "\n"
phylip_content.extend(df_tr["Concat"].values.tolist())
phylip=open(input_file[:-4]+ "_AllHybSites_noN.phy","w")
phylip.writelines(phylip_content)
phylip.close()
