
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 14:30:16 2023

@author: Erick Santiago

    This script calls an Rscript to compute phylogenetic distances and to fetch the closest parental for each hybrid.

"""
import os,csv
import subprocess
import re
import pandas as pd
import argparse


ag = argparse.ArgumentParser()
ag.add_argument("-f", "--hybrids_file", required=True, help = "List of hybrid IDs (file with one column containing one ID per row)")
ag.add_argument("-s", "--species", required=True, help = "SACE or SAPA. Specify the species you are working with, this is used to select a working directory")
args=vars(ag.parse_args())

hyb_file=args["hybrids_file"]
species=args["species"]

wd="/mnt/Timina/lmorales/Public/ymez/data/phylogeny/esantiago_HybTrees/" + species + "/"
scripts_dir="/mnt/Timina/lmorales/Public/ymez/bin/scripts/06_genotyping/"

output=open(wd + "HybridParentals_NatHybs_" + species +  "_noN.csv",'w')
output_writer=csv.writer(output)
header=["Hybrid_ID","Closest_"+ species +"_ID","Phylogenetic_distance_(avg.sust.p/site)","Total_SNPs_used_in_phylogeny","SharedSNPs_HybxParental"]
output_writer.writerow(header)

with open(hyb_file) as hybrids:
    for hybrid in hybrids.readlines():
        hyb_dir= wd + hybrid.strip() + "/"
        list2write=[]
        list2write.append(hybrid.strip())
        tree= "RAxML_bipartitionsBranchLabels.Matrix_SNPs_" + species  +"_from_CONC_gt_onlySNPs_filtered_missing_10_plus_" + hybrid.strip() + "_recode_min4_AllHybSites_noN.tree"
        process=["Rscript", scripts_dir + "getHybDists.R", hyb_dir + tree, hybrid.strip()]
        p=subprocess.run(process,stdout=subprocess.PIPE)
        stdout=str(p.stdout).split()
        closest_parental=stdout[0].split("'")[1]
        list2write.append(closest_parental)
        list2write.append(stdout[1].split("'")[0])
        phylip="Matrix_SNPs_" + species + "_from_CONC_gt_onlySNPs_filtered_missing_10_plus_" + hybrid.strip() + ".recode.min4_AllHybSites_noN.phy"
        with open(hyb_dir +phylip) as phylip_file:
            first_line=phylip_file.readline().split()[1]
        list2write.append(first_line)
        py_process=["python3.8", scripts_dir + "count_snps.py", hybrid.strip(), hyb_dir + phylip, closest_parental,species]
        py_p=subprocess.run(py_process,stdout=subprocess.PIPE)
        py_stdout=str(py_p.stdout)
        SNP_count=py_stdout.split("'")[1].strip()
        list2write.append(SNP_count[:-2])
        output_writer.writerow(list2write)
        print("################# Done: " + hybrid.strip() + " #################")
hybrids.close()
output.close()
