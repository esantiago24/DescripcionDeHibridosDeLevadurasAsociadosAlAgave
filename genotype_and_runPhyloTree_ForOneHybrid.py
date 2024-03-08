#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 14:30:16 2022

@author: luis
    modified by Erick Santiago

    This script creates a raxML tree for each hybrid in the samplesheet (-f). Each tree consists of a single hybrid with a set of parental strains (-s)

"""

#%%
import sys
import argparse
import getpass
import csv
import os
from pathlib import Path

sys.path.append("../generic_functions/")

from generic import getdata
from generic import header
from generic import bashfile
from generic import get_refs
from generic import unique
from textwrap import wrap
import uuid

#%%
ag = argparse.ArgumentParser()
ag.add_argument("-r", "--references", required = True, help = "file with the ids of the genome references to use")
ag.add_argument("-s", "--samples", required = True, help = "csv file with the samples to use")
ag.add_argument("-m", "--missing", required = True, help = "The proportion of missing data allowed in the final multisample vcf. Use only values ​​in the range 0 to 1. Multiple values ​​separated by commas are allowed (without white spaces).")
ag.add_argument("-f", "--hybrids_file", required=True, help = "List of hybrid IDs to generate an SGE for")
ag.add_argument("-d", "--delete_vcf", required=True, help = "Boolean: True or False. Delete VCF at the end of the script to save space")
ag.add_argument("-S", "--species", required=True, help = "SACE or SAPA. Specify the species you are working with, this is used to select a working directory")

args = vars(ag.parse_args())
refs = args["references"]
samples = args["samples"]
m = args["missing"]
m = m.split(",")
hyb_file=args["hybrids_file"]
del_vcf=args["delete_vcf"]
species=args["species"]

#hybrids=getdata(hyb_file,'SRA')
datos = getdata(samples, 'SRA')
ids = get_refs(refs)
hybs=getdata(hyb_file,'SRA')
#%%

np="8"
ram="2G"

#sge_dir = "/mnt/C0E023C8E023C40E/Projects/MezcalYeast/05_vcalling"
sge_dir = "/mnt/Timina/lmorales/Public/ymez/data/phylogeny/esantiago_HybTrees/" + species + "/SGE"
#ref_dir="/mnt/C0E023C8E023C40E/Projects/MezcalYeast/05_vcalling"
ref_dir = "/mnt/Timina/lmorales/Public/ymez/data/ref"
tmp_dir = "/mnt/Timina/lmorales/Public/ymez/tmp/05_vcalling"
missingVcfs_dir= "/mnt/Timina/lmorales/Public/ymez/data/phylogeny/esantiago_HybTrees/"+ species + "/SGE/MissingFiles"
#tmp_dir="/mnt/C0E023C8E023C40E/Projects/MezcalYeast/05_vcalling"


hybs=dict(hybs)
d = dict(datos)
#ids["short"]="CONC"  #This was modified in order to have the called vcf files have the format {ID}_CONC.SNP_onlychr_SAPA.g.vcf and not {ID}_SAPA.SNP_onlychr_SAPA.g.vcf
r = ids["short"]
#print("ids:" , ids)
#print("ids[short]:" , ids["short"])
#print("ids[large]:" ,ids["large"])


username = getpass.getuser()
sh_out = sge_dir + "/" + username + "_SH_getSNPmatrix_For_AllHybs.sh"
sgs = []

#%%
#Checking subgenome
for hybrid in hybs:
    for rf in range(len(ids["short"])):
        index= Path(ref_dir + "/"+ ids["large"][rf]+"_v1/fasta/"+ ids["large"][rf]+"_v1_allChr.fasta.fai")
        if index.is_file():
            with open(index, 'r') as csvfile:
                reader = csv.reader(csvfile, delimiter=' ')
                first_column = [ row[0].split("\t")[0].split("_")[0] for row in reader ]
                subgen = list(unique(first_column))

        else:
            print("Error: the file %s is missing, please prepare the corresponding index for the genome. See 03_prepare" % index)
            r.remove(ids["short"][rf])
            continue

        g = list(subgen)
        sufix= wrap(str(uuid.uuid4().hex), int(len(str(uuid.uuid4().hex))/len(g)) )
        zip_iterator = zip(g, sufix)
        sufix_dic = dict(zip_iterator)
        files_dict=dict()
        for sg in subgen:
            lista=list()
            print(sg)
            cuenta=0
            datos[hybrid]=None
            for sample in list(datos.keys()):
                my_file=Path("%s/%s_%s.SNP_onlychr_%s.g.vcf.gz" %(tmp_dir, sample, ids["short"][rf], sg) )
                if my_file.is_file():
                    lista.append(str(my_file))
                else:
                    cuenta += 1
                    log = missingVcfs_dir + "/missingVcfs_" + sg + "_from_" + ids["short"][rf] + "_" + sufix_dic[sg] + ".log"
                    with open(log,'a') as log_file:
                        print("%s/%s_%s.SNP_onlychr_%s.g.vcf.gz" %(tmp_dir, sample, ids["short"][rf], sg), file=log_file)

            if cuenta > 0:
                print("Warning: from " + ids["short"][rf] + " missing " + str(cuenta) + " vcf files for " + sg + ". A list of missing files can be found at %s" % log)
                g.remove(sg)
            else:
                files_dict.update({sg:lista})
                #print(g)

        print('')
        del datos[hybrid]
        if len(g) > 0:
            print('Working only with ' + " ".join(g) + " for reference " + ids["short"][rf])
        else:
            print('Reference %s was discarded because not all vcf files were found for the samples listed in the sampleSheet' % ids["large"][rf])
            continue

        stri = "--variant "
        stri += '% s'
        # lista = list(datos.keys())


        for sg in g:
            files_dict[sg] =  [stri % i for i in files_dict[sg]]
            sgs.append("get_" + sg + "_SNPmatrix_from_" + ids["short"][rf]+ "_plus_" + hybrid)
            file = sge_dir + "/get_" + sg + "_SNPmatrix_from_" + ids["short"][rf] + "_plus_" + hybrid + ".sge"
            with open(file,'w') as sge_file:
                out_dir="/mnt/Timina/lmorales/Public/ymez/data/phylogeny/esantiago_HybTrees/"+ species + "/" + hybrid
                header(sge_dir, ids["short"][rf], sge_file, hybrid + "_SNP_Matrix", username, np, ram)
                print("##Use node 10",file=sge_file)
                print("#$ -q all.q@compute-00-10",file=sge_file)
                print("#$ -q all.q@compute-00-11",file=sge_file)
                print("#$ -q all.q@compute-00-12",file=sge_file)
                print('module load gcc/5.1.0 python37/3.7.0 vcftools/0.1.14 gatk/4.1.1.0 samtools/1.10 htslib/1.9 bcftools/1.9 phyml/3.1 jmodeltest/2.1.10 raxml/8.2.12-avx2-pthreads \n###', file =sge_file)
                print('start=$(date +%s.%N)', file = sge_file)
                print('tmp_dir=%s' % tmp_dir, file = sge_file)
                print('ref_dir=%s' % ref_dir, file = sge_file)
                print('out_dir=%s' % out_dir, file = sge_file)

                print('cd /mnt/Timina/lmorales/Public/ymez/tmp/05_vcalling', file=sge_file)
                print("# Create output directory", file=sge_file)
                print("mkdir " + out_dir , file=sge_file)
                print('', file=sge_file)
                joinedVCF="/Matrix_SNPs_%s_from_%s_plus_%s.vcf" % (sg, ids["short"][rf],hybrid)
                if not os.path.exists(out_dir + joinedVCF):
                    print('# We merge the necessary gvcf files', file=sge_file)
                    print('gatk --java-options "-Xmx16g" CombineGVCFs -R ${ref_dir}/%s_v1_allChr.fasta %s -O ${out_dir}/Matrix_SNPs_%s_from_%s_plus_%s.vcf' % (ids["large"][rf]+"_v1/fasta/"+ ids["large"][rf], " ".join(files_dict[sg]), sg, ids["short"][rf],hybrid), file=sge_file)
                else:
                    print('# Merged GVCF file already exists: skipping gatk CombineGVCFs',file=sge_file)
                print('',file=sge_file)
                print('#  We genotype the multisample vcf file',file=sge_file)
                print('echo "############## Genotyping multisample VCF: ##############"',file=sge_file)
                print('gatk --java-options "-Xmx16g" GenotypeGVCFs -R ${ref_dir}/%s_v1_allChr.fasta -V ${out_dir}/Matrix_SNPs_%s_from_%s_plus_%s.vcf -O ${out_dir}/Matrix_SNPs_%s_from_%s_gt_plus_%s.vcf' %(ids["large"][rf]+"_v1/fasta/"+ ids["large"][rf], sg,ids["short"][rf],hybrid, sg, ids["short"][rf],hybrid), file=sge_file)
                print('',file=sge_file)
                print('#  We select only SNPs from the multisample gt.vcf',file=sge_file)
                print('echo "############## Keeping only SNPs: ##############"',file=sge_file)
                print('gatk --java-options "-Xmx16g" SelectVariants -R ${ref_dir}/%s_v1_allChr.fasta -V ${out_dir}/Matrix_SNPs_%s_from_%s_gt_plus_%s.vcf --select-type-to-include SNP -O ${out_dir}/Matrix_SNPs_%s_from_%s_gt_onlySNPs_plus_%s.vcf' % (ids["large"][rf]+"_v1/fasta/"+ ids["large"][rf], sg, ids["short"][rf], hybrid,sg, ids["short"][rf],hybrid), file=sge_file)
                print('',file=sge_file)
                print('# We filter low-quality SNPs',file=sge_file)
                print('echo "############## Filtering low-quality SNPs: ##############"',file=sge_file)
                print('gatk --java-options "-Xmx16g" VariantFiltration -R ${ref_dir}/%s_v1_allChr.fasta -V ${out_dir}/Matrix_SNPs_%s_from_%s_gt_onlySNPs_plus_%s.vcf -O ${out_dir}/Matrix_SNPs_%s_from_%s_gt_onlySNPs_filtered_plus_%s.vcf --filter-name "SNP_QD_filters" --filter-expression "QD<2.0" --filter-name "SNP_MQ_filters" --filter-expression "MQ<40.0" --filter-name "SNP_FS_filters" --filter-expression "FS>60.0" --filter-name "SNP_SOR_filters" --filter-expression "SOR>3.0" --filter-name "SNP_MQRankSum_filters" --filter-expression "MQRankSum<-12.5" --filter-name "SNP_ReadPosRankSum_filters" --filter-expression "ReadPosRankSum<-8.0"' % (ids["large"][rf]+"_v1/fasta/"+ ids["large"][rf], sg, ids["short"][rf], hybrid ,sg, ids["short"][rf],hybrid), file=sge_file)
                print('',file=sge_file)
                print('# We discard low-quality SNPs',file=sge_file)
                print('echo "############## Discarding low-quality SNPs: ##############"',file=sge_file)
                print('vcftools --vcf ${out_dir}/Matrix_SNPs_%s_from_%s_gt_onlySNPs_filtered_plus_%s.vcf --remove-filtered-all --recode --out ${out_dir}/Matrix_SNPs_%s_from_%s_gt_onlySNPs_filtered_plus_%s' % (sg, ids["short"][rf], hybrid,sg, ids["short"][rf],hybrid), file=sge_file)
                print('',file=sge_file)
                print('# Cleaning working directory',file=sge_file)
                print('mv ${out_dir}/Matrix_SNPs_%s_from_%s_gt_onlySNPs_filtered_plus_%s.recode.vcf ${out_dir}/Matrix_SNPs_%s_from_%s_gt_onlySNPs_filtered_plus_%s.vcf' % (sg, ids["short"][rf],hybrid,sg, ids["short"][rf],hybrid), file=sge_file)
                print('if [ -f ${out_dir}/Matrix_SNPs_%s_from_%s_gt_onlySNPs_filtered_plus_%s.vcf.gz ];then' % (sg, ids["short"][rf],hybrid), file=sge_file)
                print('echo "Warning: Removing previous file Matrix_SNPs_%s_from_%s_gt_onlySNPs_filtered_plus_%s.vcf.gz"'% (sg, ids["short"][rf],hybrid), file=sge_file)
                print('rm ${out_dir}/Matrix_SNPs_%s_from_%s_gt_onlySNPs_filtered_plus_%s.vcf.gz' % (sg, ids["short"][rf],hybrid), file=sge_file)
                print('fi',file=sge_file)
                print('bgzip ${out_dir}/Matrix_SNPs_%s_from_%s_gt_onlySNPs_filtered_plus_%s.vcf'% (sg, ids["short"][rf],hybrid), file=sge_file)
                print('bgzip ${out_dir}/Matrix_SNPs_%s_from_%s_plus_%s.vcf'% (sg, ids["short"][rf],hybrid), file=sge_file)
                print('',file=sge_file)
                print('#  We filter for missing data and conserved only biallelic SNPs',file=sge_file)
                for ms in m:
                    if float(ms) >=0 and float(ms) <=1:
                        print('# Missing data: %s' % ms, file=sge_file)
                        print('vcftools --gzvcf ${out_dir}/Matrix_SNPs_%s_from_%s_gt_onlySNPs_filtered_plus_%s.vcf.gz --max-missing %s --recode --out ${out_dir}/Matrix_SNPs_%s_from_%s_gt_onlySNPs_filtered_missing_%s_plus_%s' % (sg, ids["short"][rf], hybrid ,str(1- float(ms)), sg, ids["short"][rf], str(int(100*float(ms))),hybrid ), file=sge_file)
                        print('',file=sge_file)
                        print('bcftools view -M2 -v snps ${out_dir}/Matrix_SNPs_%s_from_%s_gt_onlySNPs_filtered_missing_%s_plus_%s.recode.vcf > ${out_dir}/Matrix_SNPs_%s_from_%s_gt_onlySNPs_filtered_missing_%s_biallelic_plus_%s.vcf' %(sg, ids["short"][rf], str(int(100*float(ms))), hybrid ,sg, ids["short"][rf], str(int(100*float(ms))),hybrid), file=sge_file)
                        print('',file=sge_file)
                        print('vcftools --vcf ${out_dir}/Matrix_SNPs_%s_from_%s_gt_onlySNPs_filtered_missing_%s_biallelic_plus_%s.vcf --freq --out ${out_dir}/Matrix_SNPs_%s_from_%s_gt_onlySNPs_filtered_missing_%s_biallelic_plus_%s' %(sg, ids["short"][rf], str(int(100*float(ms))), hybrid,sg, ids["short"][rf], str(int(100*float(ms))),hybrid), file=sge_file)
                        print('',file=sge_file)
                    else:
                        print("Only values ​​in the range 0 to 1 are used for missing data. %s is discarded." % ms)

                print('# Convert vcf to phylip format',file=sge_file)
                print('echo "############## Converting vcf to phylip format ##############"',file=sge_file)
                print('python3.7 vcf2phylip.py -i ${out_dir}/Matrix_SNPs_%s_from_%s_gt_onlySNPs_filtered_missing_%s_plus_%s.recode.vcf' %(sg, ids["short"][rf], str(int(100*float(ms))), hybrid), file=sge_file)
                print('',file=sge_file)
                print('# Reduce phylip to use only the sites where the sample %s has available information',file=sge_file)
                print('echo "############## Reducing phylip file ##############"',file=sge_file)
                print('python3.7 /mnt/Timina/lmorales/Public/ymez/bin/scripts/06_genotyping/reduce_phylip.py -p ${out_dir}/Matrix_SNPs_%s_from_%s_gt_onlySNPs_filtered_missing_%s_plus_%s.recode.min4.phy -s %s' %(sg, ids["short"][rf], str(int(100*float(ms))), hybrid, hybrid), file=sge_file)
                print('',file=sge_file)
                print('# Create phylogenetic tree',file=sge_file)
                print('echo "############## Running RaXML: ##############"',file=sge_file)
                if(species=="SAPA"):
                    print("raxmlHPC-PTHREADS-AVX2 -f a -x 12345 -p 12345 -N 100  -T 8 -m GTRGAMMA -s ${out_dir}/Matrix_SNPs_%s_from_%s_gt_onlySNPs_filtered_missing_%s_plus_%s.recode.min4_AllHybSites_noN.phy -n Matrix_SNPs_%s_from_%s_gt_onlySNPs_filtered_missing_%s_plus_%s_recode_min4_AllHybSites_noN.tree" %(sg, ids["short"][rf], str(int(100*float(ms))), hybrid,sg,ids["short"][rf], str(int(100*float(ms))),hybrid), file=sge_file)
                    print('',file=sge_file)
                else:
                    print("raxmlHPC-PTHREADS-AVX2 -f a -x 12345 -p 12345 -N 10 -T 8 -m GTRGAMMA -s ${out_dir}/Matrix_SNPs_%s_from_%s_gt_onlySNPs_filtered_missing_%s_plus_%s.recode.min4_AllHybSites_noN.phy -n Matrix_SNPs_%s_from_%s_gt_onlySNPs_filtered_missing_%s_plus_%s_recode_min4_AllHybSites_noN.tree" %(sg, ids["short"][rf], str(int(100*float(ms))), hybrid,sg,ids["short"][rf], str(int(100*float(ms))),hybrid), file=sge_file)
                    print('',file=sge_file)
                print('#Move RaXML output to working directory',file=sge_file)
                print('echo "############## Moving tree files to wd: ##############"',file=sge_file)
                print('mv RAxML_*_%s_*.tree ${out_dir}/' % (hybrid),file=sge_file)
                print('',file=sge_file)
                print('#Remove VCF to save space',file=sge_file)
                print('echo "############## Removing %s VCF file to save space: ##############"' % (hybrid),file=sge_file)
                print("if [ -f ${out_dir}/Matrix_SNPs_%s_from_%s_plus_%s.vcf]; then" % (sg, ids["short"][rf],hybrid),file=sge_file)
                print("rm ${out_dir}/Matrix_SNPs_%s_from_%s_plus_%s.vcf" % (sg, ids["short"][rf],hybrid),file=sge_file)
                print('fi',file=sge_file)
                print('',file=sge_file)
                print('duration=$(echo "$(date +%s.%N) - $start" | bc)', file =sge_file)
                print('execution_time=`printf "%.2f seconds" $duration`', file = sge_file)
                print('echo "Script Execution Time: $execution_time"', file =sge_file)
                print('########################################### END OF SAMPLE ###########################################',file= sge_file)

bashfile(sh_out,sgs )
#%%
