#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 20:31:50 2020

@author: luisfer
"""
#%%
import sys
import argparse
import getpass

sys.path.append("../generic_functions/")

from generic import getdata
from generic import header
from generic import bashfile
from generic import get_refs

#%%
ag = argparse.ArgumentParser()
ag.add_argument("-r", "--references", required = True, help = "file with the ids of the genome references to use")
ag.add_argument("-s", "--samples", required = True, help = "csv file with the samples to use")

args = vars(ag.parse_args())
refs = args["references"]
samples = args["samples"]

datos = getdata(samples, 'SRA')
ids = get_refs(refs)

np="4"
ram="6G"

#sge_dir = "/home/luisgarcia/projects/Mezcal_yeast/1011genomes/sge_test"
sge_dir = "/mnt/Timina/lmorales/Public/ymez/bin/SGE/05_vcalling"
bam_dir="/mnt/Timina/lmorales/Public/ymez/data/bam"
ref_dir = "/mnt/Timina/lmorales/Public/ymez/data/ref"
tmp_dir = "/mnt/Timina/lmorales/Public/ymez/tmp"
stats = "/mnt/Timina/lmorales/Public/ymez/stats/genotyping"
#%%
#path, sample, sge, app, user, np, ram
username = getpass.getuser()
sh_out = sge_dir + "/" + username + "_SH_vcalling_genotyping.sh"
sgs = []
for sample in list(datos.keys()):
    for rf in range(len(ids["short"])):
        sgs.append(sample + "_" + ids["short"][rf])
        file = sge_dir + "/" + sample + "_" + ids["short"][rf]  + ".sge"
        with open(file,'w') as sge_file:
             header(sge_dir, ids["short"][rf] + "_" + sample, sge_file, "VCalling", username, np, ram)
             print('module load gcc/5.1.0 bedops/2.4.20 bedtools/2.24 gatk/4.1.1.0 picard/2.6.0 samtools/1.9 \n###', file =sge_file)
             print('start=$(date +%s.%N)', file = sge_file)
             print('bam_dir=%s' % bam_dir, file = sge_file)
             print('tmp_dir=%s' % tmp_dir, file = sge_file)
             print('stats=%s' % stats, file = sge_file)
             print('ref_dir=%s/%s_v1/fasta' % (ref_dir, ids["large"][rf]), file = sge_file)
             print('FILE1=${bam_dir}/%s_%s.rmdup.addgp.bam' % (sample, ids["short"][rf]), file =sge_file)
             print('gatk --java-options "-Xmx16g" HaplotypeCaller -R ${ref_dir}/%s_v1_allChr.fasta -I ${FILE1} -O ${tmp_dir}/05_vcalling/$(basename ${FILE1} .rmdup.addgp.bam).g.vcf --annotation ClippingRankSumTest --annotation DepthPerSampleHC --annotation StrandBiasBySample --annotation AlleleFraction --annotation AS_FisherStrand --annotation ChromosomeCounts  --emit-ref-confidence GVCF' % ids["large"][rf], file=sge_file)
             print('gatk --java-options "-Xmx16g" GenotypeGVCFs -R ${ref_dir}/%s_v1_allChr.fasta -V ${tmp_dir}/05_vcalling/$(basename ${FILE1} .rmdup.addgp.bam).g.vcf  -O ${tmp_dir}/06_genotyping/$(basename ${FILE1} .rmdup.addgp.bam).gt.g.vcf --annotation ClippingRankSumTest --annotation DepthPerSampleHC --annotation StrandBiasBySample --annotation AlleleFraction --annotation AS_FisherStrand --annotation ChromosomeCounts' % ids["large"][rf], file= sge_file)
             print('gatk --java-options "-Xmx16g" SelectVariants -R ${ref_dir}/%s_v1_allChr.fasta -V ${tmp_dir}/06_genotyping/$(basename ${FILE1} .rmdup.addgp.bam).gt.g.vcf --select-type-to-include SNP -O ${tmp_dir}/06_genotyping/$(basename ${FILE1} .rmdup.addgp.bam).gt.SNP.g.vcf' % ids["large"][rf], file = sge_file)
             print('total_raw_snps_Saccharomyces=$(grep -P "SA\w+\t"  ${tmp_dir}/06_genotyping/$(basename ${FILE1} .rmdup.addgp.bam).gt.SNP.g.vcf | wc -l)', file = sge_file)
             print('total_raw_snps_KLMA=$(grep -P "KLMA\w+\t"  ${tmp_dir}/06_genotyping/$(basename ${FILE1} .rmdup.addgp.bam).gt.SNP.g.vcf | wc -l)', file = sge_file)
             print('total_raw_snps_PIKU=$(grep -P "PIKU\w+\t"  ${tmp_dir}/06_genotyping/$(basename ${FILE1} .rmdup.addgp.bam).gt.SNP.g.vcf | wc -l)', file = sge_file)
             print('total_raw_snps=$(($total_raw_snps_Saccharomyces + $total_raw_snps_KLMA + $total_raw_snps_PIKU))', file = sge_file)
             print('echo $(basename ${FILE1} .rmdup.addgp.bam)_total_SNPs = $total_raw_snps >> ${stats}/Isolates_total_SNPs_v1.txt', file = sge_file)
             print('duration=$(echo "$(date +%s.%N) - $start" | bc)', file =sge_file)
             print('execution_time=`printf "%.2f seconds" $duration`', file = sge_file)
             print('echo "Script Execution Time: $execution_time"', file =sge_file)

bashfile(sh_out, sgs)
