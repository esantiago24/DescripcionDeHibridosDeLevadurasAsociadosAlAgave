#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 18:30:01 2020

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
##########################################
ag = argparse.ArgumentParser()
ag.add_argument("-r", "--references", required = True, help = "file with the ids of the genome references to prepare") 
ag.add_argument("-s", "--samples", required = True, help = "csv file with the data of the samples to work") 
ag.add_argument("--hyb", required = False, help = "flag to indicate if strain is classified as hybrid. Use only when remapping SAPA reads to SAPA reference.")

args = vars(ag.parse_args())
refs = args["references"]
samples = args["samples"]
hyb = args["hyb"]

datos = getdata(samples,"SRA")
ids = get_refs(refs)

np="10"
ram="8G"


#sge_dir = "/mnt/C0E023C8E023C40E/Projects/MezcalYeast/sge_tests"
sge_dir = "/mnt/Timina/lmorales/Public/ymez/bin/SGE/03_mapping"
ref_dir = "/mnt/Timina/lmorales/Public/ymez/data/ref"
out_dir="/mnt/Timina/lmorales/Public/ymez/data/bam"
if hyb is None:
    fq_dir="/mnt/Timina/lmorales/Public/ymez/data/fastq/clean"
elif hyb == "1":
    fq_dir="/mnt/Timina/lmorales/Public/ymez/tmp/03_mapping"
else:
    fq_dir="/mnt/Timina/lmorales/Public/ymez/data/fastq/clean"

fig_dir="/mnt/Timina/lmorales/Public/ymez/data/figures/03_mapping"
stat_dir="/mnt/Timina/lmorales/Public/ymez/stats/mapping"
tmp="/mnt/Timina/lmorales/Public/ymez/tmp/03_mapping"
username = getpass.getuser()
#%%
sh_out = sge_dir + "/" + username + "_SH_map2ref.sh"
#path, sample, sge, app, user, np, ram
sgs = []
for rf in range(len(ids["short"])):
    for sample in list(datos.keys()):
    
        sgs.append(sample + "_vs_" + ids["short"][rf])
        file = sge_dir + "/" + sample + "_vs_" + ids["short"][rf] + ".sge"
        with open(file,'w') as sge_file:
             header(sge_dir, sample + "_" + ids["short"][rf], sge_file, "map", username, np, ram)
             print('module load bwa/0.7.4 htslib/1.2.1 gcc/5.1.0 samtools/1.9 picard/2.6.0 r/3.6.1 \n###', file =sge_file)
             print('start=$(date +%s.%N)', file = sge_file)
             print('ref_dir=%s/%s_v1/fasta' % (ref_dir, ids["large"][rf]), file = sge_file)
             print('out_dir=%s' % out_dir, file = sge_file)
             print('fq_dir=%s' % fq_dir, file = sge_file)
             print('fig_dir=%s' % fig_dir, file = sge_file)
             print('stat_dir=%s' % stat_dir, file = sge_file)   
             print('tmp=%s' % tmp, file = sge_file)
             if hyb == "1":
                 print('bwa mem -M -t10 ${ref_dir}/%s_v1_allChr.fasta ${fq_dir}/%s_CONC_SAPA_R1.fastq.gz ${fq_dir}/%s_CONC_SAPA_R2.fastq.gz | samtools view -hbS - | samtools sort -o ${out_dir}/%s_%s.bam - ' % (ids["large"][rf],  sample, sample, sample, ids["short"][rf]), file =sge_file)
             else:
                 print('bwa mem -M -t10 ${ref_dir}/%s_v1_allChr.fasta ${fq_dir}/%s_R1_clean.fastq.gz ${fq_dir}/%s_R2_clean.fastq.gz | samtools view -hbS - | samtools sort -o ${out_dir}/%s_%s.bam - ' % (ids["large"][rf],  sample, sample, sample, ids["short"][rf]), file =sge_file)
             print('picard ValidateSamFile I=${out_dir}/%s_%s.bam MODE=SUMMARY O=${tmp}/%s_%s_bam_0_status.txt' % (sample, ids["short"][rf], sample, ids["short"][rf]), file=sge_file)
             print('picard CollectGcBiasMetrics R=${ref_dir}/%s_v1_allChr.fasta I=${out_dir}/%s_%s.bam O=${stat_dir}/%s_%s_GCBias.txt CHART=${fig_dir}/%s_%s_GCBias.pdf ASSUME_SORTED=true SUMMARY_OUTPUT=${stat_dir}/%s_%s_summary_metrics.txt VALIDATION_STRINGENCY=LENIENT' % (ids["large"][rf], sample, ids["short"][rf], sample, ids["short"][rf], sample, ids["short"][rf], sample, ids["short"][rf] ), file =sge_file)
             print('picard MeanQualityByCycle R=${ref_dir}/%s_v1_allChr.fasta I=${out_dir}/%s_%s.bam O=${stat_dir}/%s_%s_Qcycle.txt CHART=${fig_dir}/%s_%s_Qcycle.pdf VALIDATION_STRINGENCY=LENIENT' % (ids["large"][rf], sample, ids["short"][rf], sample, ids["short"][rf], sample, ids["short"][rf]), file =sge_file)
             print('picard QualityScoreDistribution R=${ref_dir}/%s_v1_allChr.fasta I=${out_dir}/%s_%s.bam O=${stat_dir}/%s_%s_Qdist.txt CHART=${fig_dir}/%s_%s_Qdist.pdf VALIDATION_STRINGENCY=LENIENT' % (ids["large"][rf], sample, ids["short"][rf], sample, ids["short"][rf], sample, ids["short"][rf] ), file = sge_file)
             print('picard MarkDuplicates INPUT= ${out_dir}/%s_%s.bam OUTPUT= ${out_dir}/%s_%s.rmdup.bam METRICS_FILE=${stat_dir}/%s_%s_duplicateMatrix VALIDATION_STRINGENCY=LENIENT' % (sample, ids["short"][rf], sample, ids["short"][rf], sample, ids["short"][rf]), file =sge_file)
             print('picard ValidateSamFile I=${out_dir}/%s_%s.rmdup.bam MODE=SUMMARY O=${tmp}/%s_%s_bam_1_status.txt' % (sample, ids["short"][rf], sample, ids["short"][rf]), file=sge_file)
             print('picard AddOrReplaceReadGroups I=${out_dir}/%s_%s.rmdup.bam O=${out_dir}/%s_%s.rmdup.addgp.bam LB=%s PL=illumina PU=%s SM=%s VALIDATION_STRINGENCY=LENIENT' % (sample, ids["short"][rf], sample, ids["short"][rf],sample, sample, sample), file =sge_file)
             print('picard ValidateSamFile I=${out_dir}/%s_%s.rmdup.addgp.bam MODE=SUMMARY O=${tmp}/%s_%s_bam_2_status.txt' % (sample, ids["short"][rf], sample, ids["short"][rf]), file=sge_file)
             print('samtools index ${out_dir}/%s_%s.rmdup.addgp.bam' % (sample, ids["short"][rf]), file = sge_file)        
#             print('java -Xmx8g -jar $(which GenomeAnalysisTK).jar -T RealignerTargetCreator -R ${ref_dir}/%s_v1_allChr.fasta -I ${out_dir}/%s_%s.rmdup.addgp.bam -o ${out_dir}/%s_%s.rmdup.addgp.intervals' % (ids["large"][rf], sample, ids["short"][rf], sample, ids["short"][rf]), file =sge_file)
#             print('java -Xmx8g -jar $(which GenomeAnalysisTK).jar -T IndelRealigner -R ${ref_dir}/%s_v1_allChr.fasta -I ${out_dir}/%s_%s.rmdup.addgp.bam -o ${out_dir}/%s_%s.realigned.bam --maxReadsForConsensuses 500 --maxReadsForRealignment 40000 --maxConsensuses 60 --maxPositionalMoveAllowed 400 --entropyThreshold 0.2 -targetIntervals ${out_dir}/%s_%s.rmdup.addgp.intervals' % (ids["large"][rf], sample, ids["short"][rf], sample, ids["short"][rf] ,sample, ids["short"][rf] ), file =sge_file)
             print('if [[ -s ${out_dir}/%s_%s.rmdup.addgp.bam ]]; then' % (sample, ids["short"][rf]), file = sge_file)
            # print('rm -f ${out_dir}/%s_%s.bam' % (sample, ids["short"][rf]), file = sge_file)
             print('rm -f ${out_dir}/%s_%s.rmdup.bam' % (sample, ids["short"][rf]), file = sge_file)
#             print('#rm ${out_dir}/%s_%s.rmdup.addgp.bam' % (sample, ids["short"][rf]), file = sge_file)
             print('fi', file = sge_file)
             print('duration=$(echo "$(date +%s.%N) - $start" | bc)', file =sge_file)
             print('execution_time=`printf "%.2f seconds" $duration`', file = sge_file)
             print('echo "Script Execution Time: $execution_time"', file =sge_file)

bashfile(sh_out, sgs)
