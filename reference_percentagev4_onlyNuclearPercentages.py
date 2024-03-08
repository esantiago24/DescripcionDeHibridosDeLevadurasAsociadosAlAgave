"""
Created on June 2023

@author: Iván Sedeño
    @modified by: Erick Santiago

    This script obtains the number of nuclear  reads mapped to each reference, computes the proportions of mapped reads and stores it in a csv.

"""

import os,csv
import subprocess
import re
import pandas as pd
import argparse

ag=argparse.ArgumentParser()
ag.add_argument("-s","--sample",required=True, help= "Absolute path to the csv with the samples to process. Example: /mnt/Timina/lmorales/Public/ymez/data/metadata/samplesheet.csv. File must have only an ID column with no header and one sample per line.")
args=vars(ag.parse_args())
samples=args["sample"]

path = "/mnt/Timina/lmorales/Public/ymez/stats/coverage/" #Path to stats/coverage/
script_dir = "/mnt/Timina/lmorales/Public/ymez/bin/scripts/04_coverage/" #Path of the invoked Rscript
id_file =samples #Sample sheet with the samples to be processed
OutRefPercentageFilePath="RefCovPercentageSheets/RefPercentageSheet_" + os.path.basename(samples)[:-4] + ".csv" #Name of the output csv that will be generated at the end of the program, contains mapped, unmapped and read % of each sample.

#Function to extract the IDs from the sample sheet. The file should only contain ID per line.
def get_ids(file):
    with open(file,newline="\n") as f:
        reader = csv.reader(f,delimiter=",")
        ids = list(reader)
    return ids

#Function to fetch stat files for each ID.
def get_files(path):
    ids = get_ids(id_file)
    files = {}
    for id in ids:
        files[id[0]] = path + id[0] + "/CONC/Q30" + id[0] + "_perReference.csv"
    return files

#Runs samtools idxstats to get number of mapped and unmapped reads for each reference, stores the information in a dataframe and returns it for further use.
def get_stats(id):
    file = f"/mnt/Timina/lmorales/Public/ymez/data/bam/{id}_CONC.rmdup.addgp.bam"
    p = "samtools idxstats " + file
    ps_total=subprocess.Popen(p,shell=True,stdout=subprocess.PIPE)
    ps_total.wait()
    total=str(ps_total.communicate()[0]).split("'")[1]
    total=total.split('\\n')
    totalSACE=total[0:16]
    totalSAPA=total[18:34]
    print("TOTALSACE: \n", totalSACE)
    print("TOTALSAPA: \n", totalSAPA)
    SAPA=re.sub(r'\\t','.',totalSAPA).split('.')
    SACE=re.sub(r'\\t','.',totalSACE).split('.')
    print("SAPA: \n", SAPA)
    print("SACE: \n", SACE)
    Unmpd=re.sub(r'\\t','.',total[2]).split('.')

    df=pd.DataFrame(columns=['Reference','Length','Mapped','Unmapped'])
    df.loc[len(df)]=SACE
    df.loc[len(df)]=SAPA
    df.loc[len(df)]=Unmpd
    return df

#Computes the total number of reads (Mapped + Unmapped Reads)
def get_total(df):
    Mapped=df['Mapped'].astype('int').sum()
    Unmapped=df['Unmapped'].astype('int').sum()

    total_reads=Mapped+Unmapped
    return(total_reads)

#Creates a csv with "Sample","SACE_reads","SAPA_mtReads","Mapped_mtReads","Total_mtReads","%SACE","%SAPA" stored in stats/coverage/RefCovPercentageSheets
def gen_csvs(files):
    stats_csv=path+OutRefPercentageFilePath
    percentage_stats=open(stats_csv,'w')
    percentage_stats_writer=csv.writer(percentage_stats)
    i = 1
    for f in files:
        name = os.path.basename(f)
        sample = name.split("_")[0]
        row = []
        if(i==1):
            header_stats=["Sample","SACE_reads","SAPA_Reads","Mapped_Reads","Total_Reads","%SACE","%SAPA"]
            percentage_stats_writer.writerow(header_stats)
            i += 1
        if sample not in files.keys():
            pass
        else:
            stats=get_stats(sample)
            total_reads =get_total(stats)
            total_mapped_reads=float(stats['Mapped'][0])+float(stats['Mapped'][1])
            row_4stats=[sample,stats['Mapped'][0],stats['Mapped'][1],total_mapped_reads,total_reads,float(stats['Mapped'][0])/total_mapped_reads,float(stats['Mapped'][1])/total_mapped_reads]
            percentage_stats_writer.writerow(row_4stats)
    percentage_stats.close()

files = get_files(path)
gen_csvs(files)

print("Sheet with reference coverage percentage of the processed samples can be found in: ", path+OutRefPercentageFilePath)
