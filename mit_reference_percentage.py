import os,csv
import subprocess
import re
import pandas as pd
import argparse

ag=argparse.ArgumentParser()
ag.add_argument("-s","--sample",required=True, help= "Absolute path to the csv with the samples to process. Example: mnt/Timina/lmorales/Public/ymez/data/metadata/samplesheet.csv. File must have only an ID column with no header and one sample per line.")
args=vars(ag.parse_args())
samples=args["sample"]

path = "/mnt/Timina/lmorales/Public/ymez/stats/coverage/" #Path to stats/coverage/
script_dir = "/mnt/Timina/lmorales/Public/ymez/bin/scripts/04_coverage/" #Path of the invoked Rscript
id_file =samples #Sample sheet with the samples to be processed
OutRefPercentageFilePath="RefCovPercentageSheets/RefPercentageSheet_" + os.path.basename(samples)[:-4] + ".csv" #Name of the output csv that will be generated at the end of the program, contains mapped, unmapped and mt read % of each sample.

#Function to extract the IDs from the sample sheet. The file should only contain ID per line.
def get_ids(file):
    with open(file,newline="\n") as f:
        reader = csv.reader(f,delimiter=",")
        ids = list(reader)
    return ids

#Function to fetch stat files for each ID.
def get_files(path,Q):
    ids = get_ids(id_file)
    files = {}
    for id in ids:
        files[id[0]] = path + id[0] + "/mit/Q00/" + id[0] + "_mit_perReference.csv"
    return files

#Runs samtools idxstats to get number of mapped and unmapped reads for each reference, stores de information in a dataframe and returns it for further use.
def get_stats(id):
    file = f"/mnt/Timina/lmorales/Public/ymez/data/bam/{id}_mit.rmdup.addgp.bam"
    p = "samtools idxstats " + file
    ps_total=subprocess.Popen(p,shell=True,stdout=subprocess.PIPE)
    ps_total.wait()
    total=str(ps_total.communicate()[0]).split("'")[1]
    total=total.split('\\n')
    del total[-1]
    SAPA=re.sub(r'\\t','.',total[0]).split('.')
    SACE=re.sub(r'\\t','.',total[1]).split('.')
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

#Creates two files: 1) csv with "Sample",'SACE','SAPA',"Mapped","Total" for the plotting R script to use, 2) csv with "Sample","mtSACE_reads","mtSAPA_mtReads","Mapped_mtReads","Total_mtReads","%mtSACE","%mtSAPA" stored in stats/coverage/RefCovPercentageSheets
def gen_csvs(files,Q):
    if Q == "C_p":
        csv_name = path+"perRef.csv"
    else:
        csv_name = path+Q+"_perRef.csv"
    stats_csv=path+OutRefPercentageFilePath
    percentage_stats=open(stats_csv,'w')
    percentage_stats_writer=csv.writer(percentage_stats)
    summary = open(csv_name,"w")
    writer = csv.writer(summary)
    i = 1
    for f in files:
        name = os.path.basename(f)
        sample = name.split("_")[0]
        row = []
        if(i==1):
            header = ["Sample",'SACE','SAPA',"Mapped","Total"]
            writer.writerow(header)
            header_stats=["Sample","mtSACE_reads","mtSAPA_mtReads","Mapped_mtReads","Total_mtReads","%mtSACE","%mtSAPA"]
            percentage_stats_writer.writerow(header_stats)
            i += 1
        if sample not in files.keys():
            pass
        else:
            stats=get_stats(sample)
            total_reads =get_total(stats)
            total_mapped_reads=float(stats['Mapped'][0])+float(stats['Mapped'][1])
            row = [sample,stats['Mapped'][0],stats['Mapped'][1],total_mapped_reads,total_reads]
            row_4stats=[sample,stats['Mapped'][0],stats['Mapped'][1],total_mapped_reads,total_reads,float(stats['Mapped'][0])/total_mapped_reads,float(stats['Mapped'][1])/total_mapped_reads]
            percentage_stats_writer.writerow(row_4stats)
            writer.writerow(row)
    summary.close()
    percentage_stats.close()
    return csv_name

Qs = ["C_p","Q20"]
for Q in Qs: #for cycle to run the plotting process on two mapping qualities: Q00 and Q20
    files = get_files(path,Q)
    name = gen_csvs(files,Q)
    process = ["Rscript",script_dir+"mit_porcentage_reference.R",name,Q] #Calls Rscript to run the code "mit_porcentage_reference.R", which creates the barplots.
    subprocess.call(process)
print("Sheet with mitochondrial reference coverage percentage of the processed samples can be found in: ", path+OutRefPercentageFilePath)
