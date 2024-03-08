#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 18:43:50 2020

@author: luisfer
"""

import sys
import argparse
import getpass

sys.path.append("../generic_functions/")

from generic import getdata
from generic import header
from generic import bashfile


#################################################################################################################################
ag = argparse.ArgumentParser()
ag.add_argument("-s", "--samples", required = True, help = "csv file with the data of the samples to work") 


args = vars(ag.parse_args())
sample = args["samples"]
#project_dir = "/mnt/C0E023C8E023C40E/Projects/MezcalYeast"
#sge_dir = project_dir + "/sge_tests"
project_dir = "/mnt/Timina/lmorales/Public/ymez"
sge_dir = project_dir + "/bin/SGE/02_cleaning"
username = getpass.getuser()
datos = getdata(sample,'Lib_Protocol')
sh_out = sge_dir + "/" + username + "_SH_clean_fastq.sh"
proc = "fastp"
np="16"
ram="1G"
r = dict(datos)

#path, sample, sge, app, user, np, ram

for ids in datos:
    
    file = sge_dir + "/" + ids + ".sge"
    with open(file,'w') as sge_file:
        header(sge_dir, ids, sge_file, proc, username, np, ram)
        print('module load fastp/0.20.0 \n###', file =sge_file)
        print('start=$(date +%s.%N)', file = sge_file)
        print('echo "Start to clean fastq files"', file = sge_file)
        if datos[ids].startswith("PE"):
            str = ids
            lista = ["_R1.fastq.gz","_R1_clean.fastq.gz", "_R2.fastq.gz", "_R2_clean.fastq.gz", "_unpaired_clean.fastq.gz", "_fastp.json", "_fastp.html"]
            str += '% s'
            lista =  [str % i for i in lista]
            print("fastp -i /mnt/Timina/lmorales/Public/ymez/data/fastq/raw/%s -o /mnt/Timina/lmorales/Public/ymez/data/fastq/clean/%s -I /mnt/Timina/lmorales/Public/ymez/data/fastq/raw/%s -O /mnt/Timina/lmorales/Public/ymez/data/fastq/clean/%s --unpaired1 /mnt/Timina/lmorales/Public/ymez/data/fastq/clean/%s -w 16 -y -x -z 9 -j /mnt/Timina/lmorales/Public/ymez/tmp/02_cleaning/%s -h /mnt/Timina/lmorales/Public/ymez/tmp/02_cleaning/%s" % tuple(lista), file = sge_file)            
        elif datos[ids].startswith("SE"):
            str = ids
            lista = ["_R1.fastq.gz","_R1_clean.fastq.gz", "_fastp.json", "_fastp.html"]
            str += '% s'
            lista =  [str % i for i in lista]
            print("fastp -i /mnt/Timina/lmorales/Public/ymez/data/fastq/raw/%s -o /mnt/Timina/lmorales/Public/ymez/data/fastq/clean/%s -w 16 -y -x -z 9 -j /mnt/Timina/lmorales/Public/ymez/tmp/02_cleaning/%s -h /mnt/Timina/lmorales/Public/ymez/tmp/02_cleaning/%s" % tuple(lista), file = sge_file)
        else:
            print("Error in sample %s: %s  is an invalid type of Lib_Protocol. Please use PE or SE in your sampleSheet" % (ids, datos[ids]))
            del r[ids]
            continue

        print('duration=$(echo "$(date +%s.%N) - $start" | bc)', file = sge_file)
        print('execution_time=`printf "%.2f seconds" $duration`',file = sge_file)
        print('echo "Script Execution Time: $execution_time"',file = sge_file)
        print('echo "The scipt ends"', file = sge_file)

#print(list(r.keys()))
bashfile(sh_out, list(r.keys() ) )        
