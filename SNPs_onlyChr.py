import argparse
import getpass
import pandas as pd
from Bio import SeqIO
import sys

sys.path.append("../generic_functions/")

from generic import header
from generic import bashfile
from generic import get_refs


##Parse arguments
ag = argparse.ArgumentParser()
ag.add_argument("-s", "--samples", required = True, help = "File with samples to work with")
ag.add_argument("-r", "--reference", required = True, help = "Reference genome")
ag.add_argument("--subgenome", required = True, help = "Reference genome")

args = vars(ag.parse_args())
samples = args["samples"]
ref = args["reference"]
sub = args["subgenome"]
##


username = getpass.getuser()
ref_dir= "/mnt/Timina/lmorales/Public/ymez/data/ref"
sge_dir = "/mnt/Timina/lmorales/Public/ymez/bin/SGE/05_vcalling"
tmp_dir = "/mnt/Timina/lmorales/Public/ymez/tmp/05_vcalling"
sh_out = sge_dir + "/" + username + "_SH_vcfOnlychr" + ".sh"
sample_list = pd.read_csv(samples).ID
chromosomes = []
sgs = []
id = get_refs(ref)
np="2"
ram="4G"

for rf in range(len(id["short"])):
    fasta=ref_dir + "/"+ id["large"][rf]+"_v1/fasta/"+ id["large"][rf]+"_v1_allChr.fasta"
    with open(fasta, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            chromosomes.append(record.description)
    for sample in sample_list:
        file = sge_dir + "/" + sample + "_vcf_" + id["short"][rf] + "_Onlychr" + sub +".sge"
        sgs.append(sample + "_vcf_" + id["short"][rf] + "_Onlychr" + sub)
        with open(file,'w') as sge_file:
            header(sge_dir, sample, sge_file, "vcf_" + id["short"][rf] + "_Onlychr" ,username,np,ram)
            print('module load vcftools/0.1.14 htslib/1.9\n###', file =sge_file)
            print('start=$(date +%s.%N)', file = sge_file)
            print(f"vcftools --vcf {tmp_dir}/{sample}_{id['short'][rf]}.g.vcf",end=' ', file=sge_file)
            for chr in chromosomes:
                if chr.split("_")[0] == sub:
                    if chr.split("_")[4] != "2m" and chr.split("_")[4] != "mt":
                        print(f"--chr {chr}",end=' ', file=sge_file)
            print(f"--recode --recode-INFO-all --out {tmp_dir}/{sample}_{id['short'][rf]}.SNP_onlychr_{sub}.g\n",file=sge_file)
            print(f"mv {tmp_dir}/{sample}_{id['short'][rf]}.SNP_onlychr_{sub}.g.recode.vcf {tmp_dir}/{sample}_{id['short'][rf]}.SNP_onlychr_{sub}.g.vcf",file=sge_file)
            print('',file=sge_file)
            print(f"bgzip -f {tmp_dir}/{sample}_{id['short'][rf]}.SNP_onlychr_{sub}.g.vcf",file=sge_file)
            print(f"tabix -f {tmp_dir}/{sample}_{id['short'][rf]}.SNP_onlychr_{sub}.g.vcf.gz",file=sge_file)
            print('',file=sge_file)
            print('duration=$(echo "$(date +%s.%N) - $start" | bc)', file =sge_file)
            print('execution_time=`printf "%.2f seconds" $duration`', file = sge_file)
            print('echo "Script Execution Time: $execution_time"', file =sge_file)
bashfile(sh_out,sgs)
