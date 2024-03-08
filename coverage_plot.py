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
ag.add_argument("-r", "--references", required = True, help = "file with the ids of the genome references to prepare")
ag.add_argument("-s", "--samples", required = True, help = "csv file with the data of the samples to work")

args = vars(ag.parse_args())
refs = args["references"]
samples = args["samples"]

datos = getdata(samples,'SRA')
ids = get_refs(refs)

np="10"
ram="4G"

username = getpass.getuser()

sge_dir = "/mnt/Timina/lmorales/Public/ymez/bin/SGE/04_coverage"
bam_dir = "/mnt/Timina/lmorales/Public/ymez/data/bam"
out_dir = "/mnt/Timina/lmorales/Public/ymez/data/figures/04_coverage"
tmp_dir = "/mnt/Timina/lmorales/Public/ymez/tmp/04_coverage"
script_dir = "/mnt/Timina/lmorales/Public/ymez/bin/scripts/04_coverage"
stats_dir = "/mnt/Timina/lmorales/Public/ymez/stats/coverage"
#%%
sh_out = sge_dir + "/" + username + "_SH_coverage_plot.sh"

sgs = []
for sample in list(datos.keys()):
    for rf in range(len(ids["short"])):
        sgs.append(sample + "_" + ids["short"][rf])
        file = sge_dir + "/" + sample + "_" + ids["short"][rf] + ".sge"
        with open(file,'w') as sge_file:
             header(sge_dir, "04_" + sample + "_" + ids["short"][rf], sge_file, "depth",username,np, ram)
             print('module load htslib/1.2.1 gcc/5.1.0 samtools/1.9 r/3.6.1 \n###', file =sge_file)
             print('start=$(date +%s.%N)', file = sge_file)
             print('bam_dir=%s' % bam_dir, file = sge_file)
             print('out_dir=%s/%s/%s' % (out_dir,sample,ids["short"][rf]), file = sge_file)
             print('tmp_dir=%s' % tmp_dir, file = sge_file)
             print('script_dir=%s' % script_dir, file = sge_file)
             print('stats_dir=%s/%s/%s' % (stats_dir,sample,ids["short"][rf]), file = sge_file)
             print('mkdir -p ${stats_dir}',file=sge_file)
             print('mkdir -p ${out_dir}',file=sge_file)
             for q in ["Q00","Q20","Q30"]:
                 print('mkdir -p ${stats_dir}/%s' % q, file=sge_file)
                 print('mkdir -p ${out_dir}/%s' % q, file=sge_file)
             print('samtools depth -aa ${bam_dir}/%s_%s.rmdup.addgp.bam > ${tmp_dir}/%s_%s.depth' % (sample, ids["short"][rf], sample, ids["short"][rf]), file = sge_file)
             print('samtools depth -aa -Q 20 ${bam_dir}/%s_%s.rmdup.addgp.bam > ${tmp_dir}/%s_%s_Q20.depth' % (sample, ids["short"][rf], sample, ids["short"][rf]), file = sge_file)
             print('samtools depth -aa -Q 30 ${bam_dir}/%s_%s.rmdup.addgp.bam > ${tmp_dir}/%s_%s_Q30.depth' % (sample, ids["short"][rf], sample, ids["short"][rf]), file = sge_file)
             print('R CMD BATCH --no-save --no-restore "--args FILE=\'${tmp_dir}/%s_%s.depth\' SIZE=1000 STEP=1000 FILE_OUT=\'${out_dir}/Q00/%s_%s\' STATS=\'${stats_dir}/Q00\'" ${script_dir}/make_depth_plot.R ${tmp_dir}/%s_%s.Rout' % (sample,ids["short"][rf], sample, ids["short"][rf], sample, ids["short"][rf]), file = sge_file)
             print('R CMD BATCH --no-save --no-restore "--args FILE=\'${tmp_dir}/%s_%s_Q20.depth\' SIZE=1000 STEP=1000 FILE_OUT=\'${out_dir}/Q20/%s_%s_Q20\' STATS=\'${stats_dir}/Q20\'" ${script_dir}/make_depth_plot.R ${tmp_dir}/%s_%s_Q20.Rout' % (sample,ids["short"][rf],  sample, ids["short"][rf],  sample, ids["short"][rf]), file = sge_file)
             print('R CMD BATCH --no-save --no-restore "--args FILE=\'${tmp_dir}/%s_%s_Q30.depth\' SIZE=1000 STEP=1000 FILE_OUT=\'${out_dir}/Q30/%s_%s_Q30\' STATS=\'${stats_dir}/Q30\'" ${script_dir}/make_depth_plot.R ${tmp_dir}/%s_%s_Q30.Rout' % (sample, ids["short"][rf], sample, ids["short"][rf], sample, ids["short"][rf]), file = sge_file)
             print('duration=$(echo "$(date +%s.%N) - $start" | bc)', file =sge_file)
             print('execution_time=`printf "%.2f seconds" $duration`', file = sge_file)
             print('echo "Script Execution Time: $execution_time"', file =sge_file)

bashfile(sh_out, sgs)
