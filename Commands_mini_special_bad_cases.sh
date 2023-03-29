month="Lote55_08Junio2020_IBT"
out_dir="/scratch/rodrigog/01_Projects/VAMIP/04_fastq/temp/01_IBT/$month"; mkdir -p $out_dir; cd $out_dir
# Concatenate 4 lanes into a single fastq
# in_dir="/share/DiskCoViGen/VAMIP/04_fastq/01_IBT/1000Genomas/Segunda"
# in_dir="/share/DiskCoViGen/VAMIP/04_fastq/01_IBT/$month"
in_dir="/scratch/btaboada/Data/CoV2019/Vigilancia/Lote55_08Junio2020_IBT"
for sam in $(ls $in_dir/*.gz);do echo 'cp '$sam' '$out_dir/'';done >cp_files.sh
nohup bash cp_files.sh &
mv nohup.out FAILED_fastqs_missing_lanes.txt;cat FAILED_fastqs_missing_lanes.txt
cd $out_dir; mkdir -p 00_qsubout 01_Preprocess 02_derep 03_assembly
ls *.gz >files.list
# Filter reads with fastp
echo 'qsub -R y -l h_rt=23:59:59 -N preProcess_IBT_'$month' -o '$out_dir'/00_qsubout -t 1-'$(ls *.gz|wc -l)':2 /share/DiskCoViGen/VAMIP/Programs/preProcess_list_v2.sh '$out_dir' files.list '$out_dir'/01_Preprocess'|bash
a=$(ls 01_Preprocess/*.gz|wc -l);for x in $(grep -m 1 "Done\!" 00_qsubout/preProcess*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done));do echo 'grep -m 1 "\.'$x'$" <(ls 00_qsubout/preProcess*)';done|bash >00_qsubout/FAILED_preProcess.txt;cat 00_qsubout/FAILED_preProcess.txt
ls 01_Preprocess/*.gz|wc -l
# dereplicate with cdhit
for i in $(ls 01_Preprocess/*.fastq.gz);do basename $i;done >$out_dir/01_Preprocess/files_clean.list
echo 'qsub -R y -l h_rt=23:59:59 -N cdHitDup_IBT_'$month' -o '$out_dir'/00_qsubout -t 1-'$(ls 01_Preprocess/*.gz|wc -l)':2 /share/DiskCoViGen/VAMIP/Programs/cdHitDupPair_list_v2.1.sh '$out_dir/01_Preprocess' files_clean.list '$out_dir/02_derep''|bash
# Get which ones fail
a=$(ls 01_Preprocess/*.gz|wc -l);for x in $(grep -m 1 "Done\!" 00_qsubout/cdHitDup*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done));do echo "\.$x$";done|grep -f - <(ls 00_qsubout/cdHitDup*) >00_qsubout/FAILED_cdHitDup.txt;cat 00_qsubout/FAILED_cdHitDup.txt
for i in $(cat 00_qsubout/FAILED_cdHitDup.txt);do grep "Total number of sequences" $i;done|wc -l
# Assembly and analyze changes vs Original SARS-CoV-2 reference
for i in $(ls 02_derep/*.fastq.gz);do basename $i;done >$out_dir/02_derep/files_derep.list
echo 'qsub -R y -l h_rt=23:59:59 -N ConsesusVariant_IBT_'$month' -o '$out_dir'/00_qsubout -t 1-'$(ls 02_derep/*.gz|wc -l)':2 /share/DiskCoViGen/VAMIP/Programs/variantConsePairCall_List_v2.1b_for_bad_cases.sh '$out_dir/02_derep' files_derep.list /share/DiskCoViGen/VAMIP/DB/MN908947.fasta test '$out_dir/03_assembly''|bash
# Get which ones fail
a=$(ls 02_derep/*.gz|wc -l);for x in $(grep -m 1 "Done\!" 00_qsubout/ConsesusVariant*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done));do echo "\.$x$";done|grep -f - <(ls 00_qsubout/ConsesusVariant*) >00_qsubout/FAILED_ConsensusVariant.txt; cat 00_qsubout/FAILED_ConsensusVariant.txt
ls 03_assembly/*.fasta|wc -l
/share/DiskCoViGen/VAMIP/Programs/cleanConsensus_v4b_for_bad_cases.pl
