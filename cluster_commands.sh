# folder Blanca: /scratch/btaboada/Data/CoV2019/
# Started on 2021-11-25
# by Rodrigo García-López for Carlos Arias's Viromics Group at IBt, UNAM, Cuernavaca, Mexico as part of the CoViGen-Mex SARS-CoV-2 survillance in Mexico.
# Under GNU GPLv3 license
# Disclamer: This script was written for use with specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.
# The following commands were used to process all fastq from samples diagnosed for SARS-CoV-2 from CoViGen-Mex
# Files have been deposited in an hdd at /share/DiskCoViGen/VAMIP, and were separated by batch (lote)
# The contents are as follows:
# /share/DiskCoViGen/VAMIP
# ├── 04_fastq
# │   ├── 01_IBT
# │   ├── 02_INER
# │   ├── 03_LANGEBIO
# │   └── 04_CIAD
# └── Programs
#     ├── cdHitDupPair_list_v2.sh
#     ├── cleanConsensus_v4.pl
#     ├── preProcess_list_v2.sh
#     └── variantConsePairCall_List_v2.sh

# Processing will be carried in scratch, as it is a faster disk at /scratch/rodrigog/01_Projects/VAMIP/04_fastq/temp

# CIAD
# I'll test everything with the smallest sets, those from 04_CIAD
set="FASTQ_Generation_2023-02-02_21_40_57Z-651299655"
out_dir="/scratch/rodrigog/01_Projects/VAMIP/04_fastq/temp/04_CIAD/$set"
mkdir -p $out_dir; cd $out_dir
in_dir="/share/DiskCoViGen/VAMIP/04_fastq/04_CIAD/2nd_CIAD-2022-11-24/$set"
for i in $(find $in_dir -name "*.gz");do echo cp $i $out_dir/;done >cp_files.sh
nohup bash cp_files.sh &

mv nohup.out FAILED_fastqs_missing_lanes.txt; cat FAILED_fastqs_missing_lanes.txt
rename "_L001_" "_" *.gz
rename "_001.fastq" ".fastq" *.gz
# rename ".R1" '_R1' *.gz # Names from CIAD files have a uncommon .R1 separator, which is not compatible with our current scripts
# rename ".R2" '_R2' *.gz
mkdir -p 00_qsubout 01_Preprocess 02_derep 03_assembly
ls *.gz >files.list
# Filter reads with fastp
# out_dir="/scratch/rodrigog/01_Projects/VAMIP/04_fastq/temp/04_CIAD"
echo 'qsub -R y -l h_rt=23:59:59 -l h_vmem=1G -pe thread 2 -N preProcess_CIAD -o '$out_dir'/00_qsubout -t 1-'$(ls *.gz|wc -l)':2 -tc 20 /share/DiskCoViGen/VAMIP/Programs/preProcess_list_v2.sh '$out_dir' files.list '$out_dir'/01_Preprocess'|bash

a=$(ls *.gz|wc -l);for x in $(grep -m 1 "Done\!" 00_qsubout/preProcess*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done));do echo 'grep -m 1 "\.'$x'$" <(ls 00_qsubout/preProcess*)';done|bash >00_qsubout/FAILED_preProcess.txt;cat 00_qsubout/FAILED_preProcess.txt
ls 01_Preprocess/*.gz|wc -l
# dereplicate with cdhit
for i in $(ls 01_Preprocess/*.fastq.gz);do basename $i;done >$out_dir/01_Preprocess/files_clean.list
echo 'qsub -R y -l h_rt=23:59:59 -l h_vmem=700M -pe thread 1 -N cdHitDup_CIAD -o '$out_dir'/00_qsubout -t 1-'$(ls 01_Preprocess/*.gz|wc -l)':2 -tc 20 /share/DiskCoViGen/VAMIP/Programs/cdHitDupPair_list_v2.1.sh '$out_dir/01_Preprocess' files_clean.list '$out_dir/02_derep''|bash

# Get which ones fail
a=$(ls 01_Preprocess/*.gz|wc -l);grep -m 1 "Done\!" 00_qsubout/cdHitDup*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done)|grep -f - <(ls 00_qsubout/cdHitDup*) >00_qsubout/FAILED_cdHitDup.txt;cat 00_qsubout/FAILED_cdHitDup.txt
# Assembly and analyze changes vs Original SARS-CoV-2 reference
for i in $(ls 02_derep/*.fastq.gz);do basename $i;done >$out_dir/02_derep/files_derep.list
echo 'qsub -R y -l h_rt=23:59:59 -l h_vmem=1G -N ConsesusVariant_CIAD -o '$out_dir'/00_qsubout -t 1-'$(ls 02_derep/*.gz|wc -l)':2 -tc 20 /share/DiskCoViGen/VAMIP/Programs/variantConsePairCall_List_v2.1.sh '$out_dir/02_derep' files_derep.list /share/DiskCoViGen/VAMIP/DB/MN908947.fasta genome '$out_dir/03_assembly''|bash

# Get which ones fail
a=$(ls 01_Preprocess/*.gz|wc -l);grep -m 1 "Done\!" 00_qsubout/ConsesusVariant*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done)|grep -f - <(ls 00_qsubout/ConsesusVariant*) >00_qsubout/FAILED_ConsensusVariant.txt; cat 00_qsubout/FAILED_ConsensusVariant.txt
ls 03_assembly/*.fasta|wc -l
rm -R 01_Preprocess/; rm -R 02_derep/

# LANGEBIO
month="Lote63_26Julio2022_LANGEBIO"
out_dir="/scratch/rodrigog/01_Projects/VAMIP/04_fastq/temp/03_LANGEBIO/$month"; mkdir -p $out_dir; cd $out_dir
# Concatenate 4 lanes into a single fastq
in_dir="/scratch/rodrigog/01_Projects/VAMIP/04_fastq/temp/$month"
# in_dir="/scratch/rodrigog/01_Projects/VAMIP/04_fastq/temp/$month/reads"
# for sam in $(ls $in_dir|grep fastq.gz|sed "s/_L00.*//"|sort|uniq);do for pair in R1 R2;do printf "cat"; for lane in L001 L002 L003 L004; do printf ' '$in_dir/$sam'_'$lane'_'$pair'_001.fastq.gz ';done;printf ' >'$out_dir'/'$sam'_'$pair'.fastq.gz\n';done;done >cat_files.sh
for sam in $(find $in_dir -name "*.fastq.gz");do base=$(basename $sam); printf "cat $sam >>$out_dir/";echo $base|sed "s/_L00[1-4]_\(R[ls 12]\)_00[12]\.fastq\.gz/_\1.fastq.gz/";done|sort >cat_files.sh
nohup bash cat_files.sh &
mv nohup.out FAILED_fastqs_missing_lanes.txt; cat FAILED_fastqs_missing_lanes.txt
cd $out_dir; mkdir -p 00_qsubout 01_Preprocess 02_derep 03_assembly
ls *.gz|sort >files.list
# Filter reads with fastp
echo 'qsub -R y -l h_rt=23:59:59 -l h_vmem=1200M -pe thread 2 -N preProcess_LANGEBIO_'$(basename $out_dir)' -o '$out_dir'/00_qsubout -t 1-'$(ls *.gz|wc -l)':2 -tc 30 /share/DiskCoViGen/VAMIP/Programs/preProcess_list_v2.sh '$out_dir' files.list '$out_dir'/01_Preprocess'|bash
a=$(ls *.gz|wc -l);for x in $(grep -m 1 "Done\!" 00_qsubout/preProcess*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done));do echo 'grep -m 1 "\.'$x'$" <(ls 00_qsubout/preProcess*)';done|bash >00_qsubout/FAILED_preProcess.txt;cat 00_qsubout/FAILED_preProcess.txt
ls 01_Preprocess/*.gz|wc -l
# # dereplicate with cdhit
for i in $(ls 01_Preprocess/*.fastq.gz);do basename $i;done >$out_dir/01_Preprocess/files_clean.list
echo 'qsub -R y -l h_rt=23:59:59 -l h_vmem=1100M -pe thread 2 -N cdHitDup_LANGEBIO_'$(basename $out_dir)' -o '$out_dir'/00_qsubout -t 1-'$(ls 01_Preprocess/*.gz|wc -l)':2 -tc 30 /share/DiskCoViGen/VAMIP/Programs/cdHitDupPair_list_v2.1.sh '$out_dir/01_Preprocess' files_clean.list '$out_dir/02_derep''|bash
# Get which ones fail
a=$(ls 01_Preprocess/*.gz|wc -l);grep -m 1 "Done\!" 00_qsubout/cdHitDup*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done)|grep -f - <(ls 00_qsubout/cdHitDup*) >00_qsubout/FAILED_cdHitDup.txt;cat 00_qsubout/FAILED_cdHitDup.txt
# Assembly and analyze changes vs Original SARS-CoV-2 reference
for i in $(ls 02_derep/*.fastq.gz);do basename $i;done >$out_dir/02_derep/files_derep.list
echo 'qsub -R y -l h_rt=23:59:59 -l h_vmem=1G -pe thread 2 -N ConsesusVariant_LANGEBIO_'$(basename $out_dir)' -o '$out_dir'/00_qsubout -t 1-'$(ls 02_derep/*.gz|wc -l)':2 -tc 30 /share/DiskCoViGen/VAMIP/Programs/variantConsePairCall_List_v2.1.sh '$out_dir/02_derep' files_derep.list /share/DiskCoViGen/VAMIP/DB/MN908947.fasta test '$out_dir/03_assembly''|bash
# Get which ones fail
a=$(ls 02_derep/*.gz|wc -l);grep -m 1 "Done\!" 00_qsubout/ConsesusVariant*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done)|grep -f - <(ls 00_qsubout/ConsesusVariant*) >00_qsubout/FAILED_ConsensusVariant.txt; cat 00_qsubout/FAILED_ConsensusVariant.txt
ls 03_assembly/*.fasta|wc -l
rm -R 01_Preprocess/; rm -R 02_derep/

# Append the Batch at the beggining of each name
ls *fastq.gz|sed -e 's/\(.*\)/mv \1 L080_\1/'|bash

# IBT
month="Lote79_30Enero2023_IBT"
out_dir="/scratch/rodrigog/01_Projects/VAMIP/04_fastq/temp/01_IBT/$month"; mkdir -p $out_dir; cd $out_dir
# Concatenate 4 lanes into a single fastq
# in_dir="/share/DiskCoViGen/VAMIP/04_fastq/01_IBT/1000Genomas/Segunda"
in_dir="/scratch/btaboada/Data/CoV2019/Vigilancia/$month"
# in_dir="/scratch/btaboada/Data/CoV2019/Vigilancia/Lote55_08Junio2020_IBT"
for sam in $(ls $in_dir/*.gz);do echo 'cp '$sam' '$out_dir/'';done >cp_files.sh
nohup bash cp_files.sh &
mv nohup.out FAILED_fastqs_missing_lanes.txt;cat FAILED_fastqs_missing_lanes.txt
cd $out_dir; mkdir -p 00_qsubout 01_Preprocess 02_derep 03_assembly
ls *.gz >files.list
# Filter reads with fastp
echo 'qsub -R y -l h_rt=23:59:59 -l h_vmem=1G -pe thread 2 -N preProcess_IBT_'$month' -o '$out_dir'/00_qsubout -t 1-'$(ls *.gz|wc -l)':2 -tc 20 /share/DiskCoViGen/VAMIP/Programs/preProcess_list_v2.sh '$out_dir' files.list '$out_dir'/01_Preprocess'|bash
a=$(ls *.gz|wc -l);for x in $(grep -m 1 "Done\!" 00_qsubout/preProcess*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done));do echo 'grep -m 1 "\.'$x'$" <(ls 00_qsubout/preProcess*)';done|bash >00_qsubout/FAILED_preProcess.txt;cat 00_qsubout/FAILED_preProcess.txt
ls 01_Preprocess/*.gz|wc -l
# dereplicate with cdhit
for i in $(ls 01_Preprocess/*.fastq.gz);do basename $i;done >$out_dir/01_Preprocess/files_clean.list
echo 'qsub -R y -l h_rt=23:59:59 -l h_vmem=700M -pe thread 1 -N cdHitDup_IBT_'$month' -o '$out_dir'/00_qsubout -t 1-'$(ls 01_Preprocess/*.gz|wc -l)':2 -tc 20 /share/DiskCoViGen/VAMIP/Programs/cdHitDupPair_list_v2.1.sh '$out_dir/01_Preprocess' files_clean.list '$out_dir/02_derep''|bash
# Get which ones fail
a=$(ls 01_Preprocess/*.gz|wc -l);for x in $(grep -m 1 "Done\!" 00_qsubout/cdHitDup*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done));do echo "\.$x$";done|grep -f - <(ls 00_qsubout/cdHitDup*) >00_qsubout/FAILED_cdHitDup.txt;cat 00_qsubout/FAILED_cdHitDup.txt
for i in $(cat 00_qsubout/FAILED_cdHitDup.txt);do grep "Total number of sequences" $i;done|wc -l
# Assembly and analyze changes vs Original SARS-CoV-2 reference
for i in $(ls 02_derep/*.fastq.gz);do basename $i;done >$out_dir/02_derep/files_derep.list
echo 'qsub -R y -l h_rt=23:59:59 -l h_vmem=1G -N ConsesusVariant_IBT_'$month' -o '$out_dir'/00_qsubout -t 1-'$(ls 02_derep/*.gz|wc -l)':2 -tc 20 /share/DiskCoViGen/VAMIP/Programs/variantConsePairCall_List_v2.1.sh '$out_dir/02_derep' files_derep.list /share/DiskCoViGen/VAMIP/DB/MN908947.fasta test '$out_dir/03_assembly''|bash
# Get which ones fail
a=$(ls 02_derep/*.gz|wc -l);for x in $(grep -m 1 "Done\!" 00_qsubout/ConsesusVariant*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done));do echo "\.$x$";done|grep -f - <(ls 00_qsubout/ConsesusVariant*) >00_qsubout/FAILED_ConsensusVariant.txt; cat 00_qsubout/FAILED_ConsensusVariant.txt
ls 03_assembly/*.fasta|wc -l

# INER
month="Joel_basespace"
out_dir="/scratch/rodrigog/01_Projects/VAMIP/04_fastq/temp/02_INER/$month"; mkdir -p $out_dir; cd $out_dir
# Concatenate 4 lanes into a single fastq
in_dir="/share/DiskCoViGen/VAMIP/04_fastq/02_INER/$month"
for sam in $(find $in_dir -name "*.fastq.gz");do base=$(basename $sam); printf "cat $sam >>$out_dir/";echo $base|sed "s/_L00[1-4]_\(R[12]\)_00[12]\.fastq\.gz/_\1.fastq.gz/";done >cat_files.sh # this assumes there are no repeated prefixes. otherwise, it will concatenate incorrect files
nohup bash cat_files.sh &
mv nohup.out FAILED_fastqs_missing_lanes.txt; cat FAILED_fastqs_missing_lanes.txt
mkdir -p 00_qsubout 01_Preprocess 02_derep 03_assembly
ls *.gz|sort >files.list
# Filter reads with fastp
echo 'qsub -R y -l h_rt=23:59:59 -l h_vmem=1G -pe thread 2 -N preProcess_INER_'$month' -o '$out_dir'/00_qsubout -t 1-'$(ls *.gz|wc -l)':2 -tc 20 /share/DiskCoViGen/VAMIP/Programs/preProcess_list_v2.sh '$out_dir' files.list '$out_dir'/01_Preprocess'|bash
a=$(ls *.gz|wc -l);for x in $(grep -m 1 "Done\!" 00_qsubout/preProcess*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done));do echo 'grep -m 1 "\.'$x'$" <(ls 00_qsubout/preProcess*)';done|bash >00_qsubout/FAILED_preProcess.txt;cat 00_qsubout/FAILED_preProcess.txt
ls 01_Preprocess/*.gz|wc -l
# dereplicate with cdhit
for i in $(ls 01_Preprocess/*.fastq.gz);do basename $i;done >$out_dir/01_Preprocess/files_clean.list
echo 'qsub -R y -l h_rt=23:59:59 -l h_vmem=700M -pe thread 2 -N cdHitDup_INER_'$month' -o '$out_dir'/00_qsubout -t 1-'$(ls 01_Preprocess/*.gz|wc -l)':2 -tc 20 /share/DiskCoViGen/VAMIP/Programs/cdHitDupPair_list_v2.1.sh '$out_dir/01_Preprocess' files_clean.list '$out_dir/02_derep''|bash
# Get which ones fail
a=$(ls 01_Preprocess/*.gz|wc -l);for x in $(grep -m 1 "Done\!" 00_qsubout/cdHitDup*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done));do echo 'grep -m 1 "\.'$x'$" <(ls 00_qsubout/cdHitDup*)';done|bash >00_qsubout/FAILED_cdHitDup.txt;cat 00_qsubout/FAILED_cdHitDup.txt
# Assembly and analyze changes vs Original SARS-CoV-2 reference
for i in $(ls 02_derep/*.fastq.gz);do basename $i;done >$out_dir/02_derep/files_derep.list
echo 'qsub -R y -l h_rt=23:59:59 -l h_vmem=1G -pe thread 2 -N ConsesusVariant_INER_'$month' -o '$out_dir'/00_qsubout -t 1-'$(ls 02_derep/*.gz|wc -l)':2 -tc 20 /share/DiskCoViGen/VAMIP/Programs/variantConsePairCall_List_v2.1.sh '$out_dir/02_derep' files_derep.list /share/DiskCoViGen/VAMIP/DB/MN908947.fasta test '$out_dir/03_assembly''|bash
# Get which ones fail
a=$(ls 02_derep/*.gz|wc -l);for x in $(grep -m 1 "Done\!" 00_qsubout/ConsesusVariant*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done));do echo 'grep -m 1 "\.'$x'$" <(ls 00_qsubout/ConsesusVariant*)';done|bash >00_qsubout/FAILED_ConsensusVariant.txt;cat 00_qsubout/FAILED_ConsensusVariant.txt
ls 03_assembly/*.fasta|wc -l


# Special case (Joel's 379 sequences) Pending identification in GISAID
cd /scratch/btaboada/Data/CoV2019/Vigilancia/BasespacesJoel/
for i in $(find COVID* -name "*.gz");do echo cp $i /scratch/rodrigog/01_Projects/VAMIP/04_fastq/temp/02_INER/;done|bash

/scratch/rodrigog/01_Projects/VAMIP/04_fastq/temp/02_INER/

#### Fix variant calling (now: 10 or more) ###
# The script variantConsePairCall_List_v2.1.sh had the following call for ivar:
# samtools mpileup -aa -B -A -d 1000000 -Q 20 $5/$nSe".srt.bam" | ivar consensus  -q 20 -t 0 -m 5 -n N -p $5/$nSe".ivar"
# We now want to use -m 10 instead (which means repeating the ConsensusVariant step)
# To do this, we need a bam file. Due to time constraints, we rather work with the bam files onward and thus have created a patch based on the variantConsePairCall_List_v2.1.sh file but only carrying out the ivar steps
# This was all carried @ /scratch/rodrigog/01_Projects/VAMIP/04_fastq/fix

# DEPRECATED: This was originally carried out by copying bam files in diskcoViGen into a faster disk, processing them, then moving results back to DiskCoViGen. But it's faster to just do it directly on DiskCoViGen.
# # # cd /scratch/rodrigog/01_Projects/VAMIP/04_fastq/fix
# # # mkdir 01_IBt 02_INER 03_LANGEBIO 04_CIAD
# # # #### IBT:
# # # month="Lote21_IBT_25Agosto2021"
# # # out_dir="/scratch/rodrigog/01_Projects/VAMIP/04_fastq/fix/01_IBT/$month"; mkdir -p $out_dir/00_qsubout $out_dir/02_derep $out_dir/03_assembly; cd $out_dir
# # # cp /share/DiskCoViGen/VAMIP/05_variants/01_IBT/$month/02_derep/files_derep.list 02_derep
# # # nohup time cp /share/DiskCoViGen/VAMIP/05_variants/01_IBT/$month/03_assembly/*.bam 03_assembly/ &
# # # # Map and analyze changes vs Original SARS-CoV-2 reference
# # # echo 'qsub -R y -l h_rt=23:59:59 -l h_vmem=1G -N fix_ivar_LANGEBIO_'$month' -o '$out_dir'/00_qsubout -t 1-'$(cat 02_derep/files_derep.list|wc -l)':2 /share/DiskCoViGen/VAMIP/Programs/ivar_only_v1.sh '$out_dir/02_derep' files_derep.list /share/DiskCoViGen/VAMIP/DB/MN908947.fasta test '$out_dir/03_assembly''|bash
# # # # Get which ones fail
# # # a=$(cat 02_derep/files_derep.list|wc -l);for x in $(grep -m 1 "Done\!" 00_qsubout/fix_ivar*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done));do echo "\.$x$";done|grep -f - <(ls 00_qsubout/fix_ivar*) >00_qsubout/FAILED_Fixivar.txt; cat 00_qsubout/FAILED_Fixivar.txt
# # # rm 03_assembly/*.bam
# # # nohup cp 03_assembly/ /share/DiskCoViGen/VAMIP/05_variants/01_IBT/$month/03_assembly/ &


#### This is done directly in the output folders in DiskCoViGen:
month="Noviembre25"
ins="03_LANGEBIO"
out_dir="/share/DiskCoViGen/VAMIP/05_variants/$ins/$month"; cd $out_dir
rm 03_assembly/*.gz # remove bad tables
# Map and analyze changes vs Original SARS-CoV-2 reference
echo 'qsub -R y -l h_rt=23:59:59 -l h_vmem=1G -N fix_ivar_LANGEBIO_'$month' -o '$out_dir'/00_qsubout -t 1-'$(cat 02_derep/files_derep.list|wc -l)':2 /share/DiskCoViGen/VAMIP/Programs/ivar_only_v1.sh '$out_dir/02_derep' files_derep.list /share/DiskCoViGen/VAMIP/DB/MN908947.fasta test '$out_dir/03_assembly''|bash
# Get which ones fail
a=$(cat 02_derep/files_derep.list|wc -l);for x in $(grep -m 1 "Done\!" 00_qsubout/fix_ivar*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done));do echo "\.$x$";done|grep -f - <(ls 00_qsubout/fix_ivar*) >00_qsubout/FAILED_Fixivar.txt; cat 00_qsubout/FAILED_Fixivar.txt
ls 03_assembly/*.fasta|wc -l

# CIAD
out_dir="/share/DiskCoViGen/VAMIP/05_variants/04_CIAD"; cd $out_dir
rm 03_assembly/*.gz # remove bad tables
# Assembly and analyze changes vs Original SARS-CoV-2 reference
echo 'qsub -R y -l h_rt=23:59:59 -l h_vmem=1G -N fix_ivar_CIAD -o '$out_dir'/00_qsubout -t 1-'$(cat 02_derep/files_derep.list|wc -l)':2 /share/DiskCoViGen/VAMIP/Programs/ivar_only_v1.sh '$out_dir/02_derep' files_derep.list /share/DiskCoViGen/VAMIP/DB/MN908947.fasta test '$out_dir/03_assembly''|bash
# Get which ones fail
a=$(cat 02_derep/files_derep.list|wc -l);for x in $(grep -m 1 "Done\!" 00_qsubout/fix_ivar*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done));do echo "\.$x$";done|grep -f - <(ls 00_qsubout/fix_ivar*) >00_qsubout/FAILED_Fixivar.txt; cat 00_qsubout/FAILED_Fixivar.txt
ls 03_assembly/*.fasta|wc -l


cd /LUSTRE/usuario/aherrera/covid/enero17-22/reads/
rsync -hvrPtLK *.gz rodrigog@teopanzolco.ibt.unam.mx:/scratch/rodrigog/01_Projects/VAMIP/04_fastq/temp/marzo-22/reads

#### RENAME FILES ###
# Most current names are actually the original names for illumina sequencers. In order to use the metadata, we need a common identifier (We are using the IMSS [folio] id).
# 01_IBT:
# We have xref tables of seq id x folio in tables found @ /share/DiskCoViGen/VAMIP/04_fastq/01_IBT/IdSampleRelations
# I wrote a simple oneliner to process this (in a script called /share/DiskCoViGen/VAMIP/04_fastq/01_IBT/rename_template.sh, which includes the following lines:
find $dir -name "*.fast*" -or -name "*_test.*"|sort >temp.txt; paste -d " " temp.txt <(sed 's/\(.*\)\t\(.*\)/s\/\\\/\1_S\/\\\/\2_S\//' /share/DiskCoViGen/VAMIP/04_fastq/01_IBT/IdSampleRelations/$dict|sed -f - temp.txt|sed -e "s/_001\./\./" -e "s/_001_Q/_Q/" -e "s/_test//")|sed -e "s/^/mv /" >rename_samples_$dir.sh
# It is run as follows and generates a new bash script with renaming instructions, where the first is the folder for finding the files, the second is the dictionary for replacing sample names. # The ones commented have a correct name, and thus need no further changes
cd /share/DiskCoViGen/VAMIP/04_fastq/01_IBT
bash rename_template.sh Lote01_IBT_17Febrero2021 Lote01_IBT_17Feb2021_SampleRelation.txt
bash rename_template.sh Lote04_IBT_17Marzo2021 Lote04_IBT_17Mar2021_SampleRelation.txt
bash rename_template.sh Lote07_IBT_14Abril2021 Lote07_IBT_14Abril2021_SampleRelation.txt
bash rename_template.sh Lote10_IBT_12Mayo2021 Lote10_IBT_12Mayo2021_SampleRelation.txt
# bash rename_template.sh Lote10_IBT_12Mayo2021_P2 Lote10_IBT_12Mayo2021_SampleRelationP2.txt
bash rename_template.sh Lote13_IBT_09Junio2021 Lote13_IBT_09Jun2021_SampleRelation.txt
bash rename_template.sh Lote14_IBT_16Junio2021 Lote14_IBT_16Jun2021_SampleRelation.txt
# bash rename_template.sh Lote18_IBT_28Julio2021 Lote18_IBT_28Jul2021_SampleRelationMor.txt ## Requires adding date
# bash rename_template.sh Lote21_IBT_25Agosto2021 Lote21_IBT_25Agosto2021_SampleRelation.txt
bash rename_template.sh Lote24_IBT_22Septiembree2021 Lote25_IBT_22Sept2021_SampleRelation.txt
# bash rename_template.sh Lote27_IBT_20Octubre2021 Lote27_IBT_20Oct2021_SampleRelation.txt
# bash rename_template.sh Lote30_IBT_17Noviembre2021
# bash rename_template.sh Lote33_IBT_08Diciembre2021
# bash rename_template.sh Lote36_06Enero2022_IBT
bash rename_template.sh Lote39_19Enero2022_IBT Lote39_19Enero2022_IBT_SampleRelation.tsv

# Batches 30 onwards need a different strategy
for i in $(find Lote30_IBT_17Noviembre2021/ -name "*R1*.gz"|sort);do base=$(basename $i); echo $base|sed -e "s/_.*//" -e "s/\([0-9]*\)/\1\t202101\1/";done|grep -vP "^\t2021" >IdSampleRelations/Lote30_IBT_17Noviembre2021_SampleRelation.txt
for i in $(find Lote33_IBT_08Diciembre2021/ -name "*R1*.gz"|sort);do base=$(basename $i); echo $base|sed -e "s/_.*//" -e "s/\([0-9]*\)/\1\t202101\1/";done|grep -vP "^\t2021" >IdSampleRelations/Lote33_IBT_08Diciembre2021_SampleRelation.txt
for i in $(find Lote36_06Enero2022_IBT/ -name "*R1*.gz"|sort);do base=$(basename $i); echo $base|sed -e "s/IBt_//" -e "s/_.*//" -e "s/\([0-9]*\)/IBt_\1\t202101\1/" -e "s/202101000/202201000/";done|grep -vP "IBt_\t" >IdSampleRelations/Lote36_06Enero2022_IBT_SampleRelation.txt # Note that 2022 items are still in the hundreads

# The script was also adapted for processing all 02_derep and 03_assembly items
cd /share/DiskCoViGen/VAMIP/05_variants/01_IBT
bash rename_template.sh Lote01_IBT_17Febrero2021 Lote01_IBT_17Feb2021_SampleRelation.txt
bash rename_template.sh Lote04_IBT_17Marzo2021 Lote04_IBT_17Mar2021_SampleRelation.txt
bash rename_template.sh Lote07_IBT_14Abril2021 Lote07_IBT_14Abril2021_SampleRelation.txt
bash rename_template.sh Lote10_IBT_12Mayo2021 Lote10_IBT_12Mayo2021_SampleRelation.txt
# bash rename_template.sh Lote10_IBT_12Mayo2021_P2 Lote10_IBT_12Mayo2021_SampleRelationP2.txt
bash rename_template.sh Lote13_IBT_09Junio2021 Lote13_IBT_09Jun2021_SampleRelation.txt
bash rename_template.sh Lote14_IBT_16Junio2021 Lote14_IBT_16Jun2021_SampleRelation.txt
# bash rename_template.sh Lote18_IBT_28Julio2021 Lote18_IBT_28Jul2021_SampleRelationMor.txt ## Requires adding date
# bash rename_template.sh Lote21_IBT_25Agosto2021 Lote21_IBT_25Agosto2021_SampleRelation.txt
bash rename_template.sh Lote24_IBT_22Septiembre2021 Lote25_IBT_22Sept2021_SampleRelation.txt
# bash rename_template.sh Lote27_IBT_20Octubre2021 Lote27_IBT_20Oct2021_SampleRelation.txt
bash rename_template.sh Lote30_IBT_17Noviembre2021 Lote30_IBT_17Noviembre2021_SampleRelation.txt
bash rename_template.sh Lote33_IBT_08Diciembre2021 Lote33_IBT_08Diciembre2021_SampleRelation.txt
bash rename_template.sh Lote36_06Enero2022_IBT Lote36_06Enero2022_IBT_SampleRelation.txt
# Batches 37 and 43 need no change but the "test" part in 03_assembly files should be removed
paste <(ls *.gz|sed 's/_.*//') <(ls *.gz|sed 's/_.*//') >/share/DiskCoViGen/VAMIP/05_variants/IdSampleRelations/01_IBT/Lote37_IBT_12Enero2022_SampleRelations.txt
bash rename_template.sh Lote37_IBT_12Enero2022 Lote37_IBT_12Enero2022_SampleRelations.txt
paste <(ls *.gz|sed 's/_.*//') <(ls *.gz|sed 's/_.*//') >/share/DiskCoViGen/VAMIP/05_variants/IdSampleRelations/01_IBT/Lote43_IBT_16Febrero2021_SampleRelations.txt
bash rename_template.sh Lote43_IBT_16Febrero2021 Lote43_IBT_16Febrero2021_SampleRelations.txt
# Now, just execute all bash files (I did this one by one after checking everything was fine in the renaming instructions (check for missing items). Expected errors were left as part of this process to track which ones were not altered (as they weren't present in the samples).
#NOTE: To revert any change, use paste <(cut -d " " -f 1 rename_samples_Lote33_IBT_08Diciembre2021.sh) <(cut -d " " -f 3 rename_samples_Lote33_IBT_08Diciembre2021.sh) <(cut -d " " -f 2 rename_samples_Lote33_IBT_08Diciembre2021.sh)
# Lote46 and Lote47 require no changes. String "tests" originally marks unprocessed (prior to renamming) files so we need to remove this from these two batches.
cd /share/DiskCoViGen/VAMIP/05_variants/01_IBT
for i in $(find Lote46_IBT_09Marzo2022 -name "*test*"); do echo 'mv '$i' '$(echo $i|sed "s/_test//")'';done >remove_test_flag_Lote46.sh
for i in $(find Lote47_IBT_23Marzo2022 -name "*test*"); do echo 'mv '$i' '$(echo $i|sed "s/_test//")'';done >remove_test_flag_Lote47.sh

# 03_LANGEBIO
# A set of renaming tables (xref id to IMSS folio) were added to /share/DiskCoViGen/VAMIP/05_variants/IdSampleRelations/03_LANGEBIO
cd /share/DiskCoViGen/VAMIP/05_variants/03_LANGEBIO
# The script was tuned to work with LANGEBIO's files
bash rename_template.sh Marzo1 Corrida01-Marzo2021.tsv
bash rename_template.sh marzo Corrida01-marzo2021_AltB.tsv # Most of these are repeated, I will keep only those that were changed by deleting like this (from folder marzo): for i in $(find . -name "0*");do echo 'rm '$i'';done|bash
bash rename_template.sh abril Corrida02-Abril2021.tsv
bash rename_template.sh mayo Corrida03-Mayo2021.tsv
bash rename_template.sh junio Corrida04-Junio2021.tsv
bash rename_template.sh julio1 Corrida05-Julio2021.tsv
bash rename_template.sh julio2 Corrida06-Julio2021.tsv
bash rename_template.sh agosto Corrida07-Agosto2021.tsv
bash rename_template.sh septiembre Corrida08-Septiembre2021.tsv
bash rename_template.sh octubre Corrida09-Octubre2021.tsv
bash rename_template.sh noviembre Corrida10-Noviembre2021.tsv
bash rename_template.sh diciembre12 Corrida11-Diciembre2021.tsv
bash rename_template.sh diciembre22 Corrida12-Diciembre2021.tsv
bash rename_template.sh enero10-22 Corrida13-Enero2022.tsv
bash rename_template.sh enero17-22 Corrida14-Enero2022.tsv
bash rename_template.sh febrero01-22 Corrida15-Febrero2022.tsv
bash rename_template.sh febrero08-22 Corrida16-Febrero2022.tsv
bash rename_template.sh marzo-22 Corrida17-Marzo2022.tsv
bash rename_template.sh abril29-2022 Corrida18-Abril2022.tsv
bash rename_template.sh mayo31-22 Corrida19-Mayo2022.tsv # table from Corrida20 has the same items.
bash rename_template.sh AH1COV2SSr040 Corrida21-Junio22-2022.tsv
bash rename_template.sh AH1COV2SSr041 Corrida22-Junio30-2022.tsv
bash rename_template.sh AH1COV2SSr042 Corrida23-Julio29-2022.tsv # MISSING # UPDATE 2023-03-15: Ok now
bash rename_template.sh AH1COV2SSr043 Corrida24-Agosto16-2022.tsv
bash rename_template.sh AH1COV2SSr044 Corrida25-Septiembre15-2022.tsv
bash rename_template.sh AH1COV2SSr045 Corrida26-Octubre4-2022.tsv

#Some pending items:
bash rename_template.sh Lote60_LANGEBIO_29Junio2022 Corrida22-Junio30-2022.tsv



for i in $(find . -name "*.gz"|grep -v L54); do echo 'rm '$i'';done|bash # Some items are not part of the set, we remove those not renaming correctly as a proxy


# All renaming schemes are found in bash files which are used to trace back original names and should be executed to carry out renaming

# And rename batches' names in directories
mv Marzo1 Lote03_LANGEBIO_10Marzo2021
mv marzo Lote03_LANGEBIO_10Marzo2021_P2
mv abril Lote06_LANGEBIO_07Abril2021
mv mayo Lote09_LANGEBIO_05Mayo2021
mv junio Lote12_LANGEBIO_02Junio2021
mv julio1 Lote15_LANGEBIO_30Junio2021
mv julio2 Lote16_LANGEBIO_07Julio2021
mv agosto Lote20_LANGEBIO_11Agosto2021
mv septiembre Lote23_LANGEBIO_09Septiembre2021
mv octubre Lote26_LANGEBIO_06Octubre2021
mv noviembre Lote29_LANGEBIO_03Noviembre2021
mv diciembre12 Lote32_LANGEBIO_01Diciembre2021
mv diciembre22 Lote34_LANGEBIO_15Diciembre2021
mv enero10-22 Lote35_LANGEBIO_05Enero2022
mv enero17-22 Lote38_LANGEBIO_12Enero2022
mv febrero01-22 Lote40_LANGEBIO_26Enero2022 # There are some faulty names here starting with "2.022", so we fix these. These also have a bad E11 suffix as they were written as scientific notation
find . -name "*E11_*" >temp; for i in $(cat temp);do n=$(echo $i|sed -e 's/\/2\.02/\/202/' -e 's/E11_/_/');echo 'mv '$i' '$n'';done|bash
mv febrero08-22 Lote41_LANGEBIO_02Febrero2022
mv marzo-22 Lote45_LANGEBIO_02Marzo2022
mv abril29-2022 Lote50_LANGEBIO_20Abril2022
mv mayo31-22 Lote54_LANGEBIO_25Mayo2022
mv AH1COV2SSr040 Lote57_LANGEBIO_21Junio2022 # (junio 2022)
mv AH1COV2SSr041 Lote60_LANGEBIO_29Junio2022 # (junio 2022)
mv AH1COV2SSr042 Lote63_26Julio2022_LANGEBIO # (julio 2022)
mv AH1COV2SSr043 Lote66_LANGEBIO16_Agosto2022 # (agosto 2022)
mv AH1COV2SSr044 Lote69_LANGEBIO_06Septiembre2022 # (septiembre 2022)
mv AH1COV2SSr045 Lote72_LANGEBIO_04Octubre2022 # (octubre 2022) # Carried out by Nelly, no longer required







# NOTE ERROR: I missed a "/" in the sed replacement (already fixed that line above) and files were mistakenly moved to the batch main folder instead appending the name of the folder as part of the filename. Some items were missplaced as a result, having names such as 02_derep202201003188_UIMY_S144_R2_Q_Dp.fastq.gz and 03_assembly202201002368_UIMY_S91_Q_Dp.ivar2.fasta (54 files). These were placed in the correct corresponding folders manually.
# NOTE ERROR: Lote12 folder had repeated items, new commands for fixing the nomenclature are included @ /share/DiskCoViGen/VAMIP/05_variants/03_LANGEBIO/fix_Lote12_repeated_QZ.sh

# 02_INER
# A set of renaming tables (xref id to IMSS folio) were added to /share/DiskCoViGen/VAMIP/05_variants/IdSampleRelations/02_INER
cd /share/DiskCoViGen/VAMIP/05_variants/02_INER
# The script was tuned to work with INER's IMSS files (I'm still not sure how to deal with Cardio and CIENI)
bash rename_template.sh IMSS-1_2021-02-24 IMSS-1_2021-02-24.tsv
bash rename_template.sh IMSS-2_2021-03-24 IMSS-2_2021-03-24.tsv
bash rename_template.sh IMSS-3_2021-04-21 IMSS-3_2021-04-21.tsv
bash rename_template.sh IMSS-4_2021-05-19 IMSS-4_2021-05-19.tsv
bash rename_template.sh IMSS-5_2021-07-14 IMSS-5_2021-07-14.tsv
bash rename_template.sh IMSS-6_2021-08-04 IMSS-6_2021-08-04.tsv
bash rename_template.sh IMSS-7_2021-09-01 IMSS-7_2021-09-01.tsv
bash rename_template.sh IMSS-8_2021-09-29 IMSS-8_2021-09-29.tsv
bash rename_template.sh IMSS-9_2021-10-27 IMSS-9_2021-10-27.tsv
bash rename_template.sh IMSS-10_2021-11-24 IMSS-10_2021-11-24.tsv
# All renaming schemes are found in bash files which are used to trace back original names and should be executed to carry out renaming

# And rename batches' names in directories
mv IMSS-1_2021-02-24 Lote02_INER_25Febrero2021
mv IMSS-2_2021-03-24 Lote05_INER_24Marzo2021
mv IMSS-3_2021-04-21 Lote08_INER_21Abril2021
mv IMSS-4_2021-05-19 Lote11_INER_19Mayo2021
mv IMSS-5_2021-07-14 Lote17_INER_14Julio2021
mv IMSS-6_2021-08-04 Lote19_INER_04Agosto2021
mv IMSS-7_2021-09-01 Lote22_INER_01Septiembre2021
mv IMSS-8_2021-09-29 Lote25_INER_29Septiembre2021
mv IMSS-9_2021-10-27 Lote28_INER_27Octubre2021
mv IMSS-10_2021-11-24 Lote31_INER_24Noviembre2021

# # # # # # 04_CIAD DEPRECATED
# # # # # # All CIAD files were stored in the same folder as their nomenclature differs from CoViGen's. The renaming tables (dictionaries) were updated manually and stored @ /share/DiskCoViGen/VAMIP/05_variants/IdSampleRelations/04_CIAD/CIAD_2022-10-21.tsv
# # # # # # File /share/DiskCoViGen/VAMIP/05_variants/04_CIAD/rename_template.sh has the necessary commands to process a two table dictionary with old name to new name. Use as follows:
# # # # # bash rename_template.sh All_CIAD-2022-04-06 CIAD_2022-10-21.tsv
# # # # # # And run the resulting renaming commands (errors are expected whenever the old and new name are the same)
# # # # # bash rename_samples_All_CIAD-2022-04-06.sh
# # # # # # Some files were missing from the dictionary, thus, we change them manually (these belong to sample 202201038655 and its resulting files)
# # # # # for i in $(find All_CIAD-2022-04-06 -name "*test*"); do echo 'mv '$i' '$(echo $i|sed "s/_test//")'';done >remove_test_flag_CIAD-2022-04-06.sh

# 04_CIAD
# CIAD files are now complete and classified like this: FASTQ_Generation_2021-02-25_17_26_35Z-384253871
# I tracked all existing items with their upload (GISAID) names and put them in /share/DiskCoViGen/VAMIP/05_variants/IdSampleRelations/04_CIAD
# The next step will be to use these tables to create renamining commands
# This can be done with file /share/DiskCoViGen/VAMIP/05_variants/04_CIAD/rename_template.sh
# IMPORTANT NOTE: Many items cannot be identified, either because they were not uploaded or no metadata was available. These will result in renaming error as the old and new files will be the same
# IMPORTANT NOTE 2: No metadata was provided for items previous to CIAD11, which correspond to:
# FASTQ_Generation_2021-02-25_17_26_35Z-384253871
# FASTQ_Generation_2021-03-20_21_41_14Z-392576184
# FASTQ_Generation_2021-03-29_18_05_53Z-395935571
# FASTQ_Generation_2021-04-13_19_45_54Z-401452051
# FASTQ_Generation_2021-04-27_19_13_59Z-408745337
# FASTQ_Generation_2021-05-09_16_01_48Z-413894481
# FASTQ_Generation_2021-05-22_01_15_32Z-418698281
mkdir /share/DiskCoViGen/VAMIP/NON_USED/CIAD_no_metadata
mv FASTQ_Generation_2021-02-25_17_26_35Z-384253871 ../../NON_USED/CIAD_no_metadata/
mv FASTQ_Generation_2021-03-20_21_41_14Z-392576184 ../../NON_USED/CIAD_no_metadata/
mv FASTQ_Generation_2021-03-29_18_05_53Z-395935571 ../../NON_USED/CIAD_no_metadata/
mv FASTQ_Generation_2021-04-13_19_45_54Z-401452051 ../../NON_USED/CIAD_no_metadata/
mv FASTQ_Generation_2021-04-27_19_13_59Z-408745337 ../../NON_USED/CIAD_no_metadata/
mv FASTQ_Generation_2021-05-09_16_01_48Z-413894481 ../../NON_USED/CIAD_no_metadata/
mv FASTQ_Generation_2021-05-22_01_15_32Z-418698281 ../../NON_USED/CIAD_no_metadata/
mv FASTQ_Generation_2021-06-10_01_12_15Z-426349924 ../../NON_USED/CIAD_no_metadata/


# IMPORTANT NOTE 3: The ones that are actually renamed will be found with prefix CIADXX_ or LXXX_
# Create renaming schemes for each. Identation is used to identify how each folder was processed
bash rename_template.sh FASTQ_Generation_2021-06-10_01_12_15Z-426349924 CIAD11.tsv
	bash rename_samples_FASTQ_Generation_2021-06-10_01_12_15Z-426349924_CIAD11.sh
	mv FASTQ_Generation_2021-06-10_01_12_15Z-426349924 CIAD011_2021-06-10
bash rename_template.sh FASTQ_Generation_2021-06-23_17_54_26Z-431620203 CIAD12.tsv
	bash rename_samples_FASTQ_Generation_2021-06-23_17_54_26Z-431620203_CIAD12.sh
	mv FASTQ_Generation_2021-06-23_17_54_26Z-431620203 CIAD012_2021-06-23
bash rename_template.sh FASTQ_Generation_2021-07-08_15_34_29Z-436964529 CIAD13.tsv
	bash rename_samples_FASTQ_Generation_2021-07-08_15_34_29Z-436964529_CIAD13.sh
	mv FASTQ_Generation_2021-07-08_15_34_29Z-436964529 CIAD013_2021-07-08
bash rename_template.sh FASTQ_Generation_2021-07-10_14_58_12Z-437595159 CIAD14.tsv
	bash rename_samples_FASTQ_Generation_2021-07-10_14_58_12Z-437595159_CIAD14.sh
	mv FASTQ_Generation_2021-07-10_14_58_12Z-437595159 CIAD014_2021-07-10
bash rename_template.sh FASTQ_Generation_2021-07-14_22_49_42Z-438740309 CIAD15.tsv
	bash rename_samples_FASTQ_Generation_2021-07-14_22_49_42Z-438740309_CIAD15.sh
	mv FASTQ_Generation_2021-07-14_22_49_42Z-438740309 CIAD015_2021-07-14
bash rename_template.sh FASTQ_Generation_2021-07-27_22_13_50Z-443535092 CIAD16.tsv
	bash rename_samples_FASTQ_Generation_2021-07-27_22_13_50Z-443535092_CIAD16.sh
	mv FASTQ_Generation_2021-07-27_22_13_50Z-443535092 CIAD016_2021-07-27
bash rename_template.sh FASTQ_Generation_2021-08-03_16_46_40Z-445748303 CIAD17.tsv
	bash rename_samples_FASTQ_Generation_2021-08-03_16_46_40Z-445748303_CIAD17.sh
	mv FASTQ_Generation_2021-08-03_16_46_40Z-445748303 CIAD017_2021-08-03
bash rename_template.sh FASTQ_Generation_2021-08-25_22_22_40Z-453329878 CIAD18.tsv
	bash rename_samples_FASTQ_Generation_2021-08-25_22_22_40Z-453329878_CIAD18.sh
	mv FASTQ_Generation_2021-08-25_22_22_40Z-453329878 CIAD018_2021-08-25
bash rename_template.sh FASTQ_Generation_2021-09-17_16_55_31Z-461391932 CIAD19.tsv
	bash rename_samples_FASTQ_Generation_2021-09-17_16_55_31Z-461391932_CIAD19.sh
	mv FASTQ_Generation_2021-09-17_16_55_31Z-461391932 CIAD019_2021-09-17
bash rename_template.sh FASTQ_Generation_2021-09-28_18_05_38Z-465679223 CIAD20.tsv
	bash rename_samples_FASTQ_Generation_2021-09-28_18_05_38Z-465679223_CIAD20.sh
	mv FASTQ_Generation_2021-09-28_18_05_38Z-465679223 CIAD020_2021-09-28
bash rename_template.sh FASTQ_Generation_2021-10-01_20_13_44Z-467307842 CIAD21.tsv
bash rename_template.sh FASTQ_Generation_2021-10-01_20_13_44Z-467307842 CIAD22.tsv
bash rename_template.sh FASTQ_Generation_2021-10-01_20_13_44Z-467307842 CIAD23.tsv
	grep CIAD rename_samples_FASTQ_Generation_2021-10-01_20_13_44Z-467307842*|cut -d ":" -f 2|bash
	mv FASTQ_Generation_2021-10-01_20_13_44Z-467307842 CIAD021_022_023
bash rename_template.sh FASTQ_Generation_2021-10-08_17_35_11Z-470315852 CIAD24.tsv
	bash rename_samples_FASTQ_Generation_2021-10-08_17_35_11Z-470315852_CIAD24.sh
	mv FASTQ_Generation_2021-10-08_17_35_11Z-470315852 CIAD024_2021-10-08
bash rename_template.sh FASTQ_Generation_2021-11-09_20_23_59Z-485069591 CIAD25.tsv
	bash rename_samples_FASTQ_Generation_2021-11-09_20_23_59Z-485069591_CIAD25.sh
	mv FASTQ_Generation_2021-11-09_20_23_59Z-485069591 CIAD025_2021-11-09
bash rename_template.sh FASTQ_Generation_2021-11-30_17_35_43Z-494345853 CIAD26.tsv
	bash rename_samples_FASTQ_Generation_2021-11-30_17_35_43Z-494345853_CIAD26.sh
	mv FASTQ_Generation_2021-11-30_17_35_43Z-494345853 CIAD026_2021-11-30
bash rename_template.sh FASTQ_Generation_2021-12-10_21_29_00Z-499261764 CIAD27.tsv
	bash rename_samples_FASTQ_Generation_2021-12-10_21_29_00Z-499261764_CIAD27.sh
	mv FASTQ_Generation_2021-12-10_21_29_00Z-499261764 CIAD027_2021-12-10
bash rename_template.sh FASTQ_Generation_2021-12-24_18_04_57Z-505998494 CIAD28.tsv
bash rename_template.sh FASTQ_Generation_2021-12-24_18_04_57Z-505998494 CIAD29.tsv
	grep CIAD2 rename_samples_FASTQ_Generation_2021-12-24_18_04_57Z-505998494*|cut -d ":" -f 2|bash
	mv FASTQ_Generation_2021-12-24_18_04_57Z-505998494 CIAD028-029
bash rename_template.sh FASTQ_Generation_2022-01-08_20_09_05Z-510578068 CIAD30.tsv
	bash rename_samples_FASTQ_Generation_2022-01-08_20_09_05Z-510578068_CIAD30.sh
	mv FASTQ_Generation_2022-01-08_20_09_05Z-510578068 CIAD030_2022-01-08
bash rename_template.sh FASTQ_Generation_2022-01-18_14_58_58Z-514906392 CIAD31.tsv
bash rename_template.sh FASTQ_Generation_2022-01-18_14_58_58Z-514906392 CIAD32.tsv
	grep CIAD3 rename_samples_FASTQ_Generation_2022-01-18_14_58_58Z-514906392*|cut -d ":" -f 2|bash
	mv FASTQ_Generation_2022-01-18_14_58_58Z-514906392 CIAD031-032
bash rename_template.sh FASTQ_Generation_2022-02-04_16_28_02Z-523879361 CIAD33.tsv
bash rename_template.sh FASTQ_Generation_2022-02-04_16_28_02Z-523879361 CIAD36.tsv
	grep CIAD3 rename_samples_FASTQ_Generation_2022-02-04_16_28_02Z-523879361*|cut -d ":" -f 2|bash
	mv FASTQ_Generation_2022-02-04_16_28_02Z-523879361 CIAD033_036
# bash rename_template.sh FASTQ_Generation_2022-02-17_20_47_47Z-529765246 CIAD34_Lote42_09Febrero2022_CIAD.tsv # DEPRECATED
# 	bash rename_samples_FASTQ_Generation_2022-02-17_20_47_47Z-529765246_CIAD34_Lote42_09Febrero2022_CIAD.sh
# 	mv FASTQ_Generation_2022-02-17_20_47_47Z-529765246 Lote42_09Febrero2022_CIAD34.tsv
cd FASTQ_Generation_2022-02-17_20_47_47Z-529765246 # Instead, remane them separately as there are only IMSS samples here
	ls *.gz|sed -e 's/\(.*\)/mv \1 L042_\1/'|bash
	cd 03_assembly
	ls *|sed -e 's/\(.*\)/mv \1 L042_\1/'|bash
	cd ../..
	mv FASTQ_Generation_2022-02-17_20_47_47Z-529765246 Lote42_09Febrero2022_CIAD34.tsv
bash rename_template.sh FASTQ_Generation_2022-02-24_18_35_50Z-533226695 CIAD35.tsv
	bash rename_samples_FASTQ_Generation_2022-02-24_18_35_50Z-533226695_CIAD35.sh
	mv FASTQ_Generation_2022-02-24_18_35_50Z-533226695 CIAD035_2022-02-24
bash rename_template.sh FASTQ_Generation_2022-03-10_20_15_54Z-538832318 CIAD36.tsv # This one is more complicated as we have to first deal with non-IMSS samples
	bash rename_samples_FASTQ_Generation_2022-03-10_20_15_54Z-538832318_CIAD36.sh
bash rename_template.sh FASTQ_Generation_2022-03-10_20_15_54Z-538832318 CIAD36-37_Lote44_23Febrero2022_CIAD_INER.tsv
	bash rename_samples_FASTQ_Generation_2022-03-10_20_15_54Z-538832318_CIAD36-37_Lote44_23Febrero2022_CIAD_INER.sh
	mv FASTQ_Generation_2022-03-10_20_15_54Z-538832318 CIAD036_2022-03-10.tsv
cd FASTQ_Generation_2022-03-15_19_13_57Z-540710170
	ls *.gz|sed -e 's/\(.*\)/mv \1 L044_\1/'|bash
	cd 03_assembly
	ls *|sed -e 's/\(.*\)/mv \1 L044_\1/'|bash
	cd ../..
	mv FASTQ_Generation_2022-03-15_19_13_57Z-540710170 Lote44_23Febrero2022_CIAD036-037_INER # This will held all samples from this set, which are in two different folders
	# Now move the remaining L044 items to the corresponding folder
	mv CIAD036_2022-03-10.tsv/L044* Lote44_23Febrero2022_CIAD036-037_INER/
	mv CIAD036_2022-03-10.tsv/03_assembly/L044* Lote44_23Febrero2022_CIAD036-037_INER/03_assembly/
bash rename_template.sh FASTQ_Generation_2022-04-06_15_42_58Z-549015469 CIAD38.tsv # This one is more complicated as we have IMSS and non-IMSS samples. We fist deal with the non-IMSS samples
	bash rename_samples_FASTQ_Generation_2022-04-06_15_42_58Z-549015469_CIAD38.sh
	# The rest is renamed directly where 202x is found
	cd FASTQ_Generation_2022-04-06_15_42_58Z-549015469
	ls 202*.gz|sed -e 's/\(.*\)/mv \1 L048_\1/'|bash
	cd 03_assembly
	ls 202*|sed -e 's/\(.*\)/mv \1 L048_\1/'|bash
	cd ../..
	mv FASTQ_Generation_2022-04-06_15_42_58Z-549015469 Lote48_30Marzo2022_CIAD038_INER
cd FASTQ_Generation_2022-04-20_15_39_24Z-553606060
	ls 202*.gz|sed -e 's/\(.*\)/mv \1 L049_\1/'|bash
	cd 03_assembly
	ls 202*|sed -e 's/\(.*\)/mv \1 L049_\1/'|bash
	cd ../..
	mv FASTQ_Generation_2022-04-20_15_39_24Z-553606060 Lote49_13Abril2022_CIAD039_INER_IBT
cd FASTQ_Generation_2022-04-20_15_39_24Z-553606060
	ls 202*.gz|sed -e 's/\(.*\)/mv \1 L051_\1/'|bash
	cd 03_assembly
	ls 202*|sed -e 's/\(.*\)/mv \1 L051_\1/'|bash
	cd ../..
	mv FASTQ_Generation_2022-05-04_18_25_47Z-558713161 Lote51_27Abril2022_CIAD040
cd FASTQ_Generation_2022-05-26_15_28_41Z-566997431
	ls 202*.gz|sed -e 's/\(.*\)/mv \1 L053_\1/'|bash
	cd 03_assembly
	ls 202*|sed -e 's/\(.*\)/mv \1 L053_\1/'|bash
	cd ../..
	mv FASTQ_Generation_2022-05-26_15_28_41Z-566997431 Lote53_18Mayo2022_CIAD041
cd FASTQ_Generation_2022-06-22_22_44_22Z-576473905
	ls 202*.gz|sed -e 's/\(.*\)/mv \1 L056_\1/'|bash
	cd 03_assembly
	ls 202*|sed -e 's/\(.*\)/mv \1 L056_\1/'|bash
	cd ../..
	mv FASTQ_Generation_2022-06-22_22_44_22Z-576473905 Lote56_15Junio2022_CIAD043_p1
cd FASTQ_Generation_2022-06-30_17_50_49Z-579374797
	ls 202*.gz|sed -e 's/\(.*\)/mv \1 L056_\1/'|bash
	cd 03_assembly
	ls 202*|sed -e 's/\(.*\)/mv \1 L056_\1/'|bash
	cd ../..
	mv FASTQ_Generation_2022-06-30_17_50_49Z-579374797 Lote56_15Junio2022_CIAD043_p2
cd FASTQ_Generation_2022-07-08_19_30_46Z-582887308
	ls 202*.gz|sed -e 's/\(.*\)/mv \1 L059_\1/'|bash
	cd 03_assembly
	ls 202*|sed -e 's/\(.*\)/mv \1 L059_\1/'|bash
	cd ../..
	mv FASTQ_Generation_2022-07-08_19_30_46Z-582887308 Lote59_28Junio2022_CIAD045
bash rename_template.sh FASTQ_Generation_2022-07-20_18_49_25Z-587708124 CIAD47.tsv # This one includes both non-IMSS and IMSS items, we'll process them in that order
	bash rename_samples_FASTQ_Generation_2022-07-20_18_49_25Z-587708124_CIAD47.sh
	cd FASTQ_Generation_2022-07-20_18_49_25Z-587708124
	ls 202*.gz|sed -e 's/\(.*\)/mv \1 L061_\1/'|bash
	cd 03_assembly
	ls 202*|sed -e 's/\(.*\)/mv \1 L061_\1/'|bash
	cd ../..
	mv FASTQ_Generation_2022-07-20_18_49_25Z-587708124 Lote61_11Julio2022_CIAD046_plusCIAD047
bash rename_template.sh FASTQ_Generation_2022-08-18_22_34_43Z-598267673 CIAD48.tsv # This one includes both non-IMSS and IMSS items, we'll process them in that order
	bash rename_samples_FASTQ_Generation_2022-08-18_22_34_43Z-598267673_CIAD48.sh
	cd FASTQ_Generation_2022-08-18_22_34_43Z-598267673
	ls 202*.gz|sed -e 's/\(.*\)/mv \1 L065_\1/'|bash
	cd 03_assembly
	ls 202*|sed -e 's/\(.*\)/mv \1 L065_\1/'|bash
	cd ../..
	mv FASTQ_Generation_2022-08-18_22_34_43Z-598267673 Lote65_09Agosto2022_CIAD049_plusCIAD048
bash rename_template.sh FASTQ_Generation_2022-09-08_18_24_01Z-606250645 CIAD51.tsv # This one includes both non-IMSS and IMSS items, we'll process them in that order
	bash rename_samples_FASTQ_Generation_2022-09-08_18_24_01Z-606250645_CIAD51.sh
	cd FASTQ_Generation_2022-09-08_18_24_01Z-606250645
	ls 202*.gz|sed -e 's/\(.*\)/mv \1 L068_\1/'|bash
	cd 03_assembly
	ls 202*|sed -e 's/\(.*\)/mv \1 L068_\1/'|bash
	cd ../..
	mv FASTQ_Generation_2022-09-08_18_24_01Z-606250645 Lote68_30Agosto2022_CIAD050_plusCIAD051
cd FASTQ_Generation_2022-10-08_20_56_09Z-617138524
	ls 202*.gz|sed -e 's/\(.*\)/mv \1 L071_\1/'|bash
	cd 03_assembly
	ls 202*|sed -e 's/\(.*\)/mv \1 L071_\1/'|bash
	cd ../..
	mv FASTQ_Generation_2022-10-08_20_56_09Z-617138524 Lote71_27Septiembre2022_CIAD052
cd FASTQ_Generation_2022-11-25_22_12_22Z-634036406
	ls 202*.gz|sed -e 's/\(.*\)/mv \1 L074_\1/'|bash
	cd 03_assembly
	ls 202*|sed -e 's/\(.*\)/mv \1 L074_\1/'|bash
	cd ../..
	mv FASTQ_Generation_2022-11-25_22_12_22Z-634036406 Lote74_15Noviembre2022_CIAD053
mv FASTQ_Generation_2023-02-02_21_40_57Z-651299655/ Lote77_16Enero2023_CIAD56 # Files within this folder had already been renamed by CIAD members.

# RENAME FILES TO ADD BATCH NUMBER
cd /share/DiskCoViGen/VAMIP/05_variants
# P1000, Cardio and CIENI were transfered to folder: /share/DiskCoViGen/VAMIP/NON_USED
mkdir /share/DiskCoViGen/VAMIP/NON_USED
mv /share/DiskCoViGen/VAMIP/05_variants/01_IBT/P1000_1 /share/DiskCoViGen/VAMIP/NON_USED
mv /share/DiskCoViGen/VAMIP/05_variants/01_IBT/P1000_2 /share/DiskCoViGen/VAMIP/NON_USED
mv /share/DiskCoViGen/VAMIP/05_variants/02_INER/Cardio /share/DiskCoViGen/VAMIP/NON_USED
mv /share/DiskCoViGen/VAMIP/05_variants/02_INER/CIENI /share/DiskCoViGen/VAMIP/NON_USED
mv /share/DiskCoViGen/VAMIP/05_variants/03_LANGEBIO/Noviembre25 /share/DiskCoViGen/VAMIP/NON_USED
# UPDATE 2023-03-14: We have now processed these (except Noviembre25, which was deleted) and they are now found at 05_variants/00_Extra

# Gzip all fastas and tables
cd /share/DiskCoViGen/VAMIP/05_variants
find . -name "*.fasta" -exec gzip "{}" \;
find . -name "*_Q_Dp.ivar.tsv" -exec gzip "{}" \;

# 2022-04-27: We decided to add an additional identifier to the filenames like L001 for batch1, and so on
for i in $(find .|grep "fastq\|_Q_Dp");do echo 'mv '$i' '$(echo $i|sed "s/Lote\(\w\w\)\(.*\)\//Lote\1\2\/L0\1_/")'';done >rename_add_batches.sh

# Summarize contents
# Make sure the string "test" is not present in the file names (this means names have not been updated)
for i in $(find . -name "*test*"); do echo 'mv '$i' '$(echo $i|sed "s/_test//")'';done >2023-03-15_rm_test_flag.sh
# wc -l 2023-03-15_rm_test_flag.sh
# 56445 2023-03-15_rm_test_flag.sh
# We now need a list of all files having an ivar table
# for i in $(cat 2022-02-25.list);do base=$(basename $i);echo ''${base%_Q_Dp.ivar2.fasta}'';done >2022-02-25_samples_with_ivar.list
cd /share/DiskCoViGen/VAMIP/05_variants
# This will exclude all items whose names have not been processed (those flagged with string "test")
find ~+ -name "*.ivar.tsv.gz"|grep -v test >2023-03-14_ivar_tables.list # The ~+ expands the current path to full path
mkdir /share/DiskCoViGen/VAMIP/resp_2023-03-14_ivar_tables
cat 2023-03-14_ivar_tables.list| grep CIAD >2023-03-14_CIAD_only.tsv # This is latter used for a filter
# We now need the nucleotide contents
# File 2023-03-14_CIAD_wFolio.tsv has those new id folios. With this, we can remove the unwanted CIAD ones DEPRECATED
# for i in $(cat 2023-03-14_CIAD_wFolio.tsv);do echo 'grep -m 1 "03_assembly/'$i'_S" 2023-03-14_CIAD_only.tsv';done|bash|grep -vwFf - 2023-03-14_CIAD_only.tsv >2023-03-14_ignore_CIAD_files.tsv

cd /share/DiskCoViGen/VAMIP/05_variants
find ~+ -name "*.ivar2.fasta.gz"|grep -v test >2023-03-14_fastas.list # List all target consensus (ignore those flagged tests)
# grep "04_CIAD\|/202" 2023-03-14_ivar.list >2023-03-14_ivar.list.fix # Keep only those identified with a folio id or CIAD id
# Using the list of files, we get which ones were the original (raw, just after joinning lanes) files
find ~+ -name "*fastq.gz"|grep -v "01_Preprocess\|02_derep\|03_assembly" >2023-03-14_raw_fastqs.list
# We can now extract out of the original files, how many sequences they had
grep -v "/L0" 2023-03-14_raw_fastqs.list|grep -v 00_Extra >2023-03-14_Names_wo_batchIDs.txt # Get those without batch identifier
cut -d "/" -f 6-7 2023-03-14_Names_wo_batchIDs.txt|sort|uniq|grep Lote # Now print the folders so that samples can be renamed

# 03_LANGEBIO/Lote60_LANGEBIO_29Junio2022 # This one has a problem with files in 03_assembly, which were not renamed NOTE

# 01_IBT/Lote55_08Junio2020_IBT
# 01_IBT/Lote58_23Junio2022_IBT
# 01_IBT/Lote62_19Julio2022_IBT
# 01_IBT/Lote64_02Agosto2022_IBT
# 01_IBT/Lote67_23Agosto2022_IBT
# 01_IBT/Lote70_20Septiembre2022_IBT
# 01_IBT/Lote73_25Octubre2022_IBT
# 01_IBT/Lote76_12Enero2023_IBT
# 01_IBT/Lote79_30Enero2023_IBT

# 03_LANGEBIO/Lote50_LANGEBIO_20Abril2022 # These were renamed but with two digit names (e.g. L50 instead of L050)
# 03_LANGEBIO/Lote54_LANGEBIO_25Mayo2022
# 04_CIAD/Lote77_16Enero2023_CIAD56

# 03_LANGEBIO/Lote63_26Julio2022_LANGEBIO # Pending renaming NOTE
# 03_LANGEBIO/Lote69_LANGEBIO_06Septiembre2022

# 04_CIAD/Lote48_30Marzo2022_CIAD038_INER # These ones are not bad but have CIAD files (not IMSS) in them
# 04_CIAD/Lote61_11Julio2022_CIAD046_plusCIAD047
# 04_CIAD/Lote65_09Agosto2022_CIAD049_plusCIAD048
# 04_CIAD/Lote68_30Agosto2022_CIAD050_plusCIAD051

# All these problems were addressed according to each case's requirements

# UPDATE 2023-03-14: I created a script for renaming files within batches with their corresponding batch number:
bash rename_batches.sh 01_IBT/Lote55_08Junio2020_IBT L055 # The output was send to rename_cmd.sh
bash rename_batches.sh 01_IBT/Lote58_23Junio2022_IBT L058
bash rename_batches.sh 01_IBT/Lote62_19Julio2022_IBT L062
bash rename_batches.sh 01_IBT/Lote64_02Agosto2022_IBT L064
bash rename_batches.sh 01_IBT/Lote67_23Agosto2022_IBT L067
bash rename_batches.sh 01_IBT/Lote70_20Septiembre2022_IBT L070
bash rename_batches.sh 01_IBT/Lote73_25Octubre2022_IBT L073
bash rename_batches.sh 01_IBT/Lote76_12Enero2023_IBT L076
bash rename_batches.sh 01_IBT/Lote79_30Enero2023_IBT L079

# in the end, apparently only the following have any non-renamed items (although this is no error and we can proceed now):
# 03_LANGEBIO/Lote60_LANGEBIO_29Junio2022
# 04_CIAD/Lote48_30Marzo2022_CIAD038_INER
# 04_CIAD/Lote61_11Julio2022_CIAD046_plusCIAD047
# 04_CIAD/Lote65_09Agosto2022_CIAD049_plusCIAD048
# 04_CIAD/Lote68_30Agosto2022_CIAD050_plusCIAD051

for i in $(cat 2023-03-14_raw_fastqs.list);do echo ls -lh $i;done|bash >2023-03-14_rawfastq_size.tsv
# This was carried fo R1 and R2 because there are very rare cases where they don't match in total seqs
grep _R1.fastq.gz 2023-03-14_rawfastq_size.tsv >2023-03-14_sample_sizeR1.tsv
# For the raw files, we'll create a catalog with batch, Institute and size. Id is the folio IMSS ID (e.g. 202201000977) which is appendend as the las column (8)
paste <(cut -f 2 2023-03-14_sample_sizeR1.tsv) <(sed -e 's/\t.*academico /\t/' -e 's/\//\t/g' -e 's/ /\t/' -e 's/\t\(202.*\)\(_S.*\)/\t\1\2\t\1/' 2023-03-14_sample_sizeR1.tsv) >2023-03-14_sample_ok.tsv
# DEPRECATED: START
# # # # # Now, we discard those not having a valid IMSS ID
# grep "2021/2021\|2022/2022\|2022/2021\|P2/202\|2022/2.\|CIAD" 2023-03-14_sample_ok.tsv|sort >2023-03-14_sample_ok_target.list.tsv
# # # # # And keep track of those that didn't
# grep -v "2021/2021\|2022/2022\|2022/2021\|P2/202\|2022/2\|CIAD." 2023-03-14_sample_ok.tsv|sort >2023-03-14_sample_ok_target.bad.tsv
# DEPRECATED: END
cp 2023-03-14_sample_ok.tsv 2023-03-14_sample_ok_target.list.tsv

# The following commands should be optimized to run in parallel (maybe later)
#IMPORTANT: divisions by 0 fail and are seen as empty rows with only the filename
for i in $(cat 2023-03-14_fastas.list);do echo 'printf "'$i'\t" >>2023-03-14_nucleotide_contents.tsv;perl -ne '\''{next if $_=~ /^>/; chomp($_); $length=length($_); $Ns=tr/[N]//;$As=tr/[A]//; $Ts=tr/[T]//; $Cs=tr/[C]//; $Gs=tr/[G]//; $ATs=tr/[AT]//; $CGs=tr/[CG]//; $nonNs=tr/[ATCG]//; $Np=$Ns*100/$length; $nonNp=$nonNs*100/$length; $Ap=$As*100/$length; $Tp=$Ts*100/$length; $Cp=$Cs*100/$length; $Gp=$Gs*100/$length; $ATp=$ATs*100/$nonNs; $CGp=$CGs*100/$nonNs; $out=sprintf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f", $length, $nonNs, $Ns, $As, $Ts, $Cs, $Gs, $ATs, $CGs, $nonNp, $Np, $Ap, $Tp, $Cp, $Gp, $ATp, $CGp); print $out}'\'' <(zcat '$i') >>2023-03-14_nucleotide_contents.tsv; printf "\n" >>2023-03-14_nucleotide_contents.tsv';done >summarize_contents.sh
nohup bash summarize_contents.sh &


# Extract only those with non-empty items
awk -F'\t' '$2!=""' 2023-03-14_nucleotide_contents.tsv >2023-03-14_nucleotide_contents_nonEmtpy.tsv
# Now, append the IMSS folio as an extra column
paste 2023-03-14_nucleotide_contents_nonEmtpy.tsv <(sed -e "s/\//\t/g" 2023-03-14_nucleotide_contents_nonEmtpy.tsv -e 's/_S.*//' -e 's/_Q_Dp.*//'|cut -f 9) >2023-03-14_nucleotide_contents_wID.tsv
# sed -e "s/\//\t/g" 2023-03-14_nucleotide_contents_nonEmtpy.tsv|cut -f 9|sed -e 's/_Q_Dp.*//' -e 's/_S[0-9]*//'|grep -v "^L0"|less # check for missing batch info in the filenames



paste 2023-03-14_nucleotide_contents_nonEmtpy.tsv <(sed -e "s/\//\t/g" 2023-03-14_nucleotide_contents_nonEmtpy.tsv|cut -f 9|sed -e 's/_Q_Dp.*//' -e 's/_S[0-9]*//') >2023-03-14_nucleotide_contents_wID.tsv
# Additionally, we create a list of those having >25%, >15% and >10% Ns that are targeted for removal later on
awk -F'\t' '$12>=25' 2023-03-14_nucleotide_contents_wID.tsv|cut -f 1|sed -e 's/.*03_assembly\///' -e 's/\(_Q_Dp\).*//' >2023-03-14_for_removal_gt_25perc_Ns.tsv
awk -F'\t' '$12>=15' 2023-03-14_nucleotide_contents_wID.tsv|cut -f 1|sed -e 's/.*03_assembly\///' -e 's/\(_Q_Dp\).*//' >2023-03-14_for_removal_gt_15perc_Ns.tsv
awk -F'\t' '$12>=10' 2023-03-14_nucleotide_contents_wID.tsv|cut -f 1|sed -e 's/.*03_assembly\///' -e 's/\(_Q_Dp\).*//' >2023-03-14_for_removal_gt_10perc_Ns.tsv
awk -F'\t' '$12>=5' 2023-03-14_nucleotide_contents_wID.tsv|cut -f 1|sed -e 's/.*03_assembly\///' -e 's/\(_Q_Dp\).*//' >2023-03-14_for_removal_gt_05perc_Ns.tsv
for i in $(seq 0 100);do echo "awk -F'\t' '\$12>=$i' 2023-03-14_nucleotide_contents_wID.tsv|cut -f 1|sed -e 's/.*03_assembly\///' -e 's/\(_Q_Dp\).*//'|wc -l";done|bash >temp
printf "MaxNsCutoff\tBadGenomes\n" >2023-03-14_for_removal_gt_perc_Ns_ALL_percentages.tsv
paste <(seq 0 100) temp >>2023-03-14_for_removal_gt_perc_Ns_ALL_percentages.tsv

# Now, from the list of target samples, create a sublist of files to ignore (accept all from main IMSS batches [L0X], CIAD, INC, P1000 and RETRO"
cut -f 10 2023-03-14_sample_ok_target.list.tsv|grep -v "^L0\|^CIAD\|^P1000\|^RETRO\|^INC" >2023-03-14_IGNORED_samples.list
# wc -l 2023-03-14_IGNORED_samples.list
# 402 2023-03-14_IGNORED_samples.list
sed -e 's/_R1.fastq.*/_Q_Dp/' -e 's/^/\//' 2023-03-14_IGNORED_samples.list|grep -Ff - 2023-03-14_nucleotide_contents_wID.tsv >2023-03-14_IGNORED_samples_nucleotide_contents.tsv
# With this, we create a list of all available genomes that are part of the surveillance batches and other target samples (samples only, no controls) # IMPORTANT NOTE: This list contains the actual working samples for the downstream analyses
sed -e 's/_R1.fastq.*/_Q_Dp/' -e 's/^/\//' 2023-03-14_IGNORED_samples.list|grep -Fvf - 2023-03-14_nucleotide_contents_wID.tsv >2023-03-14_Main_mapped-Vigilancia.list
# wc -l 2023-03-14_Main_mapped-Vigilancia.list
# 32824 2023-03-14_Main_mapped-Vigilancia.list
cut -f 1 2023-03-14_Main_mapped-Vigilancia.list|sed 's/\.ivar2.*/.ivar.tsv.gz/' >2023-03-14_filtered_ivar_tables.list

# With both the Main_mapped_-Vigilancia.list and the 2023-03-14_sample_ok_target.list.tsv tables, we can now append size information to the nucleotide contents table using 2023-03-14_nucleotide_contents_wID.tsv as xref table
printf "File\tgenome_size\tNonNs\tNs\tAs\tTs\tCs\tGs\tATs\tCGs\tnonN%%\tN%%\tA%%\tT%%\tC%%\tG%%\tAT%%\tCG%%\tGenome\tPerm\tFile_size\n" >2023-03-14_nucleotide_contents_wSize.tsv
for i in $(cat 2023-03-14_Main_mapped-Vigilancia.list);do paste <(grep -m 1 "/$i\_S" 2023-03-14_nucleotide_contents_wID.tsv) <(grep -m 1 "/$i\_S" 2023-03-14_sample_ok_target.list.tsv|cut -f 2,3|cut -d " " -f 4);done >>2023-03-14_nucleotide_contents_wSize.tsv
# To decide upon this matter, we graph some an Rscript to graph some features


# ### Summarize ivar tables ###
# VERY IMPORTANT NOTE: we must first make sure all ivar tables are in gzip format for the following lines to work correctly
find . -name "*ivar.tsv" -exec gzip "{}" \;
cd /share/DiskCoViGen/VAMIP/05_variants
mkdir -p /scratch/rodrigog/01_Projects/VAMIP/06_VAMIP_tables/01_per_sample

# TEST: INIT
# Alt 1:
awk 'NR>=2 {new_var=">"$2"/"$3"/"$4"/"$11"/"$12"/"$5"/"$9" "; print new_var}' 202101077169_S38_Q_Dp.ivar.tsv|tr -d '\n';printf "\n"
# Example: >186/C/T/0.05/40/38/1
awk 'NR>=2 {new_var=">"$2"M"$3"2"$4"F"$11"T"$12"R"$5"A"$9" "; print new_var}' 202101077169_S38_Q_Dp.ivar.tsv|tr -d '\n';printf "\n"
# Example: >186MC2TF0.05T40R38A1
# # # # In order to summarize all tsv files in a folder do:
# # # for i in $(ls *.tsv);do echo 'printf "'$(echo $i|sed "s/_.*//")'\t";awk '\''NR>=2 {new_var=">"$2"/"$3"/"$4"/"$11"/"$12"/"$5"/"$9" "; print new_var}'\'' '$i'|tr -d "\n";printf "\n"';done
# # # for i in $(ls *.tsv);do echo 'printf "'$(echo $i|sed "s/_.*//")'\t";awk '\''NR>=2 {new_var=">"$2"M"$3"N"$4"F"$11"X"$12"R"$5"L"$9" "; print new_var}'\'' '$i'|tr -d "\n";printf "\n"';done
# # # # now,searches by mutation can be done by like this:
# # # egrep -o ">14933M\w*" Lote32.tsv
# TEST: END

# again, but this time, create a single file per ivar table:
# outdir="/scratch/rodrigog/01_Projects/VAMIP/06_VAMIP_tables/01_per_sample";for i in $(cat 2023-03-14_ivar_tables.list);do base=$(basename $i);echo 'printf "'$(echo $i|sed -e "s/.*05_variants\///" -e "s/03_assembly\///" -e "s/\//\\\\t/g" -e "s/_S[0-9]*_Q_Dp.ivar.tsv.gz//" -e "s/_Q_Dp_genome.ivar.tsv.gz//")'\t" >'$outdir'/'${base%.ivar.tsv.gz}'_vamip.tsv; awk '\''NR>=2 {new_var=">"$2"M"$3"N"$4"F"$11"X"$12"R"$5"L"$9" "; print new_var}'\'' <(zcat '$i')|tr . d|tr - m|tr + p|tr -d "\n" >>'$outdir/${base%.ivar.tsv.gz}'_vamip.tsv;printf "\n" >>'$outdir/${base%.ivar.tsv.gz}'_vamip.tsv';done >create_single_sample_vamip_tables.sh
outdir="/scratch/rodrigog/01_Projects/VAMIP/06_VAMIP_tables/01_per_sample";for i in $(cat 2023-03-14_filtered_ivar_tables.list);do base=$(basename $i);echo 'printf "'$(echo $i|sed -e "s/.*05_variants\///" -e "s/03_assembly\///" -e "s/\//\\\\t/g" -e "s/_Q_Dp.ivar.tsv.gz//" -e "s/_Q_Dp_genome.ivar.tsv.gz//")'\t" >'$outdir'/'${base%_Q_Dp.ivar.tsv.gz}'_vamip.tsv; awk '\''NR>=2 {new_var=">"$2"M"$3"N"$4"F"$11"X"$12"R"$5"L"$9" "; print new_var}'\'' <(zcat '$i')|tr . d|tr - m|tr + p|tr -d "\n" >>'$outdir/${base%_Q_Dp.ivar.tsv.gz}'_vamip.tsv;printf "\n" >>'$outdir/${base%_Q_Dp.ivar.tsv.gz}'_vamip.tsv';done >create_single_sample_vamip_tables.sh
# qsub script qsub_command_per_line.sh was created to run all lines in parallel (any bash script for one command per line)
mkdir -p $outdir/00_qsubout
qsub -R y -l h_rt=23:59:59 -l h_vmem=100M -N vamip_tables -o $outdir/00_qsubout -t 1-$(cat create_single_sample_vamip_tables.sh|wc -l) qsub_command_per_line.sh /share/DiskCoViGen/VAMIP/05_variants/create_single_sample_vamip_tables.sh
mv 01_per_sample/00_qsubout/ $outdir/.. # I decided to move the logs afterwards
# ls $outdir|grep "tsv$"|wc -l
# 32815
cd $outdir
rename "_Q_Dp_genome.ivar.tsv.gz_vamip" "vamip" * # Some items carry out an unwanted profix
cd /share/DiskCoViGen/VAMIP/05_variants
ls $outdir|grep "tsv$" >2022-03-28_All_current_VAMIP_coded_tables.tsv

# # IMPORTANT: FILTER UNWANTED ITEMS HERE

# DEPRECATED: START
# We use a list to filter vamip files
echo 482_S112_vamip.tsv >2023-03-14_bad_vamip.list

# sed -e 's/_Q_Dp.ivar.tsv.gz//' -e 's/.*\///' -e 's/$/_vamip.tsv/' 2023-03-14_ignore_CIAD_files.tsv >2023-03-14_bad_vamip.list
# ls 01_per_sample >files.list
# grep "_CN" files.list >>2023-03-14_bad_vamip.list
# grep "egativ" files.list >>2023-03-14_bad_vamip.list
# grep "EGATIV" files.list >>2023-03-14_bad_vamip.list
# grep "ILLUMINA" files.list >>2023-03-14_bad_vamip.list
# grep "_A_" files.list >>2023-03-14_bad_vamip.list
# grep "_B_" files.list >>2023-03-14_bad_vamip.list
# grep "_C_" files.list >>2023-03-14_bad_vamip.list
# grep "_D_" files.list >>2023-03-14_bad_vamip.list
# grep "_E_" files.list >>2023-03-14_bad_vamip.list
# grep "_F_" files.list >>2023-03-14_bad_vamip.list
# grep "_G_" files.list >>2023-03-14_bad_vamip.list
# grep "_H_" files.list >>2023-03-14_bad_vamip.list
# grep "_I_" files.list >>2023-03-14_bad_vamip.list
# grep "_J_" files.list >>2023-03-14_bad_vamip.list
# grep "_K_" files.list >>2023-03-14_bad_vamip.list
# grep "_L_" files.list >>2023-03-14_bad_vamip.list
# grep "^EX" files.list >>2023-03-14_bad_vamip.list
# grep "_INER-5-" files.list >>2023-03-14_bad_vamip.list
# grep "202201004764_S6_vamip.tsv" files.list >>2023-03-14_bad_vamip.list
# grep "_CP._" files.list >>2023-03-14_bad_vamip.list
# DEPRECATED: START

# Repeated files were also excluded
cd /scratch/rodrigog/01_Projects/VAMIP/06_VAMIP_tables
mkdir 01_per_sample-ignored 01_per_sample_gt_25perc_Ns 01_per_sample_gt_15perc_Ns 01_per_sample_gt_10perc_Ns
# Now, filter those that were selected
curdir="/share/DiskCoViGen/VAMIP/05_variants"
for i in $(cat $curdir/2023-03-14_bad_vamip.list);do echo 'mv 01_per_sample/'$i' 01_per_sample-ignored';done|bash
# And those with more than 25%, 15% and 10% Ns # I kept the errors on purpose to track which are no longer present
for i in $(cat $curdir/2023-03-14_for_removal_gt_25perc_Ns.tsv);do echo 'mv 01_per_sample/'$i'_vamip.tsv 01_per_sample_gt_25perc_Ns';done|bash 2>gt25_no_longer_present.tsv # Use the first list (>25%Ns to filter items with the lowest quality)
for i in $(cat $curdir/2023-03-14_for_removal_gt_25perc_Ns.tsv|grep -wvFf - $curdir/2023-03-14_for_removal_gt_15perc_Ns.tsv);do echo 'mv 01_per_sample/'$i'_vamip.tsv 01_per_sample_gt_15perc_Ns';done|bash 2>gt15_no_longer_present.tsv # use the remainder for the second tier (gt >15%Ns)
for i in $(cat $curdir/2023-03-14_for_removal_gt_15perc_Ns.tsv|grep -wvFf - $curdir/2023-03-14_for_removal_gt_10perc_Ns.tsv);do echo 'mv 01_per_sample/'$i'_vamip.tsv 01_per_sample_gt_10perc_Ns';done|bash 2>gt10_no_longer_present.tsv # Finally, do the same for those with >10%Ns
# IMPOTANT NOTE: Only genomes bearing <10% Ns were kept for dowstream analyses

# CIAD's files have not been renamed to include the batch, use the following to change this (this was done locally) DEPRECATED: START
# cd /home/rod/Documents/01_Projects/SARS/VAMIP
# # Fist get the list of CIAD files and keep only the column labeled "folio"
# grep "ciad\|Ciad\|CIAD" 2022_04_27_DB_MexCov2_v1.tsv|cut -f 13|sed '/^$/d' >2022_04_27_DB_MexCov2_v1-OnlyCIAD.tsv
# # Then, trim the L0XX part and keep only those that are not repeated (since we cannot know which ones match these.
# sed 's/L0[0-9][0-9]_//'  2022_04_27_DB_MexCov2_v1-OnlyCIAD.tsv|sort|uniq -c|tr -s " "|sed "s/ /\t/g"|awk -F'\t' '$2==1'|cut -f 3 >2022-04-27_unique_CIAD.txt
# # Now, get those of these that are present our actual set (as vamip files) and should be appended a batch name
# for i in $(cat 2022-04-27_unique_CIAD.txt);do echo 'printf "\n'$i'\t"; grep -m 1 "^'$i'_" <(ls 06_VAMIP_tables/01_per_sample)';done|bash|sed '/^$/d'|cut -f 2|sed '/^$/d' >2022-04-27_actual_CIAD_vamips.tsv
# # Do both searches and output to temporary file:
# for i in $(sed 's/_S[0-9]*_vamip.tsv//' 2022-04-27_actual_CIAD_vamips.tsv);do echo 'grep "^'$i'_S" 2022-04-27_actual_CIAD_vamips.tsv; grep "_'$i'$" 2022_04_27_DB_MexCov2_v1-OnlyCIAD.tsv';done|bash|paste - - >temp
# paste -d'\0' temp <(sed "s/.*\(S[0-9]*_vamip.tsv\).*/_\1/" temp)|sed -e 's/^/mv /' -e 's/\t/ /' >2022-04-27_add_Batches_to_CIAD.sh
# cd /home/rod/Documents/01_Projects/SARS/VAMIP/06_VAMIP_tables/01_per_sample # Locally
# DEPRECATED: START

# Now survey each position
cd /scratch/rodrigog/01_Projects/VAMIP/06_VAMIP_tables/
mkdir -p /scratch/rodrigog/01_Projects/VAMIP/06_VAMIP_tables/02_mut_by_pos
# for i in {1..30000};do echo 'egrep -m 1 -o ">'$i'M\w*" 01_per_sample/*|sed -e "s/.*01_per_sample\///" -e "s/_vamip.tsv:>/\t/" -e "s/_S[0-9]*\t/\t/" >02_mut_by_pos/pos'$i.tsv'';done >mutations_by_position.sh
# for i in {1..30000};do echo 'egrep -m 1 -o ">'$i'M\w*" 01_per_sample/*|sed -e "s/.*01_per_sample\///" -e ""s/_R2_Q/_Q/ -e "s/_Q_Dp_.*vamip.tsv:>/\t/" >02_mut_by_pos/pos'$i.tsv'';done >mutations_by_position.sh
for i in {1..30000};do echo 'grep -m 1 -o ">'$i'M\w*" 01_per_sample/*|sed -e "s/.*01_per_sample\///" -e "s/_R2_/_/" -e "s/_vamip.tsv:>/\t/" >02_mut_by_pos/pos'$i.tsv'';done >mutations_by_position.sh
# cd /share/DiskCoViGen/VAMIP/05_variants

# For use with a SGE queue system:
qsub -R y -l h_rt=23:59:59 -l h_vmem=100M -pe thread 1 -tc 60 -N SNP_positions -o /scratch/rodrigog/01_Projects/VAMIP/06_VAMIP_tables/00_qsubout -t 1-$(cat mutations_by_position.sh|wc -l) ~/bin/qsub_command_per_line.sh /scratch/rodrigog/01_Projects/VAMIP/06_VAMIP_tables/mutations_by_position.sh /scratch/rodrigog/01_Projects/VAMIP/06_VAMIP_tables/
# watch 'qacct -j 48684 | grep --no-group-separator maxvmem | cut -f7 -d" " | cut -f1 -d "."| sort -r -n | head -n 20'
# perl /usr/local/bin/calcul-eficiencia.pl -f 48684

# Alt for local run:
cd /home/rod/Documents/01_Projects/SARS/VAMIP/06_VAMIP_tables/
split -a 2 -d -l 5000 --additional-suffix .sh mutations_by_position.sh run_parallel
for i in run_parallel*;do echo nohup bash $i \&\> nohup_$i.out \&;done
# UPDATE 2022-04-29: The largest genome size in the current set is 29,909, so we delete non-extant items
for i in {29910..30000};do echo rm 02_mut_by_pos/pos$i.tsv;done|bash
# the refence genome, Wuhan-Hu-1 (https://www.ncbi.nlm.nih.gov/nuccore/1798174254) is 29,903 nt long but the last item having any mutation is pos29901.tsv, so the total genome length will be considered as the maximum genome size for percentage calculations.
wc -l 02_mut_by_pos/*|tr -s " "|sort -h|sed -e 's/^ //' -e 's/ /\t/' -e 's/02_mut_by_pos\/pos//' -e 's/.tsv//' >total_mutations_by_position.tsv


# Remove those with no mutations
find 02_mut_by_pos/ -type f -empty -delete
ls 02_mut_by_pos/|wc -l
# As of 2023-03-14, there were 28,684 (95.92%) positions that report any mutation
awk -F '\t' '$1==0' total_mutations_by_position.tsv |wc -l
# Only 1,225 (4.08) positions had no evidence of mutations
# These empty files are thus removed
mkdir 02_mut_by_pos/zeromut
for i in $(awk -F '\t' '{if($1==0) print $2}' total_mutations_by_position.tsv);do echo 'mv 02_mut_by_pos/pos'$i'.tsv 02_mut_by_pos/zeromut';done|bash
# Now remove those with only 1 observation (only one genome has them
mkdir 02_mut_by_pos/singletons
for i in $(wc -l 02_mut_by_pos/*.tsv|tr -s " "|sort -h|sed -e 's/^ //' -e 's/ /\t/'|awk '$1==1'|cut -f 2);do echo mv $i 02_mut_by_pos/singletons/;done|bash
ls 02_mut_by_pos/singletons/|wc -l
# 2,341 (7.83%) items were singletons
ls 02_mut_by_pos/*.tsv|wc -l
# There are 26,343 (88.09%) positions with more than 1 observation
awk -F '\t' '$1==2' total_mutations_by_position.tsv |wc -l
# from the remaining items, 2,881 (9.63%) had at least 2 observations (IMPORTANT NOTE: These will be the working set)
awk -F '\t' '$1>=5' total_mutations_by_position.tsv |wc -l
# and 18,083 (60.47%) had at least 5 observations
awk -F '\t' '$1<5' total_mutations_by_position.tsv |wc -l
# whereas 11,827 (39.55%) had less than 5 (and will not be considered for downstream analyses)
mkdir 02_mut_by_pos/items_2-4 # I further filtered it with to those not in at least 5 observations
for i in $(wc -l 02_mut_by_pos/*.tsv|tr -s " "|sort -h|sed -e 's/^ //' -e 's/ /\t/'|awk '$1<5'|cut -f 2);do echo mv $i 02_mut_by_pos/items_2-4;done|bash

# Now, create a global table and decode signs
cat 02_mut_by_pos/*.tsv >temp; paste <(cut -f 1 temp) <(cut -f 2 temp|sed -e 's/d/\./' -e 's/p/+/' -e 's/m/-/' -e 's/[MNFXRL]/\t/g') >All_mut_NoSingl.tsv;rm temp
wc -l All_mut_NoSingl.tsv
# 1,212,141 mutations will be considered for the current set

# This will be used to parse mutations with R (commands found locally @ /home/rod/Documents/01_Projects/SARS/VAMIP/explore_vamips.R)). Two extra tables are required:
# 1.- xref table for nt mutations to aa mutations @ UNIQUE_def_mutations_DeltaOmicron.tsv
# 2.- xref table for sample names to lineage and other tags @ 01_data/IMPORT_All_found_xref_with_names_fixed.tsv


# Get AA into ivar tables:
in_dir="/share/DiskCoViGen/VAMIP/05_variants/05_ivar_tables/bam"
out_dir="/share/DiskCoViGen/VAMIP/05_variants/05_ivar_tables/ivar"
mkdir -p $in_dir $out_dir
cd /share/DiskCoViGen/VAMIP/05_variants/05_ivar_tables/bam
for i in $(find /share/DiskCoViGen/VAMIP/05_variants/ -name "*.bam");do echo 'ln -s '$i' .';done|bash
ls >bams.list
cd ..
echo 'qsub -R y -l h_rt=23:59:59 -l h_vmem=1G -N ivar_AA -o /scratch/rodrigog/qsub/ -t 1-'$(ls bam/*.bam|wc -l)' /share/DiskCoViGen/VAMIP/Programs/ivar_variant_table_AA.sh '$in_dir' bams.list /share/DiskCoViGen/VAMIP/DB/MN908947.fasta /share/DiskCoViGen/VAMIP/DB/Wuhan-Hu-1_MN908947_fix.gff3 '$out_dir''

cd /run/media/rod/Bicho/Tables_ivar/ivar/
gzip *; cd ..
tar -cvf ivar.tar ivar/ # This was then copied to my local host
# from elsewhere (/run/media/rod/Bicho/Tables_ivar), we can now create a single table with everything
for i in $(ls ivar/*.gz);do base=$(basename $i);echo 'zcat '$i'|sed "s/^/'${base%_Q_Dp.ivar.tsv.gz}'\t/"';done|bash > single_ivar_table.tsv

#################################################################################################################################################   UPDATE   ######################################################### ###########################################################################################################
# IMPORTANT UPDATE 2022-11-24: Previously, we had all samples up to Lote54_LANGEBIO_25Mayo2022, with missing batches as follows: Lote48_30Marzo2022_CIAD_INER
# Lote48_30Marzo2022_CIAD_INER
# Lote49_13Abril2022_CIAD_INER_IBT
# Lote51_27Abril2022_CIAD_EpiCoV
# Lote53_18Mayo2022_CIAD_EpiCoV
# We added the missing fastq sets and the rest up until Lote73_25Octubre2022_IBT.

# Also, the IBt sets for the P1000 project, INER sets for Cardio and CENDI xref identified with their metadata. These were all in /share/DiskCoViGen/VAMIP/NON_USED/ and they were mostly from the first year (2020).
# First, process P1000_1. The xref file was deposited in the same folder as 1000Genomas_PrimerSampleSeqRelation.txt and 1000Genomas_PrimerMeta_v3.tsv
cd /share/DiskCoViGen/VAMIP/NON_USED/P1000_1
for i in $(ls *R1*.fastq.gz);do printf ''$(echo $i|sed "s/_.*//")'\t'$(echo $i|sed "s/_R1.*//")'\t'$i'\n';done >sample_list_p1000_1.tsv
paste <(cut -f 27 1000Genomas_PrimerMeta_v3.tsv) <(cut -f 3 1000Genomas_PrimerMeta_v3.tsv)|sed -e 's/SARS-CoV-2_//' -e 's/_/\t/' >sample_list_p1000_1.tsv
# The resulting file was manually curated into file sample_list_p1000_1_xref.tsv to append a code containing the new name that identifies each sample as belonging to the P1000 set and their laboratory
# Next, rename the files
grep -v "Illumina" sample_list_p1000_1_xref.tsv |paste <(cut -f 3,7)|sed -e 's/^/mv /' -e 's/$/_R1.fastq.gz/' -e 's/\t/ /'|grep -v " _R1.fastq.gz" >rename.sh
grep -v "Illumina" sample_list_p1000_1_xref.tsv |paste <(cut -f 3,7)|sed -e 's/^/mv /' -e 's/$/_R1.fastq.gz/' -e 's/\t/ /' -e 's/_R1/_R2/g'|grep -v " _R2.fastq.gz" >> rename.sh
bash rename.sh
# and the items in the folder 03_assembly
grep "_R1.*\.fastq.gz" rename.sh|sed -e 's/_R1.fastq.gz/_Q_Dp_test.ivar2.fasta.gz/g' -e 's/ / 03_assembly\//g' >rename2.sh
grep "_R1.*\.fastq.gz" rename.sh|sed -e 's/_R1.fastq.gz/_Q_Dp_test.ivar.qual.txt.gz/g' -e 's/ / 03_assembly\//g'>>rename2.sh
grep "_R1.*\.fastq.gz" rename.sh|sed -e 's/_R1.fastq.gz/_Q_Dp_test.ivar.tsv/g' -e 's/ / 03_assembly\//g'>>rename2.sh
grep "_R1.*\.fastq.gz" rename.sh|sed -e 's/_R1.fastq.gz/_Q_Dp_test.srt.bam/g' -e 's/ / 03_assembly\//g'>>rename2.sh
grep "_R1.*\.fastq.gz" rename.sh|sed -e 's/_R1.fastq.gz/_Q_Dp_test.srt.bam.bai/g' -e 's/ / 03_assembly\//g'>>rename2.sh
grep "_R1.*\.fastq.gz" rename.sh|sed -e 's/_R1.fastq.gz/_Q_Dp_test.vcf.gz/g' -e 's/ / 03_assembly\//g'>>rename2.sh
# and remove the "test" flag (used to track new items)
for i in $(find . -name "*test*"); do echo 'mv '$i' '$(echo $i|sed "s/_test//")'';done >remove_test_flag.sh

# Now, for the second part (P1000_2)
cd /share/DiskCoViGen/VAMIP/NON_USED/P1000_2
# First, eliminate those labeled RNA as they are not part of this study
for i in $(find . -name "RNA_*");do echo rm $i;done|bash
# Now, process them all with the same protocol, first, get a list of samples and their id within the Illumina flow-cell.
for i in $(ls *R1.fastq.gz);do printf ''$(echo $i|sed -e "s/RNA_//" -e "s/_.*//")'\t'$(echo $i|sed "s/_R1.*//")'\t'$i'\n';done >sample_list_p1000_2.tsv
paste <(cut -f 27 1000Genomas_SegundoMeta.tsv) <(cut -f 3 1000Genomas_PrimerMeta_v3.tsv)|sed -e 's/SARS-CoV-2_//' -e 's/_/\t/' >sample_list_p1000_2.tsv
# The resulting file was manually curated into file sample_list_p1000_1_xref.tsv to append a code containing the new name that identifies each sample as belonging to the P1000 set and their laboratory
# Next, remane items
grep -v "Illumina" sample_list_p1000_2_xref.tsv |paste <(cut -f 3,7)|sed -e 's/^/mv /' -e 's/$/_R1.fastq.gz/' -e 's/\t/ /'|grep -v " _R1.fastq.gz" > rename.sh
grep -v "Illumina" sample_list_p1000_2_xref.tsv |paste <(cut -f 3,7)|sed -e 's/^/mv /' -e 's/$/_R1.fastq.gz/' -e 's/\t/ /' -e 's/_R1.fastq.gz/_R2.fastq.gz/g'|grep -v " _R2.fastq.gz" >> rename.sh
bash rename.sh
# and the items in the folder 03_assembly
grep "_R1\.fastq.gz" rename.sh|sed -e 's/_R1.fastq.gz/_Q_Dp_test.ivar2.fasta.gz/g' -e 's/ / 03_assembly\//g' >rename2.sh
grep "_R1\.fastq.gz" rename.sh|sed -e 's/_R1.fastq.gz/_Q_Dp_test.ivar.qual.txt.gz/g' -e 's/ / 03_assembly\//g'>>rename2.sh
grep "_R1\.fastq.gz" rename.sh|sed -e 's/_R1.fastq.gz/_Q_Dp_test.ivar.tsv/g' -e 's/ / 03_assembly\//g'>>rename2.sh
grep "_R1\.fastq.gz" rename.sh|sed -e 's/_R1.fastq.gz/_Q_Dp_test.srt.bam/g' -e 's/ / 03_assembly\//g'>>rename2.sh
grep "_R1\.fastq.gz" rename.sh|sed -e 's/_R1.fastq.gz/_Q_Dp_test.srt.bam.bai/g' -e 's/ / 03_assembly\//g'>>rename2.sh
grep "_R1\.fastq.gz" rename.sh|sed -e 's/_R1.fastq.gz/_Q_Dp_test.vcf.gz/g' -e 's/ / 03_assembly\//g'>>rename2.sh
# and remove the "test" flag (used to track new items)
for i in $(find . -name "*test*"); do echo 'mv '$i' '$(echo $i|sed "s/_test//")'';done >remove_test_flag.sh


# Next, for the Cardio set
for i in $(ls *R1*);do name=$(echo $i|sed -e "s/.*Cardio-1-//" -e "s/-.*//" -e "s/^/hCoV-19\/Mexico\/CMX-INER-INC-/" -e "s/$/\/2021/");new_name=$(echo $i|sed -e "s/.*Cardio-1-//" -e "s/-.*//" -e "s/^/INC-/" -e "s/$/_INER/"); printf "$i\t$name\t$new_name\n";done >xref_INER_metadata.tsv
# These were all from 2021.
cut -f 1,3 xref_INER_metadata.tsv|sed -e 's/\t/ /' -e 's/^/mv /' -e 's/$/_R1.fastq.gz/' >rename_samples.sh
cut -f 1,3 xref_INER_metadata.tsv|sed -e 's/\t/ /' -e 's/^/mv /' -e 's/$/_R1.fastq.gz/' -e 's/_R1\.fastq/_R2\.fastq/g' >>rename_samples.sh


# Now, for the CIENI set (which are actually from INER and CIENI). These will be renamed with the prefix RETRO as they are from retrospective sets.
cd /home/rod/Documents/01_Projects/SARS/VAMIP/00_docs/
paste <(cut -f 1 2022-12-12_Celia_CIENI_xref.tsv) <(cut -f 2 2022-12-12_Celia_CIENI_xref.tsv|sed -e "s/CoV-\([A-Z]*\)-\([0-9]*\)/RETRO-\2_\1/") <(cut -f 2- 2022-12-12_Celia_CIENI_xref.tsv) >2023-01-16_Celia_xref.tsv
grep -v Nombre 2023-01-16_Celia_xref.tsv|cut -f 1-2|sed -e 's/^/mv /' -e 's/\t/ /' -e 's/$/_R1.fastq.gz/' >rename.sh
grep -v Nombre 2023-01-16_Celia_xref.tsv|cut -f 1-2|sed -e 's/^/mv /' -e 's/\t/ /' -e 's/$/_R1.fastq.gz/' -e 's/_R1\./_R2\./g' >>rename.sh
# Now, for the 03_assembly folder:
grep -v Nombre 2023-01-16_Celia_xref.tsv|cut -f 1-2|sed -e 's/_R1.fastq.gz\t/ /' -e 's/$/ /' -e 's/ /_Q_Dp_test.ivar2.fasta.gz /g' -e 's/ $//' -e 's/^/mv /' -e 's/ / 03_assembly\//g'




paste <(find 03_assembly/ -name "*Dp*") <(find 03_assembly/ -name "*Dp*") >temp
grep -v Nombre 2023-01-16_Celia_xref.tsv|cut -f 1,2|sed -e 's/_R1.fastq.gz\t/\//' -e 's/^/sed -i "s\//' -e 's/$/\/" temp/' |bash
paste <(cut -f 2 temp) <(cut -f 1 temp)|sed -e 's/^/mv /' -e 's/\t/ /' >rename2.sh
rm temp

