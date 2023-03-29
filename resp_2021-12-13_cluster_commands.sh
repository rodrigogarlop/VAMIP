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
out_dir="/scratch/rodrigog/01_Projects/VAMIP/04_fastq/temp/04_CIAD"
mkdir $out_dir cd $out_dir
in_dir="/share/DiskCoViGen/VAMIP/04_fastq/04_CIAD"
for i in $(find $in_dir -name "*.gz");do echo cp $i /scratch/rodrigog/01_Projects/VAMIP/04_fastq/temp/04_CIAD/;done >cp_files.sh
nohup bash cp_files.sh &
mv nohup.out FAILED_fastqs_missing_lanes.txt
rename ".R1" '_R1' *.gz # Names from CIAD files have a uncommon .R1 separator, which is not compatible with our current scripts
rename ".R2" '_R2' *.gz
mkdir -p 00_qsubout 01_Preprocess 02_derep 03_assembly
ls *.gz >files.list
# Filter reads with fastp
out_dir="/scratch/rodrigog/01_Projects/VAMIP/04_fastq/temp/04_CIAD"
echo 'qsub -R y -l h_rt=23:59:59 -N preProcess -o '$out_dir'/00_qsubout -t 1-808:2 /share/DiskCoViGen/VAMIP/Programs/preProcess_list_v2.sh '$out_dir' files.list '$out_dir'/01_Preprocess'|bash
# dereplicate with cdhit
for i in $(ls 01_Preprocess/*.fastq.gz);do basename $i;done >$out_dir/01_Preprocess/files_clean.list
echo 'qsub -R y -l h_rt=23:59:59 -N cdHitDup_CIAD -o '$out_dir'/00_qsubout -t 1-'$(ls 01_Preprocess/*.gz|wc -l)':2 /share/DiskCoViGen/VAMIP/Programs/cdHitDupPair_list_v2.1.sh '$out_dir/01_Preprocess' files_clean.list '$out_dir/02_derep''|bash
# # Test
# for i in $(ls 01_Preprocess/*.gz);do basename $i;done|sed 's/_R.*//'|sort|uniq >01_Preprocess_out.txt
# for i in $(ls 02_derep/*.gz);do basename $i;done|sed 's/_R.*//'|sort|uniq >02_derep_out.txt
# Get which ones fail
a=$(ls 01_Preprocess/*.gz|wc -l);grep -m 1 "Done\!" 00_qsubout/cdHitDup*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done)|grep -f - <(ls 00_qsubout/cdHitDup*) >00_qsubout/FAILED_cdHitDup.txt
# Assembly and analyze changes vs Original SARS-CoV-2 reference
for i in $(ls 02_derep/*.fastq.gz);do basename $i;done >$out_dir/02_derep/files_derep.list
echo 'qsub -R y -l h_rt=23:59:59 -N ConsesusVariant_CIAD -o '$out_dir'/00_qsubout -t 1-'$(ls 02_derep/*.gz|wc -l)':2 /share/DiskCoViGen/VAMIP/Programs/variantConsePairCall_List_v2.1.sh '$out_dir/02_derep' files_derep.list /share/DiskCoViGen/VAMIP/DB/MN908947.fasta genome '$out_dir/03_assembly''|bash
# Get which ones fail
a=$(ls 01_Preprocess/*.gz|wc -l);grep -m 1 "Done\!" 00_qsubout/ConsesusVariant*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done)|grep -f - <(ls 00_qsubout/ConsesusVariant*) >00_qsubout/FAILED_ConsensusVariant.txt

# LANGEBIO
month="Marzo1"
out_dir="/scratch/rodrigog/01_Projects/VAMIP/04_fastq/temp/03_LANGEBIO/$month"; mkdir -p $out_dir; cd $out_dir
# Concatenate 4 lanes into a single fastq
in_dir="/share/DiskCoViGen/VAMIP/04_fastq/03_LANGEBIO/$month/reads"
for sam in $(ls $in_dir|grep fastq.gz|sed "s/_L00.*//"|sort|uniq);do for pair in R1 R2;do printf "cat"; for lane in L001 L002 L003 L004; do printf ' '$in_dir/$sam'_'$lane'_'$pair'_001.fastq.gz ';done;printf ' >'$out_dir'/'$sam'_'$pair'.fastq.gz\n';done;done >cat_files.sh
nohup bash cat_files.sh &
mv nohup.out FAILED_fastqs_missing_lanes.txt
cd $out_dir; mkdir -p 00_qsubout 01_Preprocess 02_derep 03_assembly
ls *.gz|sort >files.list
# Filter reads with fastp
echo 'qsub -R y -l h_rt=23:59:59 -N preProcess_LANGEBIO_'$month' -o '$out_dir'/00_qsubout -t 1-'$(ls *.gz|wc -l)':2 /share/DiskCoViGen/VAMIP/Programs/preProcess_list_v2.sh '$out_dir' files.list '$out_dir'/01_Preprocess'|bash
ls 01_Preprocess/*.gz|wc -l
# dereplicate with cdhit
for i in $(ls 01_Preprocess/*.fastq.gz);do basename $i;done >$out_dir/01_Preprocess/files_clean.list
echo 'qsub -R y -l h_rt=23:59:59 -N cdHitDup_LANGEBIO_'$month' -o '$out_dir'/00_qsubout -t 1-'$(ls 01_Preprocess/*.gz|wc -l)':2 /share/DiskCoViGen/VAMIP/Programs/cdHitDupPair_list_v2.1.sh '$out_dir/01_Preprocess' files_clean.list '$out_dir/02_derep''|bash
# Get which ones fail
a=$(ls 01_Preprocess/*.gz|wc -l);grep -m 1 "Done\!" 00_qsubout/cdHitDup*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done)|grep -f - <(ls 00_qsubout/cdHitDup*) >00_qsubout/FAILED_cdHitDup.txt;cat 00_qsubout/FAILED_cdHitDup.txt
# Assembly and analyze changes vs Original SARS-CoV-2 reference
for i in $(ls 02_derep/*.fastq.gz);do basename $i;done >$out_dir/02_derep/files_derep.list
echo 'qsub -R y -l h_rt=23:59:59 -N ConsesusVariant_LANGEBIO_'$month' -o '$out_dir'/00_qsubout -t 1-'$(ls 02_derep/*.gz|wc -l)':2 /share/DiskCoViGen/VAMIP/Programs/variantConsePairCall_List_v2.1.sh '$out_dir/02_derep' files_derep.list /share/DiskCoViGen/VAMIP/DB/MN908947.fasta test '$out_dir/03_assembly''|bash
# Get which ones fail
a=$(ls 01_Preprocess/*.gz|wc -l);grep -m 1 "Done\!" 00_qsubout/ConsesusVariant*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done)|grep -f - <(ls 00_qsubout/ConsesusVariant*) >00_qsubout/FAILED_ConsensusVariant.txt; cat 00_qsubout/FAILED_ConsensusVariant.txt

# IBT
month="Lote01_IBT_17Febrero2021"
out_dir="/scratch/rodrigog/01_Projects/VAMIP/04_fastq/temp/01_IBT/$month"; mkdir -p $out_dir; cd $out_dir
# Concatenate 4 lanes into a single fastq
in_dir="/share/DiskCoViGen/VAMIP/04_fastq/01_IBT/$month"
for sam in $(ls $in_dir/*.gz);do echo 'cp '$sam' '$out_dir/'';done >cp_files.sh
nohup bash cp_files.sh &
mv nohup.out FAILED_fastqs_missing_lanes.txt
cd $out_dir; mkdir -p 00_qsubout 01_Preprocess 02_derep 03_assembly
ls *.gz >files.list
# Filter reads with fastp
echo 'qsub -R y -l h_rt=23:59:59 -N preProcess_IBT_'$month' -o '$out_dir'/00_qsubout -t 1-'$(ls *.gz|wc -l)':2 /share/DiskCoViGen/VAMIP/Programs/preProcess_list_v2.sh '$out_dir' files.list '$out_dir'/01_Preprocess'|bash
ls 01_Preprocess/*.gz|wc -l
# dereplicate with cdhit
for i in $(ls 01_Preprocess/*.fastq.gz);do basename $i;done >$out_dir/01_Preprocess/files_clean.list
echo 'qsub -R y -l h_rt=23:59:59 -N cdHitDup_IBT_'$month' -o '$out_dir'/00_qsubout -t 1-'$(ls 01_Preprocess/*.gz|wc -l)':2 /share/DiskCoViGen/VAMIP/Programs/cdHitDupPair_list_v2.1.sh '$out_dir/01_Preprocess' files_clean.list '$out_dir/02_derep''|bash
# Get which ones fail
a=$(ls 01_Preprocess/*.gz|wc -l);for x in $(grep -m 1 "Done\!" 00_qsubout/cdHitDup*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done));do echo ".$x";done|grep -f - <(ls 00_qsubout/cdHitDup*) >00_qsubout/FAILED_cdHitDup.txt;cat 00_qsubout/FAILED_cdHitDup.txt
# Assembly and analyze changes vs Original SARS-CoV-2 reference
for i in $(ls 02_derep/*.fastq.gz);do basename $i;done >$out_dir/02_derep/files_derep.list
echo 'qsub -R y -l h_rt=23:59:59 -N ConsesusVariant_IBT_'$month' -o '$out_dir'/00_qsubout -t 1-'$(ls 02_derep/*.gz|wc -l)':2 /share/DiskCoViGen/VAMIP/Programs/variantConsePairCall_List_v2.1.sh '$out_dir/02_derep' files_derep.list /share/DiskCoViGen/VAMIP/DB/MN908947.fasta test '$out_dir/03_assembly''|bash
# Get which ones fail
a=$(ls 02_derep/*.gz|wc -l);for x in $(grep -m 1 "Done\!" 00_qsubout/ConsesusVariant*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done));do echo ".$x";done|grep -f - <(ls 00_qsubout/ConsesusVariant*) >00_qsubout/FAILED_ConsensusVariant.txt; cat 00_qsubout/FAILED_ConsensusVariant.txt
ls 03_assembly/*.fasta|wc -l








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
# # # # Assembly and analyze changes vs Original SARS-CoV-2 reference
# # # echo 'qsub -R y -l h_rt=23:59:59 -N fix_ivar_LANGEBIO_'$month' -o '$out_dir'/00_qsubout -t 1-'$(cat 02_derep/files_derep.list|wc -l)':2 /share/DiskCoViGen/VAMIP/Programs/ivar_only_v1.sh '$out_dir/02_derep' files_derep.list /share/DiskCoViGen/VAMIP/DB/MN908947.fasta test '$out_dir/03_assembly''|bash
# # # # Get which ones fail
# # # a=$(cat 02_derep/files_derep.list|wc -l);for x in $(grep -m 1 "Done\!" 00_qsubout/fix_ivar*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done));do echo ".$x";done|grep -f - <(ls 00_qsubout/fix_ivar*) >00_qsubout/FAILED_Fixivar.txt; cat 00_qsubout/FAILED_Fixivar.txt
# # # rm 03_assembly/*.bam
# # # nohup cp 03_assembly/ /share/DiskCoViGen/VAMIP/05_variants/01_IBT/$month/03_assembly/ &


#### This is done directly in the output folders in DiskCoViGen:
month="septiembre"
ins="03_LANGEBIO"
out_dir="/share/DiskCoViGen/VAMIP/05_variants/$ins/$month"; cd $out_dir
rm 03_assembly/*.gz # remove bad tables
# Assembly and analyze changes vs Original SARS-CoV-2 reference
echo 'qsub -R y -l h_rt=23:59:59 -N fix_ivar_LANGEBIO_'$month' -o '$out_dir'/00_qsubout -t 1-'$(cat 02_derep/files_derep.list|wc -l)':2 /share/DiskCoViGen/VAMIP/Programs/ivar_only_v1.sh '$out_dir/02_derep' files_derep.list /share/DiskCoViGen/VAMIP/DB/MN908947.fasta test '$out_dir/03_assembly''|bash
# Get which ones fail
a=$(cat 02_derep/files_derep.list|wc -l);for x in $(grep -m 1 "Done\!" 00_qsubout/fix_ivar*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done));do echo ".$x";done|grep -f - <(ls 00_qsubout/fix_ivar*) >00_qsubout/FAILED_Fixivar.txt; cat 00_qsubout/FAILED_Fixivar.txt
ls 03_assembly/*.fasta|wc -l

# CIAD
out_dir="/share/DiskCoViGen/VAMIP/05_variants/04_CIAD"; cd $out_dir
rm 03_assembly/*.gz # remove bad tables
# Assembly and analyze changes vs Original SARS-CoV-2 reference
echo 'qsub -R y -l h_rt=23:59:59 -N fix_ivar_CIAD -o '$out_dir'/00_qsubout -t 1-'$(cat 02_derep/files_derep.list|wc -l)':2 /share/DiskCoViGen/VAMIP/Programs/ivar_only_v1.sh '$out_dir/02_derep' files_derep.list /share/DiskCoViGen/VAMIP/DB/MN908947.fasta test '$out_dir/03_assembly''|bash
# Get which ones fail
a=$(cat 02_derep/files_derep.list|wc -l);for x in $(grep -m 1 "Done\!" 00_qsubout/fix_ivar*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done));do echo ".$x";done|grep -f - <(ls 00_qsubout/fix_ivar*) >00_qsubout/FAILED_Fixivar.txt; cat 00_qsubout/FAILED_Fixivar.txt
ls 03_assembly/*.fasta|wc -l
