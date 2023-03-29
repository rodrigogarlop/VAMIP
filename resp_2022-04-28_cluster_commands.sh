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
out_dir="/scratch/rodrigog/01_Projects/VAMIP/04_fastq/temp/All_CIAD-2022-04-06/"
mkdir -p $out_dir; cd $out_dir
in_dir="/scratch/rodrigog/01_Projects/VAMIP/04_fastq/temp/04_CIAD"
for i in $(find $in_dir -name "*.gz");do echo cp $i $out_dir;done >cp_files.sh
nohup bash cp_files.sh &
mv nohup.out FAILED_fastqs_missing_lanes.txt
rename "_L001_" "_" *.gz
rename "_001.fastq" ".fastq" *.gz
# rename ".R1" '_R1' *.gz # Names from CIAD files have a uncommon .R1 separator, which is not compatible with our current scripts
# rename ".R2" '_R2' *.gz
mkdir -p 00_qsubout 01_Preprocess 02_derep 03_assembly
ls *.gz >files.list
# Filter reads with fastp
# out_dir="/scratch/rodrigog/01_Projects/VAMIP/04_fastq/temp/04_CIAD"
echo 'qsub -R y -l h_rt=23:59:59 -N preProcess_CIAD -o '$out_dir'/00_qsubout -t 1-'$(ls *.gz|wc -l)':2 /share/DiskCoViGen/VAMIP/Programs/preProcess_list_v2.sh '$out_dir' files.list '$out_dir'/01_Preprocess'|bash
### echo 'qsub -R y -l h_rt=23:59:59 -N preProcess -o '$out_dir'/00_qsubout -t 1-808:2 /share/DiskCoViGen/VAMIP/Programs/preProcess_list_v2.sh '$out_dir' files.list '$out_dir'/01_Preprocess'|bash
# dereplicate with cdhit
for i in $(ls 01_Preprocess/*.fastq.gz);do basename $i;done >$out_dir/01_Preprocess/files_clean.list
echo 'qsub -R y -l h_rt=23:59:59 -N cdHitDup_CIAD -o '$out_dir'/00_qsubout -t 1-'$(ls 01_Preprocess/*.gz|wc -l)':2 /share/DiskCoViGen/VAMIP/Programs/cdHitDupPair_list_v2.1.sh '$out_dir/01_Preprocess' files_clean.list '$out_dir/02_derep''|bash
a=$(ls 01_Preprocess/*.gz|wc -l);for x in $(grep -m 1 "Done\!" 00_qsubout/preProcess*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done));do echo 'grep -m 1 "\.'$x'$" <(ls 00_qsubout/preProcess*)';done|bash >00_qsubout/FAILED_preProcess.txt;cat 00_qsubout/FAILED_preProcess.txt
ls 01_Preprocess/*.gz|wc -l
# # Test
# for i in $(ls 01_Preprocess/*.gz);do basename $i;done|sed 's/_R.*//'|sort|uniq >01_Preprocess_out.txt
# for i in $(ls 02_derep/*.gz);do basename $i;done|sed 's/_R.*//'|sort|uniq >02_derep_out.txt
# Get which ones fail
a=$(ls 01_Preprocess/*.gz|wc -l);grep -m 1 "Done\!" 00_qsubout/cdHitDup*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done)|grep -f - <(ls 00_qsubout/cdHitDup*) >00_qsubout/FAILED_cdHitDup.txt;cat 00_qsubout/FAILED_cdHitDup.txt
# Assembly and analyze changes vs Original SARS-CoV-2 reference
for i in $(ls 02_derep/*.fastq.gz);do basename $i;done >$out_dir/02_derep/files_derep.list
echo 'qsub -R y -l h_rt=23:59:59 -N ConsesusVariant_CIAD -o '$out_dir'/00_qsubout -t 1-'$(ls 02_derep/*.gz|wc -l)':2 /share/DiskCoViGen/VAMIP/Programs/variantConsePairCall_List_v2.1.sh '$out_dir/02_derep' files_derep.list /share/DiskCoViGen/VAMIP/DB/MN908947.fasta genome '$out_dir/03_assembly''|bash
# Get which ones fail
a=$(ls 01_Preprocess/*.gz|wc -l);grep -m 1 "Done\!" 00_qsubout/ConsesusVariant*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done)|grep -f - <(ls 00_qsubout/ConsesusVariant*) >00_qsubout/FAILED_ConsensusVariant.txt; cat 00_qsubout/FAILED_ConsensusVariant.txt

# LANGEBIO
month="marzo-22"
out_dir="/scratch/rodrigog/01_Projects/VAMIP/04_fastq/temp/03_LANGEBIO/$month"; mkdir -p $out_dir; cd $out_dir
# Concatenate 4 lanes into a single fastq
in_dir="/scratch/rodrigog/01_Projects/VAMIP/04_fastq/temp/$month/reads"
for sam in $(ls $in_dir|grep fastq.gz|sed "s/_L00.*//"|sort|uniq);do for pair in R1 R2;do printf "cat"; for lane in L001 L002 L003 L004; do printf ' '$in_dir/$sam'_'$lane'_'$pair'_001.fastq.gz ';done;printf ' >'$out_dir'/'$sam'_'$pair'.fastq.gz\n';done;done >cat_files.sh
nohup bash cat_files.sh &
mv nohup.out FAILED_fastqs_missing_lanes.txt; cat FAILED_fastqs_missing_lanes.txt
cd $out_dir; mkdir -p 00_qsubout 01_Preprocess 02_derep 03_assembly
ls *.gz|sort >files.list
# Filter reads with fastp
echo 'qsub -R y -l h_rt=23:59:59 -N preProcess_LANGEBIO_'$month' -o '$out_dir'/00_qsubout -t 1-'$(ls *.gz|wc -l)':2 /share/DiskCoViGen/VAMIP/Programs/preProcess_list_v2.sh '$out_dir' files.list '$out_dir'/01_Preprocess'|bash
a=$(ls 01_Preprocess/*.gz|wc -l);for x in $(grep -m 1 "Done\!" 00_qsubout/preProcess*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done));do echo 'grep -m 1 "\.'$x'$" <(ls 00_qsubout/preProcess*)';done|bash >00_qsubout/FAILED_preProcess.txt;cat 00_qsubout/FAILED_preProcess.txt
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
ls 03_assembly/*.fasta|wc -l

# IBT
month="Lote47_IBT_23Marzo2022"
out_dir="/scratch/rodrigog/01_Projects/VAMIP/04_fastq/temp/01_IBT/$month"; mkdir -p $out_dir; cd $out_dir
# Concatenate 4 lanes into a single fastq
# in_dir="/share/DiskCoViGen/VAMIP/04_fastq/01_IBT/1000Genomas/Segunda"
in_dir="/share/DiskCoViGen/VAMIP/04_fastq/01_IBT/$month"
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
echo 'qsub -R y -l h_rt=23:59:59 -N ConsesusVariant_IBT_'$month' -o '$out_dir'/00_qsubout -t 1-'$(ls 02_derep/*.gz|wc -l)':2 /share/DiskCoViGen/VAMIP/Programs/variantConsePairCall_List_v2.1.sh '$out_dir/02_derep' files_derep.list /share/DiskCoViGen/VAMIP/DB/MN908947.fasta test '$out_dir/03_assembly''|bash
# Get which ones fail
a=$(ls 02_derep/*.gz|wc -l);for x in $(grep -m 1 "Done\!" 00_qsubout/ConsesusVariant*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done));do echo "\.$x$";done|grep -f - <(ls 00_qsubout/ConsesusVariant*) >00_qsubout/FAILED_ConsensusVariant.txt; cat 00_qsubout/FAILED_ConsensusVariant.txt
ls 03_assembly/*.fasta|wc -l

# INER
month="IMSS-5_2021-07-14"
out_dir="/scratch/rodrigog/01_Projects/VAMIP/04_fastq/temp/02_INER/$month"; mkdir -p $out_dir; cd $out_dir
# Concatenate 4 lanes into a single fastq
in_dir="/share/DiskCoViGen/VAMIP/04_fastq/02_INER/$month"
for sam in $(find $in_dir -name "*.fastq.gz");do base=$(basename $sam); printf "cat $sam >>$out_dir/";echo $base|sed "s/_L00[1-4]_\(R[12]\)_00[12]\.fastq\.gz/_\1.fastq.gz/";done >cat_files.sh # this assumes there are no repeated prefixes. otherwise, it will concatenate incorrect files
nohup bash cat_files.sh &
mv nohup.out FAILED_fastqs_missing_lanes.txt; cat FAILED_fastqs_missing_lanes.txt
mkdir -p 00_qsubout 01_Preprocess 02_derep 03_assembly
ls *.gz|sort >files.list
# Filter reads with fastp
echo 'qsub -R y -l h_rt=23:59:59 -N preProcess_INER_'$month' -o '$out_dir'/00_qsubout -t 1-'$(ls *.gz|wc -l)':2 /share/DiskCoViGen/VAMIP/Programs/preProcess_list_v2.sh '$out_dir' files.list '$out_dir'/01_Preprocess'|bash
a=$(ls 01_Preprocess/*.gz|wc -l);for x in $(grep -m 1 "Done\!" 00_qsubout/preProcess*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done));do echo 'grep -m 1 "\.'$x'$" <(ls 00_qsubout/preProcess*)';done|bash >00_qsubout/FAILED_preProcess.txt;cat 00_qsubout/FAILED_preProcess.txt
ls 01_Preprocess/*.gz|wc -l
# dereplicate with cdhit
for i in $(ls 01_Preprocess/*.fastq.gz);do basename $i;done >$out_dir/01_Preprocess/files_clean.list
echo 'qsub -R y -l h_rt=23:59:59 -N cdHitDup_INER_'$month' -o '$out_dir'/00_qsubout -t 1-'$(ls 01_Preprocess/*.gz|wc -l)':2 /share/DiskCoViGen/VAMIP/Programs/cdHitDupPair_list_v2.1.sh '$out_dir/01_Preprocess' files_clean.list '$out_dir/02_derep''|bash
# Get which ones fail
a=$(ls 01_Preprocess/*.gz|wc -l);for x in $(grep -m 1 "Done\!" 00_qsubout/cdHitDup*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done));do echo 'grep -m 1 "\.'$x'$" <(ls 00_qsubout/cdHitDup*)';done|bash >00_qsubout/FAILED_cdHitDup.txt;cat 00_qsubout/FAILED_cdHitDup.txt
# Assembly and analyze changes vs Original SARS-CoV-2 reference
for i in $(ls 02_derep/*.fastq.gz);do basename $i;done >$out_dir/02_derep/files_derep.list
echo 'qsub -R y -l h_rt=23:59:59 -N ConsesusVariant_INER_'$month' -o '$out_dir'/00_qsubout -t 1-'$(ls 02_derep/*.gz|wc -l)':2 /share/DiskCoViGen/VAMIP/Programs/variantConsePairCall_List_v2.1.sh '$out_dir/02_derep' files_derep.list /share/DiskCoViGen/VAMIP/DB/MN908947.fasta test '$out_dir/03_assembly''|bash
# Get which ones fail
a=$(ls 02_derep/*.gz|wc -l);for x in $(grep -m 1 "Done\!" 00_qsubout/ConsesusVariant*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done));do echo 'grep -m 1 "\.'$x'$" <(ls 00_qsubout/ConsesusVariant*)';done|bash >00_qsubout/FAILED_ConsensusVariant.txt;cat 00_qsubout/FAILED_ConsensusVariant.txt
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
# # # a=$(cat 02_derep/files_derep.list|wc -l);for x in $(grep -m 1 "Done\!" 00_qsubout/fix_ivar*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done));do echo "\.$x$";done|grep -f - <(ls 00_qsubout/fix_ivar*) >00_qsubout/FAILED_Fixivar.txt; cat 00_qsubout/FAILED_Fixivar.txt
# # # rm 03_assembly/*.bam
# # # nohup cp 03_assembly/ /share/DiskCoViGen/VAMIP/05_variants/01_IBT/$month/03_assembly/ &


#### This is done directly in the output folders in DiskCoViGen:
month="Noviembre25"
ins="03_LANGEBIO"
out_dir="/share/DiskCoViGen/VAMIP/05_variants/$ins/$month"; cd $out_dir
rm 03_assembly/*.gz # remove bad tables
# Assembly and analyze changes vs Original SARS-CoV-2 reference
echo 'qsub -R y -l h_rt=23:59:59 -N fix_ivar_LANGEBIO_'$month' -o '$out_dir'/00_qsubout -t 1-'$(cat 02_derep/files_derep.list|wc -l)':2 /share/DiskCoViGen/VAMIP/Programs/ivar_only_v1.sh '$out_dir/02_derep' files_derep.list /share/DiskCoViGen/VAMIP/DB/MN908947.fasta test '$out_dir/03_assembly''|bash
# Get which ones fail
a=$(cat 02_derep/files_derep.list|wc -l);for x in $(grep -m 1 "Done\!" 00_qsubout/fix_ivar*|sed -e 's/.*\.//' -e 's/:.*//'|grep -vwFf - <(for i in $(seq 1 2 $a);do echo $i;done));do echo "\.$x$";done|grep -f - <(ls 00_qsubout/fix_ivar*) >00_qsubout/FAILED_Fixivar.txt; cat 00_qsubout/FAILED_Fixivar.txt
ls 03_assembly/*.fasta|wc -l

# CIAD
out_dir="/share/DiskCoViGen/VAMIP/05_variants/04_CIAD"; cd $out_dir
rm 03_assembly/*.gz # remove bad tables
# Assembly and analyze changes vs Original SARS-CoV-2 reference
echo 'qsub -R y -l h_rt=23:59:59 -N fix_ivar_CIAD -o '$out_dir'/00_qsubout -t 1-'$(cat 02_derep/files_derep.list|wc -l)':2 /share/DiskCoViGen/VAMIP/Programs/ivar_only_v1.sh '$out_dir/02_derep' files_derep.list /share/DiskCoViGen/VAMIP/DB/MN908947.fasta test '$out_dir/03_assembly''|bash
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

# 04_CIAD
# All CIAD files were stored in the same folder as their nomenclature differs from CoViGen's. The renaming tables (dictionaries) were updated manually and stored @ /share/DiskCoViGen/VAMIP/05_variants/IdSampleRelations/04_CIAD/CIAD_2022-10-21.tsv
# File /share/DiskCoViGen/VAMIP/05_variants/04_CIAD/rename_template.sh has the necessary commands to process a two table dictionary with old name to new name. Use as follows:
bash rename_template.sh All_CIAD-2022-04-06 CIAD_2022-10-21.tsv
# And run the resulting renaming commands (errors are expected whenever the old and new name are the same)
bash rename_samples_All_CIAD-2022-04-06.sh
# Some files were missing from the dictionary, thus, we change them manually (these belong to sample 202201038655 and its resulting files)
for i in $(find All_CIAD-2022-04-06 -name "*test*"); do echo 'mv '$i' '$(echo $i|sed "s/_test//")'';done >remove_test_fag_CIAD-2022-04-06.sh

# RENAME FILES TO ADD BATCH NUMBER
cd /share/DiskCoViGen/VAMIP/05_variants
# P1000, Cardio and CIENI were transfered to folder: /share/DiskCoViGen/VAMIP/NON_USED
mkdir /share/DiskCoViGen/VAMIP/NON_USED
mv /share/DiskCoViGen/VAMIP/05_variants/01_IBT/P1000_1 /share/DiskCoViGen/VAMIP/NON_USED
mv /share/DiskCoViGen/VAMIP/05_variants/01_IBT/P1000_2 /share/DiskCoViGen/VAMIP/NON_USED
mv /share/DiskCoViGen/VAMIP/05_variants/02_INER/Cardio /share/DiskCoViGen/VAMIP/NON_USED
mv /share/DiskCoViGen/VAMIP/05_variants/02_INER/CIENI /share/DiskCoViGen/VAMIP/NON_USED
mv /share/DiskCoViGen/VAMIP/05_variants/03_LANGEBIO/Noviembre25 /share/DiskCoViGen/VAMIP/NON_USED

# Gzip all fastas and tables
cd /share/DiskCoViGen/VAMIP/05_variants
find . -name "*.fasta" -exec gzip "{}" \;
find . -name "*_Q_Dp.ivar.tsv" -exec gzip "{}" \;

# 2022-04-27: We decided to add an additional identifier to the filenames like L001 for batch1, and so on
for i in $(find .|grep "fastq\|_Q_Dp");do echo 'mv '$i' '$(echo $i|sed "s/Lote\(\w\w\)\(.*\)\//Lote\1\2\/L0\1_/")'';done >rename_add_batches.sh

# Summarize contents
# Make sure the string "test" is not present in the file names (this means names have not been updated)
# We now need a list of all files having an ivar table
# for i in $(cat 2022-02-25.list);do base=$(basename $i);echo ''${base%_Q_Dp.ivar2.fasta}'';done >2022-02-25_samples_with_ivar.list
cd /share/DiskCoViGen/VAMIP/05_variants
# This will exclude all items whose names have not been processed (those flagged with string "test")
find ~+ -name "*.ivar.tsv.gz"|grep -v test >2022-04-28_ivar_tables.list # The ~+ expands the current path to full path
cat 2022-04-28_ivar_tables.list| grep CIAD >2022-04-28_CIAD_only.tsv # This is latter used for a filter
# We now need the nucleotide contents
# File 2022-04-28_CIAD_wFolio.tsv has those new id folios. With this, we can remove the unwanted CIAD ones
for i in $(cat 2022-04-28_CIAD_wFolio.tsv);do echo 'grep -m 1 "03_assembly/'$i'_S" 2022-04-28_CIAD_only.tsv';done|bash|grep -vwFf - 2022-04-28_CIAD_only.tsv >2022-04-28_ignore_CIAD_files.tsv

find ~+ -name "*.ivar2.fasta.gz"|grep -v test >2022-04-28_fastas.list # List all target consensus (ignore those flagged tests)
# grep "04_CIAD\|/202" 2022-03-03_ivar.list >2022-03-03_ivar.list.fix # Keep only those identified with a folio id or CIAD id
# Using the list of files, we get which ones were the original (raw, just after joinning lanes) files
find ~+ -name "*fastq.gz"|grep -v "01_Preprocess\|02_derep\|03_assembly" >2022-04-28_raw_fastqs.list
# We can now extract out of the original files, how many sequences they had
for i in $(cat 2022-04-28_raw_fastqs.list);do echo ls -lh $i;done|bash >2022-04-28_rawfastq_size.tsv
# This was carried fo R1 and R2 because there are very rare cases where they don't match in total seqs
grep _R1.fastq.gz 2022-04-28_rawfastq_size.tsv >2022-04-28_sample_sizeR1.tsv
# For the raw files, we'll create a catalog with batch, Institute and size. Id is the folio IMSS ID (e.g. 202201000977) which is appendend as the las column (8)
paste <(cut -f 2 2022-04-28_sample_sizeR1.tsv) <(sed -e 's/\t.*academico /\t/' -e 's/\//\t/g' -e 's/ /\t/' -e 's/\t\(202.*\)\(_S.*\)/\t\1\2\t\1/' 2022-04-28_sample_sizeR1.tsv) >2022-04-28_sample_ok.tsv
# DEPRECATED: START
# # # # # Now, we discard those not having a valid IMSS ID
# grep "2021/2021\|2022/2022\|2022/2021\|P2/202\|2022/2.\|CIAD" 2022-04-28_sample_ok.tsv|sort >2022-04-28_sample_ok_target.list.tsv
# # # # # And keep track of those that didn't
# grep -v "2021/2021\|2022/2022\|2022/2021\|P2/202\|2022/2\|CIAD." 2022-04-28_sample_ok.tsv|sort >2022-04-28_sample_ok_target.bad.tsv
# DEPRECATED: END
cp 2022-04-28_sample_ok.tsv 2022-04-28_sample_ok_target.list.tsv

# The following commands should be optimized to run in parallel (maybe later)
#IMPORTANT: divisions by 0 fail and are seen as empty rows with only the filename
for i in $(cat 2022-04-28_fastas.list);do echo 'printf "'$i'\t" >>2022-04-28_nucleotide_contents.tsv;perl -ne '\''{next if $_=~ /^>/; chomp($_); $length=length($_); $Ns=tr/[N]//;$As=tr/[A]//; $Ts=tr/[T]//; $Cs=tr/[C]//; $Gs=tr/[G]//; $ATs=tr/[AT]//; $CGs=tr/[CG]//; $nonNs=tr/[ATCG]//; $Np=$Ns*100/$length; $nonNp=$nonNs*100/$length; $Ap=$As*100/$length; $Tp=$Ts*100/$length; $Cp=$Cs*100/$length; $Gp=$Gs*100/$length; $ATp=$ATs*100/$nonNs; $CGp=$CGs*100/$nonNs; $out=sprintf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f", $length, $nonNs, $Ns, $As, $Ts, $Cs, $Gs, $ATs, $CGs, $nonNp, $Np, $Ap, $Tp, $Cp, $Gp, $ATp, $CGp); print $out}'\'' <(zcat '$i') >>2022-04-28_nucleotide_contents.tsv; printf "\n" >>2022-04-28_nucleotide_contents.tsv';done >summarize_contents.sh
# Extract only those with non-empty items
awk -F'\t' '$2!=""' 2022-04-28_nucleotide_contents.tsv >2022-04-28_nucleotide_contents_nonEmtpy.tsv
# Now, append the IMSS folio as an extra column
paste 2022-04-28_nucleotide_contents_nonEmtpy.tsv <(sed -e "s/\//\t/g" 2022-04-28_nucleotide_contents_nonEmtpy.tsv -e 's/_S.*//'|cut -f 9) >2022-04-28_nucleotide_contents_wID.tsv
# Additionally, we create a list of those having >25% Ns that are targeted for removal later on
awk -F'\t' '$12>=25' 2022-04-28_nucleotide_contents_wID.tsv|cut -f 1|sed -e 's/.*03_assembly\///' -e 's/\(_Q_Dp\).*//' >2022-04-28_for_removal_gt_25perc_Ns.tsv
awk -F'\t' '$12>=15' 2022-04-28_nucleotide_contents_wID.tsv|cut -f 1|sed -e 's/.*03_assembly\///' -e 's/\(_Q_Dp\).*//' >2022-04-28_for_removal_gt_15perc_Ns.tsv
awk -F'\t' '$12>=10' 2022-04-28_nucleotide_contents_wID.tsv|cut -f 1|sed -e 's/.*03_assembly\///' -e 's/\(_Q_Dp\).*//' >2022-04-28_for_removal_gt_10perc_Ns.tsv
# With this, we create a list of all available genomes that are part of the surveillance batches (samples only, no controls) # IMPORTANT: This list contains the working samples for the downstream analyses
cut -f 10 2022-04-28_sample_ok_target.list.tsv|grep -v "^$"|sed -e 's/_R1.fastq.gz//' -e 's/^/03_assembly\//'|grep -Ff - 2022-04-28_nucleotide_contents_wID.tsv|cut -f 19 >2022-04-28_Main_mapped-Vigilancia.list
# BAD cut -f 19 2022-04-28_nucleotide_contents_wID.tsv|grep -v "^$"|grep -wFf - 2022-04-28_sample_ok_target.list.tsv|cut -f 8 >2022-04-28_Main_mapped-Vigilancia.list

# With both the Main_mapped_-Vigilancia.list and the 2022-04-28_sample_ok_target.list.tsv tables, we can now append size information to the nucleotide contents table using 2022-04-28_nucleotide_contents_wID.tsv as xref table
printf "File\tgenome_size\tNonNs\tNs\tAs\tTs\tCs\tGs\tATs\tCGs\tnonN%%\tN%%\tA%%\tT%%\tC%%\tG%%\tAT%%\tCG%%\tGenome\tPerm\tFile_size\n" >2022-04-28_nucleotide_contents_wSize.tsv
for i in $(cat 2022-04-28_Main_mapped-Vigilancia.list);do paste <(grep -m 1 "/$i\_S" 2022-04-28_nucleotide_contents_wID.tsv) <(grep -m 1 "/$i\_S" 2022-04-28_sample_ok_target.list.tsv|cut -f 2,3|cut -d " " -f 4);done >>2022-04-28_nucleotide_contents_wSize.tsv
# To decide upon this matter, we graph some an Rscript to graph some features


# ### Summarize ivar tables ###
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
# outdir="/scratch/rodrigog/01_Projects/VAMIP/06_VAMIP_tables/01_per_sample";for i in $(cat 2022-04-28_ivar_tables.list);do base=$(basename $i);echo 'printf "'$(echo $i|sed -e "s/.*05_variants\///" -e "s/03_assembly\///" -e "s/\//\\\\t/g" -e "s/_S[0-9]*_Q_Dp.ivar.tsv.gz//" -e "s/_Q_Dp_genome.ivar.tsv.gz//")'\t" >'$outdir'/'${base%.ivar.tsv.gz}'_vamip.tsv; awk '\''NR>=2 {new_var=">"$2"M"$3"N"$4"F"$11"X"$12"R"$5"L"$9" "; print new_var}'\'' <(zcat '$i')|tr . d|tr - m|tr + p|tr -d "\n" >>'$outdir/${base%.ivar.tsv.gz}'_vamip.tsv;printf "\n" >>'$outdir/${base%.ivar.tsv.gz}'_vamip.tsv';done >create_single_sample_vamip_tables.sh
outdir="/scratch/rodrigog/01_Projects/VAMIP/06_VAMIP_tables/01_per_sample";for i in $(cat 2022-04-28_ivar_tables.list);do base=$(basename $i);echo 'printf "'$(echo $i|sed -e "s/.*05_variants\///" -e "s/03_assembly\///" -e "s/\//\\\\t/g" -e "s/_Q_Dp.ivar.tsv.gz//" -e "s/_Q_Dp_genome.ivar.tsv.gz//")'\t" >'$outdir'/'${base%_Q_Dp.ivar.tsv.gz}'_vamip.tsv; awk '\''NR>=2 {new_var=">"$2"M"$3"N"$4"F"$11"X"$12"R"$5"L"$9" "; print new_var}'\'' <(zcat '$i')|tr . d|tr - m|tr + p|tr -d "\n" >>'$outdir/${base%_Q_Dp.ivar.tsv.gz}'_vamip.tsv;printf "\n" >>'$outdir/${base%_Q_Dp.ivar.tsv.gz}'_vamip.tsv';done >create_single_sample_vamip_tables.sh
# qsub script qsub_command_per_line.sh was created to run all lines in parallel (any bash script for one command per line)
mkdir -p $outdir/00_qsubout
qsub -R y -l h_rt=23:59:59 -N vamip_tables -o $outdir/00_qsubout -t 1-$(cat create_single_sample_vamip_tables.sh|wc -l) qsub_command_per_line.sh /share/DiskCoViGen/VAMIP/05_variants/create_single_sample_vamip_tables.sh
mv 01_per_sample/00_qsubout/ # I decided to move the logs afterwards

# IMPORTANT: FILTER UNWANTED ITEMS HERE
# /scratch/rodrigog/01_Projects/VAMIP/06_VAMIP_tables/2022-04-28_ignore_CIAD_files.tsv (now locally) has a list of files to ignore from CIAD.
sed -e 's/_Q_Dp.ivar.tsv.gz//' -e 's/.*\///' -e 's/$/_vamip.tsv/' 2022-04-28_ignore_CIAD_files.tsv >bad_vamip.list
ls 01_per_sample >files.list
grep "_CN" files.list >>bad_vamip.list
grep "egativ" files.list >>bad_vamip.list
grep "EGATIV" files.list >>bad_vamip.list
grep "ILLUMINA" files.list >>bad_vamip.list
grep "_A_" files.list >>bad_vamip.list
grep "_B_" files.list >>bad_vamip.list
grep "_C_" files.list >>bad_vamip.list
grep "_D_" files.list >>bad_vamip.list
grep "_E_" files.list >>bad_vamip.list
grep "_F_" files.list >>bad_vamip.list
grep "_G_" files.list >>bad_vamip.list
grep "_H_" files.list >>bad_vamip.list
grep "_I_" files.list >>bad_vamip.list
grep "_J_" files.list >>bad_vamip.list
grep "_K_" files.list >>bad_vamip.list
grep "_L_" files.list >>bad_vamip.list
grep "^EX" files.list >>bad_vamip.list
grep "_INER-5-" files.list >>bad_vamip.list
grep "202201004764_S6_vamip.tsv" files.list >>bad_vamip.list
grep "_CP._" files.list >>bad_vamip.list

# Repeated files were also excluded
mkdir 01_per_sample-ignored 01_per_sample_gt_25perc_Ns 01_per_sample_gt_15perc_Ns 01_per_sample_gt_10perc_Ns
# Now, filter those that were selected
for i in $(cat bad_vamip.list);do echo 'mv 01_per_sample/'$i' 01_per_sample-ignored';done|bash
# And those with more than 25%, 15% and 10% Ns
for i in $(cat 2022-04-28_for_removal_gt_25perc_Ns.tsv);do echo 'mv 01_per_sample/'$i'_vamip.tsv 01_per_sample_gt_25perc_Ns';done|bash
for i in $(cat 2022-04-28_for_removal_gt_15perc_Ns.tsv);do echo 'mv 01_per_sample/'$i'_vamip.tsv 01_per_sample_gt_15perc_Ns';done|bash
for i in $(cat 2022-04-28_for_removal_gt_10perc_Ns.tsv);do echo 'mv 01_per_sample/'$i'_vamip.tsv 01_per_sample_gt_10perc_Ns';done|bash
# IMPOTANT NOTE: Only genomes bearing <10% Ns were kept for dowstream analyses

# patch the CIAD files as some still bear an R2 string
for i in $(find 06_VAMIP_tables/01_per_sample -name "*R2_vamip.tsv");do echo 'mv '$i' '$(echo $i|sed "s/_R2_vamip/_vamip/")'';done|bash
# CIAD's files have not been renamed to include the batch, use the following to change this (this was done locally)
cd /home/rod/Documents/01_Projects/SARS/VAMIP
# Fist get the list of CIAD files and keep only the column labeled "folio"
grep "ciad\|Ciad\|CIAD" 2022_04_27_DB_MexCov2_v1.tsv|cut -f 13|sed '/^$/d' >2022_04_27_DB_MexCov2_v1-OnlyCIAD.tsv
# Then, trim the L0XX part and keep only those that are not repeated (since we cannot know which ones match these.
sed 's/L0[0-9][0-9]_//'  2022_04_27_DB_MexCov2_v1-OnlyCIAD.tsv|sort|uniq -c|tr -s " "|sed "s/ /\t/g"|awk -F'\t' '$2==1'|cut -f 3 >2022-04-27_unique_CIAD.txt
# Now, get those of these that are present our actual set (as vamip files) and should be appended a batch name
for i in $(cat 2022-04-27_unique_CIAD.txt);do echo 'printf "\n'$i'\t"; grep -m 1 "^'$i'_" <(ls 06_VAMIP_tables/01_per_sample)';done|bash|sed '/^$/d'|cut -f 2|sed '/^$/d' >2022-04-27_actual_CIAD_vamips.tsv
# Do both searches and output to temporary file:
for i in $(sed 's/_S[0-9]*_vamip.tsv//' 2022-04-27_actual_CIAD_vamips.tsv);do echo 'grep "^'$i'_S" 2022-04-27_actual_CIAD_vamips.tsv; grep "_'$i'$" 2022_04_27_DB_MexCov2_v1-OnlyCIAD.tsv';done|bash|paste - - >temp
paste -d'\0' temp <(sed "s/.*\(S[0-9]*_vamip.tsv\).*/_\1/" temp)|sed -e 's/^/mv /' -e 's/\t/ /' >2022-04-27_add_Batches_to_CIAD.sh
cd /home/rod/Documents/01_Projects/SARS/VAMIP/06_VAMIP_tables/01_per_sample # Locally


# Now survey each position
cd /scratch/rodrigog/01_Projects/VAMIP/06_VAMIP_tables/
mkdir -p /scratch/rodrigog/01_Projects/VAMIP/06_VAMIP_tables/02_mut_by_pos
# for i in {1..30000};do echo 'egrep -m 1 -o ">'$i'M\w*" 01_per_sample/*|sed -e "s/.*01_per_sample\///" -e "s/_vamip.tsv:>/\t/" -e "s/_S[0-9]*\t/\t/" >02_mut_by_pos/pos'$i.tsv'';done >mutations_by_position.sh
# for i in {1..30000};do echo 'egrep -m 1 -o ">'$i'M\w*" 01_per_sample/*|sed -e "s/.*01_per_sample\///" -e ""s/_R2_Q/_Q/ -e "s/_Q_Dp_.*vamip.tsv:>/\t/" >02_mut_by_pos/pos'$i.tsv'';done >mutations_by_position.sh
for i in {1..30000};do echo 'egrep -m 1 -o ">'$i'M\w*" 01_per_sample/*|sed -e "s/.*01_per_sample\///" -e "s/_R2_/_/" -e "s/_vamip.tsv:>/\t/" >02_mut_by_pos/pos'$i.tsv'';done >mutations_by_position.sh
# cd /share/DiskCoViGen/VAMIP/05_variants
qsub -R y -l h_rt=23:59:59 -N SNP_positions -o 00_qsubout -t 1-$(cat mutations_by_position.sh|wc -l) ~/bin/qsub_command_per_line.sh /scratch/rodrigog/01_Projects/VAMIP/06_VAMIP_tables/mutations_by_position.sh

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
# As of 2022-04-28, there were 28,684 (95.92%) positions that report any mutation
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
# from the remaining items, 2,881 (9.63%) had at least 2 observations
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


