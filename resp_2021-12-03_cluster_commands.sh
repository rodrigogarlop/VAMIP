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
mkdir /scratch/rodrigog/01_Projects/VAMIP/04_fastq/temp
cp -r /share/DiskCoViGen/VAMIP/04_fastq/04_CIAD/ ./04_fastq/temp/
for i in $(find . -name *.gz);do echo mv $i /scratch/rodrigog/01_Projects/VAMIP/04_fastq/temp/04_CIAD/;done|bash # move all fastqs to the root CIAD dir
for i in $(ls|grep -v ".gz");do echo rmdir $i;done|grep -v original|bash # Cleanse empty folders
rename ".R1" '_R1' *.gz # Names from CIAD files have a uncommon .R1 separator, which is not compatible with our current scripts
rename ".R2" '_R2' *.gz
mkdir -p 00_qsubout 01_Preprocess 02_derep 03_assembly
# Filter reads with fastp
out_dir="/scratch/rodrigog/01_Projects/VAMIP/04_fastq/temp/04_CIAD"
echo 'qsub -R y -l h_rt=23:59:59 -N preProcess -o '$out_dir'/00_qsubout -t 1-808:2 /share/DiskCoViGen/VAMIP/Programs/preProcess_list_v2.sh '$out_dir' files.list '$out_dir'/01_Preprocess'
# dereplicate with cdhit
for i in $(ls 01_Preprocess/*.fastq.gz);do basename $i;done >$out_dir/01_Preprocess/files_clean.list
echo 'qsub -R y -l h_rt=23:59:59 -N cdHitDup_CIAD -o '$out_dir'/00_qsubout -t 1-'$(ls 01_Preprocess/*.gz|wc -l)':2 /share/DiskCoViGen/VAMIP/Programs/cdHitDupPair_list_v2.1.sh '$out_dir/01_Preprocess' files_clean.list '$out_dir/02_derep''
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
month="julio1"
out_dir="/scratch/rodrigog/01_Projects/VAMIP/04_fastq/temp/03_LANGEBIO/$month"; mkdir -p $out_dir; cd $out_dir
# Concatenate 4 lanes into a single fastq
in_dir="/share/DiskCoViGen/VAMIP/04_fastq/03_LANGEBIO/$month/reads"
for sam in $(ls $in_dir|grep fastq.gz|sed "s/_L00.*//"|sort|uniq);do for pair in R1 R2;do printf "cat"; for lane in L001 L002 L003 L004; do printf ' '$in_dir/$sam'_'$lane'_'$pair'_001.fastq.gz ';done;printf ' >'$out_dir'/'$sam'_'$pair'.fastq.gz\n';done;done >cat_files.sh
cd $out_dir; mkdir -p 00_qsubout 01_Preprocess 02_derep 03_assembly
ls *.gz >files.list
# Filter reads with fastp
echo 'qsub -R y -l h_rt=23:59:59 -N preProcess_LANGEBIO_'$month' -o '$out_dir'/00_qsubout -t 1-'$(ls *.gz|wc -l)':2 /share/DiskCoViGen/VAMIP/Programs/preProcess_list_v2.sh '$out_dir' files.list '$out_dir'/01_Preprocess'|bash
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
