# Started 2022-04-29 by Rodrigo García-López for Carlos Arias's Viromics Group at IBt, UNAM, Cuernavaca, Mexico as part of the CoViGen-Mex SARS-CoV-2 surveillance in Mexico.
# "VAMIP" is a CONACYT project that is framed within the SARS-CoV-2 surveillance initiative CoViGen Mex aiming at determining the role of intra-patient minor allellic variants (IPMAVs in English, VAMIP in Spanish) using complete viral genome sequences produced by our group and collaborators. Fastq or bam files are required to assess the different basecalls per position where no consensus is met in the genomic sequence.
# This is an upgraded version containing local commands. Previous commands were executed in cluster teopanzolco.ibt.unam.mx by usr rodrigog

# Deprecated exploration commands (folders 01-03 in /home/rod/Documents/01_Projects/SARS/VAMIP/) are found at /home/rod/Documents/01_Projects/SARS/VAMIP/00_docs/IPMAV_Commands_First_Year_exploration.sh
# Actual commands that should be used before this file (most were executed in our cluster) /home/rod/Documents/01_Projects/SARS/VAMIP/00_docs/cluster_commands.sh
# A full backup of the CoViGen files can be found in the cluster, at /share/DiskCoViGen/VAMIP/05_variants
# Briefly, raw fastq files from amplicon high-throuput sequencing data were cleansed, derreplicated and reads were mapped against the primal Wuhan SARS-CoV-2 sequence. ivar tables were extracted to trace minor allellic variants. These were further processed into coded versions per sample and per position, which was then integrated into a single table at /home/rod/Documents/01_Projects/SARS/VAMIP/06_VAMIP_tables/All_mut_NoSingl.tsv which excludes positions with less than 5 mutation ovservations (more than 5 genomes in the set have them).


# ### PART 1: Input files ###
# A) MUTATIONS TABLE
# Input file location:
# Mutations are found per genome in files found at /home/rod/Documents/01_Projects/SARS/VAMIP/06_VAMIP_tables/01_per_sample/. These include only those that were not filtered in previous steps (some other, like controls, repeated ones, those with >10% Ns, and some missing metadata were not inclued in the main folder)
# Ignored items
# file:///home/rod/Documents/01_Projects/SARS/VAMIP/06_VAMIP_tables/01_per_sample-ignored
# file:///home/rod/Documents/01_Projects/SARS/VAMIP/06_VAMIP_tables/01_per_sample_gt_10perc_Ns
# file:///home/rod/Documents/01_Projects/SARS/VAMIP/06_VAMIP_tables/01_per_sample_gt_15perc_Ns
# file:///home/rod/Documents/01_Projects/SARS/VAMIP/06_VAMIP_tables/01_per_sample_gt_25perc_Ns
cd /home/rod/Documents/01_Projects/SARS/VAMIP/06_VAMIP_tables/01_per_sample


cd /home/rod/Documents/01_Projects/SARS/VAMIP/
ls 06_VAMIP_tables/01_per_sample/*.tsv|wc -l
# 18427
# Mutations per position were calculated from the passing files (see above). These are found at /home/rod/Documents/01_Projects/SARS/VAMIP/06_VAMIP_tables/02_mut_by_pos/ Positions with 0 genomes with mutations, single observations or less than 5 items were separated into their own folders (non-overlapping)

# Ignored items
# file:///home/rod/Documents/01_Projects/SARS/VAMIP/06_VAMIP_tables/02_mut_by_pos/items_2-4
# file:///home/rod/Documents/01_Projects/SARS/VAMIP/06_VAMIP_tables/02_mut_by_pos/singletons
# file:///home/rod/Documents/01_Projects/SARS/VAMIP/06_VAMIP_tables/02_mut_by_pos/zeromut

# The largest genome size in the current set has a length of 29,909nt, so we deleted non-extant items with
# for i in {29910..30000};do echo rm 02_mut_by_pos/pos$i.tsv;done|bash

# The resulting positions contain only positions with >=5 mutations. This is summarized in the following table as follows:
cat 02_mut_by_pos/*.tsv >temp; paste <(cut -f 1 temp) <(cut -f 2 temp|sed -e 's/d/\./' -e 's/p/+/' -e 's/m/-/' -e 's/[MNFXRL]/\t/g') >All_mut_NoSingl.tsv;rm temp
# All_mut_NoSingl.tsv was the input file for downstream analyses

# b) GENOME METADATA (MEXICO)
# IMPORTANT 2022-04-29: THIS TABLE IS NOT THE USUAL 27 COLUMN METADATA THAT I USE, SO THE FOLLOWING COMMANDS SHOULD BE ADJUSTED ACCORIDNG TO THE ACTUAL COLUMNS
# Process (CoViGen's) Mauricio's metadata (reported weekly). Please update accordingly
# File was downloaded into this folder and headers were fixed to remove rare characters
wc -l 2022_04_27_DB_MexCov2_v1.tsv
# 58728
# Among our genome vamip tables, there are some files not having an assigned batch (from CIAD). These cannot be matched as metadata does not include any kind of IMSS folio. We will first fix this accordingly. To do this, we search all items with no batch (starting with string L0XX_), then remove the SXXX identifier from the sequencing run.
ls 06_VAMIP_tables/01_per_sample|grep -v "^L0"|sed 's/_S[0-9]*_vamip.tsv//' >06_VAMIP_tables/2022-04-28_Vamip_files_with_no_batch.list
# Now extract the ones that match in the metadata:
for i in $(cat 06_VAMIP_tables/2022-04-28_Vamip_files_with_no_batch.list);do echo 'printf "\n'$i'\t";grep "'$i'" 2022_04_27_DB_MexCov2_v1.tsv';done|bash|sed '/^$/d'|sed 's/^EPI_I/\tEPI_I/' >06_VAMIP_tables/2022_04_27_DB_MexCov2_v1-unmatched.tsv
# This file is used to trace all items matching those genome names (all from CIAD). Since other matches point towards InDRE, these are ignored.
for i in $(cat 06_VAMIP_tables/2022-04-28_Vamip_files_with_no_batch.list);do echo 'printf "\n'$i'\t";grep "'$i'" <(grep -v "InDRE" 2022_04_27_DB_MexCov2_v1.tsv)';done|bash|sed '/^$/d'|sed 's/^EPI_I/\tEPI_I/' >06_VAMIP_tables/2022_04_27_DB_MexCov2_v1-unmatched-fix.tsv
# The result should now have a single hit per genome.
paste <(cut -f 2-13 06_VAMIP_tables/2022_04_27_DB_MexCov2_v1-unmatched-fix.tsv) <(cut -f 1 06_VAMIP_tables/2022_04_27_DB_MexCov2_v1-unmatched-fix.tsv) <(cut -f 15- 06_VAMIP_tables/2022_04_27_DB_MexCov2_v1-unmatched-fix.tsv) >temp
cat <(cut -f 1 temp|grep -vwFf - 2022_04_27_DB_MexCov2_v1.tsv) temp >2022_04_27_DB_MexCov2_v1_xtras.tsv
# Next, preprocess the table:
mkdir 07_metadata
# File EstadoRegion_v4.tsv was stored in that folder
Rscript preFilter_raw_metatable.R
Rscript Filter_preTable.R
# Folios L011_202101008588 and L035_2840 are repeated and will thus be ignored in dowstream analyses

# The next commands are in R at /home/rod/Documents/01_Projects/SARS/VAMIP/explore_vamips.R
