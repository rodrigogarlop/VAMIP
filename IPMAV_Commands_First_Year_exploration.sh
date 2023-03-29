# Started: 2021-10-28
# by Rodrigo García-López for Carlos Arias's Viromics Group at IBt, UNAM, Cuernavaca, Mexico as part of the CoViGen-Mex SARS-CoV-2 surveillance in Mexico.
# "VAMIP" is a CONACYT project that is framed within the SARS-CoV-2 surveillance initiative CoViGen Mex aiming at determining the role of intra-patient minor allellic variants (IPMAVs in English, VAMIP in Spanish) using complete viral genome sequences produced by our group and collaborators. Fastq or bam files are required to assess the different basecalls per position where no consensus is met in the genomic sequence.

# ### 01: PRE-PROCESSING AVAILABLE FILES ###
# LOCAL PC @ /home/rod/Documents/01_Projects/SARS/VAMIP/
# Update 2021-10-28: The starting point will be determining which genomes are available for use in this project or may be provided by our collaborators
# The first task wil be to assess which sequence sets (fastq or bam) we have in storage (IBt and collaborators) and which of these can be used to explore IPMAVs
# We've been provided a list of potential sequences from two main projects:
# 1) P1000 - Genomes from the first year of COVID-19 pandemic
# 2) Vigilancia - Genomes from the second year in the pandemic (ongoing)
# * only complete genomes (>97% coverage) will be used
#	P1000: First year genomes
# 	Two files were provided (downloaded from dropbox's VAMIP/DATA), located at 01_GISAID_Tables/2021-10-28_P1000
# 		20210319_P1000_CoV_Mex_P1Complete97.xls (2020-03-18 to 2020-09-08)
# 		20210319_P1000_CoV_Mex_P2Complete97.xls (2020-08-04 to 2021-02-03)
# 	None of this has lineage but this have been submitted to GISAID. Thus, we downloaded the metadata from GISAID (https://www.epicov.org/epi3/) -> Search -> Custom Selection -> Virus name: hCoV-19/Mexico/, collection dates from 2003-10-31 to 2021-02-03 (4229 results with 123.8 MiB, max is 5k) -> Select them -> Download (right icon) -> Input for Augur pipeline -> Download (gets a tar file containing both a fasta and a tsv table.
# PART I: 01_GISAID_Tables/2021-10-28_P1000/gisaid_auspice_input_hcov-19_2021_10_28_19/1635448214485.metadata.tsv (metadata table)
# With this, we now get the lineaje for the items that are included
cd /home/rod/Documents/01_Projects/SARS/VAMIP/
in_dir="01_GISAID_Tables/2021-10-28_P1000"
mkdir $in_dir
for i in 1 2;do echo 'touch '$in_dir'/P1000_P'$i'_query.tsv';done|bash # Create the list containers
# The Virus Name column was copied to their corresponding tsv files we've just created
# hCoV-19/Mexico/CMX-INER-IBT-9/2020's date was fixed from 2004-04-03 to 2020-04-03
for i in 1 2;do echo 'tail +3 '$in_dir'/P1000_P'$i'_query.tsv|wc -l';done|bash # Now, count the total items
# 295 01_GISAID_Tables/2021-10-28_P1000/P1000_P1_query.tsv
# 259 01_GISAID_Tables/2021-10-28_P1000/P1000_P2_query.tsv
# Now, search in the metadata for the actual items
for i in 1 2;do echo 'grep -wFf <(cut -f 3 '$in_dir'/P1000_P'$i'_query.tsv) '$in_dir'/gisaid_auspice_input_hcov-19_2021_10_28_19/1635448214485.metadata.tsv|cut -f 1,3,19 > '$in_dir'/gisaid_auspice_input_hcov-19_2021_10_28_19/P1000_P'$i'_found.tsv; wc -l '$in_dir'/gisaid_auspice_input_hcov-19_2021_10_28_19/P1000_P'$i'_found.tsv';done|bash # Search for the items in the tables, and summarize those that were found
# 290 01_GISAID_Tables/2021-10-28_P1000/gisaid_auspice_input_hcov-19_2021_10_28_19/P1000_P1_found.tsv
# 250 01_GISAID_Tables/2021-10-28_P1000/gisaid_auspice_input_hcov-19_2021_10_28_19/P1000_P2_found.tsv
# Extra: We can then get those that fail for testing (14 of them failed)
for i in 1 2;do echo 'grep -vwFf <(cut -f 1 '$in_dir'/gisaid_auspice_input_hcov-19_2021_10_28_19/P1000_P'$i'_found.tsv) '$in_dir'/P1000_P'$i'_query.tsv >'$in_dir'/gisaid_auspice_input_hcov-19_2021_10_28_19/P1000_P'$i'_notfound.tsv';done|bash
# Now, subset the gisaid table based on those items that were found and append the lineage column
for i in 1 2;do echo 'tail +2 '$in_dir'/P1000_P'$i'_query.tsv|head -n 1|sed "s/Comment\tComment Icon/Virus Name(should be the same as col 3)\tAccession ID\tLineage/" >'$in_dir'/P1000_P'$i'_GISAID_w_lineage.tsv; grep -wFf <(cut -f 1 '$in_dir'/gisaid_auspice_input_hcov-19_2021_10_28_19/P1000_P'$i'_found.tsv) '$in_dir'/P1000_P'$i'_query.tsv|cut -f 1-28|sort -k 3|paste - <(sort '$in_dir'/gisaid_auspice_input_hcov-19_2021_10_28_19/P1000_P'$i'_found.tsv) >>'$in_dir'/P1000_P'$i'_GISAID_w_lineage.tsv';done|bash
# The resulting table is the same as the GISAID original but has lineage in the last column (instead of comment). I decided to leave the virus name identifier in the previous column as well in case match verification is required
# Update 2021-11-03: We now concatenate both tables
in_dir="01_GISAID_Tables/2021-10-28_P1000"
outfile="01_GISAID_Tables/2021-10-28_First_year_All_GISDAID_metadata.tsv"
printf "Part (file)\t" >$outfile
head -n 1 $in_dir/P1000_P1_GISAID_w_lineage.tsv >>$outfile
for i in $(ls $in_dir/*lineage.tsv);do base=$(basename $i);echo 'tail +2 '$i'|awk '\''BEGIN{FS=OFS="\t"}{print "'${base%_w_lineage.tsv}'" OFS $0}'\''';done|bash >>$outfile

# and create a subset for plotting
infile="01_GISAID_Tables/2021-10-28_First_year_All_GISDAID_metadata.tsv"
printf "Accession ID\tCollection date\tLineage\tLocation\n" >01_GISAID_Tables/2021-10-28_First_year_for_areaplots.tsv
paste <(cut -f 31 $infile) <(cut -f 7 $infile) <(cut -f 32 $infile) <(cut -f 8 $infile)|sed -e 's/North America \/ Mexico \/ //' -e 's/\(202[012]\)-\(\w\w\)-\(\w\w\)/\3\/\2\/\1/'|tail +2 >>01_GISAID_Tables/2021-10-28_First_year_for_areaplots.tsv

# Update 2021-10-28: Next, filter out the different tables fron the surveillance project:
# PART II: 01_GISAID_Tables/2021-10-28_Vigilancia
in_dir="01_GISAID_Tables/2021-10-28_Vigilancia"
mkdir $in_dir
# This project is divided in batches. As of 2021-10-29, there are 26 that have been processed by our laboratory or collaborators.
# First, create an empty container for the metadata
touch $in_dir/2021-10-25_Metadata.tsv
# Now copy contents of table 2021-10-28_CoViGen-Mex/TodosGenomasGISAID/25Octubre2021_Metadata_BDMexCov2_FINAL.xlsx into the new file (manually)
# Most batches can be identified using the fastas in two folders
vigilancia="2021-10-28_CoViGen-Mex/Secuencias"
for i in $(ls $vigilancia/*.fasta);do base=$(basename $i);echo 'grep -a "^>" '$i'|grep -av "NC"|sed "s/^>//" >'$in_dir/${base%.fasta}'.list';done|bash # a sed was added to remove accents in "México"
# The following will replace the old sets with new corrected sets
vigilancia="2021-10-28_CoViGen-Mex/Secuencias/Secuencias_CORRECCION"
for i in $(ls $vigilancia/*.fasta);do base=$(basename $i);echo 'grep -a "^>" '$i'|grep -av "NC"|sed "s/^>//" >'$in_dir/${base%.fasta}'.list';done|bash
# files 21Abril2021_EpiCoV_Lote8_INER.list and 21Abril2021_EpiCoV_Lote8_INER_seq.list were repeated, thus it will be replaced in the following commands
# We'll be changing the original names, so we save a backup and create a template to work manually
ls $in_dir/*.list >$in_dir/original_names.txt; cp $in_dir/original_names.txt $in_dir/new_names.txt
# The new_names.txt was manually edited
paste $in_dir/original_names.txt $in_dir/new_names.txt|sed 's/^/mv /'|bash
# Yet, as of 2021-10-29, some files are missing, including batches 18, 20, 21, 23, and 24
# These were instead obtained from tables in 2021-10-28_CoViGen-Mex/Metadata/LinajesConSINAVE/
# Create empty containers for the missing lists:
touch $in_dir/Lote18_2021-07-28_IBT.list
touch $in_dir/Lote20_2021-08-11_LANGEBIO.list
touch $in_dir/Lote21_2021-08-25_IBT.list
touch $in_dir/Lote23_2021-09-09_LANGEBIO.list
touch $in_dir/Lote24_2021-09-22_IBT.list
# Files were filled with the Virus Name information and "non-complete" (NC) items were removed manually (there were fewer than 10 in each)
# NOTE: Batch 25 was included in the CONACYT report on October 11th (these reports [@ 2021-10-28_BD_CONACYT/] include non-CoViGen samples as well). Batch 26 should come around Nov 1st for our lab data and latter in the CONACYT report.
# There was an extra problem with Lote13 because the fastas (368) and the metadata in LinajesConSINAVE (500) did not match. Mauricio explained to me that this difference was due to a late update by INER and CIAD 11 and 12 that were included later. We now have fastas for them @ 2021-10-28_CoViGen-Mex/Secuencias/Secuencias_CORRECCION/Extras_mauricio_Lote13/.
# Test some of these
mauricio="2021-10-28_CoViGen-Mex/Secuencias/Secuencias_CORRECCION/Extras_mauricio_Lote13"
for i in $(ls $mauricio/*.fasta);do base=$(basename $i);echo 'grep "^>" '$i'|sed "s/^>//" >'$mauricio/${base%.fasta}.list'';done|bash
grep -wFf $mauricio/CIAD11.list $in_dir/Lote13_2021-06-09_IBT.list|wc -l # Test those from CIAD 11
# 34
grep -wFf $mauricio/CIAD12.list $in_dir/Lote13_2021-06-09_IBT.list|wc -l # Test those from CIAD 12
# 27
grep -wFf $mauricio/gisaid_hcov-19_2021_06_24_21_INER-IBT.list $in_dir/Lote13_2021-06-09_IBT.list|wc -l # Test those from INER
# 71

# A list of the extra sequences from CIAD and INER in batch 13 were stored @ 01_GISAID_Tables/2021-10-28_Vigilancia/Correct_but_replaced/Lote13_LinajesConSINAVE_not_in_fastas, along with the original IBT-only set
# TROUBLESHOOTING: Also, there may have been some problem with the carriage return in some files created from the fastas as I had to trim the EOL character in order for grep to work [e.g. grep -vwFf <(sed "s/.$//" mini.list) grande.txt].

# Now, I'll create subsets of the metadata table using the expected genomes per batch
out_dir="$in_dir/01_new_metadata_tables"
mkdir $out_dir $out_dir/failed
# Search for the items found in each batch
# TROUBLESHOOTING NOTE: We detected there are several mislabeled items. This includes (but is not restricted to):
# are accents in Lote10
# some are _LANGEBIO instead of -LANGEBIO in Lote16 (but there are also _LANGEBIO in GISAID so no single solution is available)
# _INER instead of -INER in Lote17
# There is no single fix that can deal with all exeptions so it might be faster to edit the corresponding list files instead
# To cope with most simple errors, all - and _ are taken as _ now and Lote10 was fixed manually (remove accents)
# Several failed, we'll try to track those as well
# Most were issues with swapped - and _.
# This was fixed by creating middle files to use all as "-" but keeping track of the origina lnames.
# Still, there was a second issue with some files having a bad carriage return.
# This was addressed with a sed removing the last character but I had to patch those that didn't have them.
# NOTE: this requires some changes where 2022 may be included instead of 2021 as this exception is not considered
in_dir="01_GISAID_Tables/2021-10-28_Vigilancia"
out_dir="$in_dir/01_new_metadata_tables"
# for i in $(ls $in_dir/Lote*.list);do base=$(basename $i);echo 'paste <(sed -e "s/.$//" -e "s/\/202$/\/2021/" '$i'|cut -f 1) <(sed -e "s/.$//" -e "s/\/202$/\/2021/" -e "s/_/-/g" '$i'|cut -f 1)|grep -v "^$" >temp; grep -Ff <(cut -f 2 temp) <(sed -e "s/.$//" -e "s/_/-/g" '$in_dir'/2021-10-25_Metadata.tsv)|sort|sed "s/EPI-ISL-/EPI_ISL_/g" >temp2;grep -Ff <(cut -f 1 temp2) temp|sort -k 2 >temp3;paste <(cut -f 1 temp3) <(cut -f 2- temp2) >temp4; grep -wFf <(cut -f 2 temp4) '$in_dir'/2021-10-25_Metadata.tsv|sort|paste <(cut -f 1 temp4) - >'$out_dir/${base%.list}.tsv'; grep -vFf <(cut -f 1 temp2) temp|sort -k 2|cut -f 1 >'$out_dir/failed/${base%.list}.failed';rm temp temp2 temp3 temp4';done|bash
# There were 739 missing items now
# Update 2021-11-03: This was updated to use a more recent table (2021-11-01_Metadata_GISAID.tsv)
for i in $(ls $in_dir/Lote*.list);do base=$(basename $i);echo 'paste <(sed -e "s/.$//" -e "s/\/202$/\/2021/" '$i'|cut -f 1) <(sed -e "s/.$//" -e "s/\/202$/\/2021/" -e "s/_/-/g" '$i'|cut -f 1)|grep -v "^$" >temp; grep -Ff <(cut -f 2 temp) <(sed -e "s/.$//" -e "s/_/-/g" '$in_dir'/2021-11-01_Metadata_GISAID.tsv)|sort|sed "s/EPI-ISL-/EPI_ISL_/g" >temp2;grep -Ff <(cut -f 1 temp2) temp|sort -k 2 >temp3;paste <(cut -f 1 temp3) <(cut -f 2- temp2) >temp4; grep -wFf <(cut -f 2 temp4) '$in_dir'/2021-11-01_Metadata_GISAID.tsv|sort|paste <(cut -f 1 temp4) - >'$out_dir/${base%.list}.tsv'; grep -vFf <(cut -f 1 temp2) temp|sort -k 2|cut -f 1 >'$out_dir/failed/${base%.list}.failed';rm temp temp2 temp3 temp4';done|bash
# Missing items are still the same 739

#DEPRECATED: These two commands were the previously used approach, first creating matching tables, then the fail report.
# for i in $(ls $in_dir/Lote*.list);do base=$(basename $i);echo 'grep -Ff <(sed "s/.$//" '$i'|cut -f 1) '$in_dir/2021-10-25_Metadata.tsv' >'$out_dir/${base%.list}.tsv'';done|bash
#for i in $(ls $out_dir/Lote*.tsv);do base=$(basename $i);echo 'grep -vwFf <(cut -f 1 '$i') '$in_dir/${base%.tsv}.list' >'$out_dir/failed/${base%.tsv}.failed'';done|bash
# It seems some items may have been removed from the current GISAID metadata (2021-10-29) whathever the reason may be due to CoViGen removing them or GISAID itself.
# Update 2021-11-03: now collate all tables
# Create a single table for the surveillance project:
outfile="01_GISAID_Tables/2021-10-28_All_batches_Full_GISDAID_metadata.tsv"
printf "Lote\tCoViGen - Nombre Interno\tGISAID - Actualizacion " >$outfile
head -n 1 $in_dir/2021-11-01_Metadata_GISAID.tsv >>$outfile
for i in $(ls $out_dir/Lote*.tsv);do base=$(basename $i);echo 'awk '\''BEGIN{FS=OFS="\t"}{print "'${base%.tsv}'" OFS $0}'\'' '$i'';done|bash >>$outfile
# And a small subset containing the id, date, lineage and State
printf "Accession ID\tCollection date\tLineage\tLocation\n" >01_GISAID_Tables/2021-10-28_All_batches_for_areaplots.tsv
paste <(cut -f 4,7,26 $outfile) <(cut -f 9 $outfile)|tail +2|sed  -e 's/\t\([1-9]\)\//\t0\1\//' -e 's/\/\([1-9]\)\//\/0\1\//'|sed 's/\t\(\w\w\)\/\(\w\w\)\//\t\2\/\1\//' >>01_GISAID_Tables/2021-10-28_All_batches_for_areaplots.tsv

# Finally, create a single table of all items from both the P1000 and the surveillance projects
outfile="01_GISAID_Tables/2021-10-28_All_items_for_areaplots.tsv"
printf "Accession ID\tCollection date\tLineage\tLocation\n" >$outfile
tail +2 01_GISAID_Tables/2021-10-28_First_year_for_areaplots.tsv >>$outfile
tail +2 01_GISAID_Tables/2021-10-28_All_batches_for_areaplots.tsv >>$outfile

# NOTE: One item had bad State information " Mexico City" (with that space), and was fixed manually

# ### 02: SPLIT AND PLOT TABLES BY DATE ###
# Update 2021-11-03: We now have a table for first-year sequences with location and lineage and one for 25 batches in the surveillance project. We want to plot per specific times
# Each table has 4 items: Accession ID	Collection date	Lineage	Location
# I first created a script to split these tables by date:
mkdir 02_Split_dates
Rscript extract_dates.R 01_GISAID_Tables/2021-10-28_First_year_for_areaplots.tsv 2020-03-18 2020-06-30 02_Split_dates/First_year
Rscript extract_dates.R 01_GISAID_Tables/2021-10-28_First_year_for_areaplots.tsv 2020-07-01 2020-10-31 02_Split_dates/First_year
Rscript extract_dates.R 01_GISAID_Tables/2021-10-28_First_year_for_areaplots.tsv 2020-11-01 2021-02-03 02_Split_dates/First_year
Rscript extract_dates.R 01_GISAID_Tables/2021-10-28_First_year_for_areaplots.tsv 2020-03-18 2021-02-03 02_Split_dates/First_year
Rscript extract_dates.R 01_GISAID_Tables/2021-10-28_All_batches_for_areaplots.tsv 2021-02-03 2021-05-31 02_Split_dates/Second_year
Rscript extract_dates.R 01_GISAID_Tables/2021-10-28_All_batches_for_areaplots.tsv 2021-06-01 2021-09-30 02_Split_dates/Second_year

# For this, I modified some of Mauricio's scripts (also from GAL) to process any number of tables
mkdir 03_Lineages

# Update 2021-11-25: We now have most of the fastq sequences, which will be processed in teopanzolco (cluster). Commands can be found @ 00_docs/cluster_commands.sh
# For renaming CIAD files, we first requiere the corresponding xref tables @ 00_docs/CIAD_xref_lotes_muestras.tsv, which have the following items: Sample  Genome  CoViGen CIAD
sed -e 's/hCoV-19\/Mexico\///' -e 's/\/20/-20/' -e 's/\(\t[A-Z]*\)_/\1-/' 00_docs/CIAD_xref_lotes_muestras.tsv |awk '{print $0 "\t""Lote"$3"_CIAD"$4}'|cut -f 1,2,5 >00_docs/CIAD_xref_for_rename.tsv

