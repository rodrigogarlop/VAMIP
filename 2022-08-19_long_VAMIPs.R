# Started: 2022-08-19
# by Rodrigo García-López for Carlos Arias's Viromics Group at IBt, UNAM, Cuernavaca, Mexico as part of the CoViGen-Mex SARS-CoV-2 survillance in Mexico.
# Under GNU GPLv3 license
# Disclamer: This script was written for use with specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.

# The script is a modified version of the analysis scripts ML_sequences.R and 2022-06-26_Analyze_vamip_contingency.R, and is intended to analyze VAMIPs in longitudinal data from a contingency table containing genomes (columns) x particular mutations (rows). The former include some metadata in their names, which can be parsed previously.

# ### LOAD LIBRARIES AND FUNCTIONS ###
library("pheatmap") # This is used for heatmaps (better than R base's
library("stringr") # This is required to extract the metadata easily
save_pheatmap_pdf <- function(x, filename, width, height) { # This is a printing function for pheatmaps with custom filename, width and height. X is the input heatmap made with the pheatmap package.
  pdf(filename,width,height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
strip <- function(string){ # Takes a string object and strips rare characters
	string <- iconv(string,from="UTF-8",to="ASCII//TRANSLIT") # Remove accents and rare chars
# 	string <- tolower(string) # change to lower
	string <- gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", string, perl=TRUE) # Remove multiple spaces
	string <- gsub("\\+","",string) # remove + characters
	return(string)
}
find_items <- function(queryStrings, dataf, colname){ # IMPORTANT. This is used to avoid partial index matching. Use it to hash a table and use rownames to match a vector. colname may be numeric
	searchVect <- dataf[,colname]
	names(searchVect) <- rownames(dataf)
	out <- searchVect[unlist(queryStrings)]
	return(out)
}
group_map <- function(vect) { # Get a mutation table (columns = genomes; rows = mutations), and a vector of the same rowsize. Then use the vector to define how groups should be splitted and build a list of df with the table observations per group (a map of positions)
	group_names <- levels(as.factor(vect))
	group_positions <- sapply(group_names,function(x){grep(paste0("^",strip(x),"$"),strip(vect))})
	return(group_positions)
}
group_prevalence <- function(gmap,inmat) { # Get the group map from the corresponding function and get the prevalence per row item. This is calculated based on the total NAs.
	nam <- names(gmap) # Save names for indices
	out_raw <- matrix(NA, nrow=nrow(inmat), ncol=length(nam)); rownames(out_raw) <- rownames(inmat); colnames(out_raw) <- nam # Create a void matrix container
	out_rel <- out_raw # make them two
	for(i in nam){ # now, for each table
# 		print(i)
		single_grp <- inmat[,gmap[[i]]] # First, subset the input table using the group map for a single group
		out_raw[,i] <- apply(single_grp,1,function(x){length(x[!is.na(x)])})  # and calculate the total items
		out_rel[,i] <- out_raw[,i]*100/ncol(single_grp)
	}
	out <- list(out_raw,out_rel)
	names(out) <- c("raw","rel")
	return(out)
}

group_freq_mean <- function(gmap,inmat) { # Get the group map from the corresponding function and get the prevalence per row item. This is calculated based on the total non-NAs.
	nam <- names(gmap) # Save names for indices
	out_raw <- matrix(NA, nrow=nrow(inmat), ncol=length(nam)); rownames(out_raw) <- rownames(inmat); colnames(out_raw) <- nam # Create a void matrix container
	for(i in nam){ # now, for each table
# 		print(i)
		single_grp <- inmat[,gmap[[i]]] # First, subset the input table using the group map for a single group
		temp <- apply(single_grp, 1, function(x){mean(x[which(!is.na(x))])})  # and calculate the total items
		temp[is.nan(temp)] <- NA
		out_raw[,i] <- temp
# 		out_rel[,i] <- out_raw[,i]*100/ncol(single_grp)
	}
	return(out_raw)
}
group_freq_centile <- function(gmap,inmat,cent=50) { # Get the group map from the corresponding function and get the prevalence per row item. This is calculated based on the total non-NAs.
	cent <- cent*0.01
	nam <- names(gmap) # Save names for indices
	out_raw <- matrix(NA, nrow=nrow(inmat), ncol=length(nam)); rownames(out_raw) <- rownames(inmat); colnames(out_raw) <- nam # Create a void matrix container
	for(i in nam){ # now, for each table
# 		print(i)
		single_grp <- inmat[,gmap[[i]]] # First, subset the input table using the group map for a single group
		out_raw[,i] <- apply(single_grp, 1, function(x){quantile(x[which(!is.na(x))],cent)})  # and calculate the quantile (centile 0-100)
# 		temp[is.nan(temp)] <- NA
# 		out_raw[,i] <- temp
# 		out_rel[,i] <- out_raw[,i]*100/ncol(single_grp)
	}
	return(out_raw)
}

# ### LOAD INPUTS ###
df <- read.table("Alt_variant_Freq.tsv", sep="\t",header=T, skip=0, comment.char='',quote="",fill=F,row.names=1, check.names=FALSE) # Load input table
# Each column was named as follows: "2021-03-04|L004_202101083013|B.1.1.222(20B_Others)|046|M|Puebla" and they have been sorted before
dim(df)
# [1]  7340 17927
df <- df[,order(colnames(df))] # Sort the columns (they start with the date)
sum(!is.na(df))
# [1] 1034529
df[df==1e-06] <- NA # First off, remove items that are below the cutoff (this was set when the input table was created to make them 1e-06)
sum(!is.na(df)) # we test the total number of items
# [1] 872199
# PREPARE METADATA
# Process the name to extract basic info and the Folio ID
# names are processed into a table: Date|Folio ID|Pango(clade)|age|gender|State
# Example: "2021-03-04|L004_202101083013|B.1.1.222(20B_Others)|046|M|Puebla"
minimetatab <- t(as.data.frame(str_split(colnames(df), "\\|"))) # Extract a mini metadata table from the names (including
rownames(minimetatab) <- NULL # row names are removed
meta <- read.table("07_metadata/Metadata_all.tsv", header=T, sep ='\t',stringsAsFactors = FALSE, check.names=F, fill=T, quote="") # Load metadata with id folio for xref. IDs should be unique
rownames(meta) <- meta[,"ID Folio"]
dim(meta)
# [1] 19024    48
# Now, explore basic dataset composition per genome
# Extract dates from the column names (each genome name includes the sample collection date)
# The rest of the metadata needs to be extracted from the external metadata table using the ID Folio for xref
meta <- meta[unique(minimetatab[,2]),]
dim(meta)
# [1] 17927		48
write.table(meta, "07_metadata/Metadata_matched_only.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
# This table may be a larger than the mutations table (columns) but only matching items will be used
tot_per_mut <- apply(df,1,function(x){sum(!is.na(x))})

# ### MAIN ###
# Create two copies to compare only vamips (<0.05) and only major variants
df_vamip <- df; df_vamip[df_vamip>0.5] <- NA # this will hold vamips (>0.5 are now NAs)
df_novamip <- df; df_novamip[df_novamip<=0.5] <- NA # this will hold vamips (>0.5 are now NAs)

# Some metadata are already found in the genome's name and will be processed directly, others must be extracted from the external meta table
genomes <- minimetatab[,2] # Extract the genome names (this should include their corresponding batches
# DATE
all_Days <- as.Date(minimetatab[,1]) # First, the day
all_Months <- format(all_Days,"%Y-%m") # Then, the month
months_distro <- all_only_Month <- format(all_Days,"%b") # we also get the month only in an abbreviated version
# sapply(names(table(all_Months)),function(x){grep(x,all_Months)})
# For the week, it's more complicated, as it depends on the epidemiological weeks
week <- data.frame("week"=c(paste0("20W",sprintf('%0.2d', rep(1:52,each=7))),paste0("21W",sprintf('%0.2d', rep(1:52,each=7))),paste0("22W",sprintf('%0.2d', rep(1:52,each=7)))))
# Now rename the rows to use them as a hash (dictionary)
rownames(week) <- seq(as.Date("2020/01/05"), as.Date("2022/12/31"), by="day")[1:nrow(week)]
# write.table(week, "week_calendar.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
# Now used the newly created dictionary to append the corresponding week
all_Weeks <- as.character(week[as.character(all_Days),1]) # match days with their weeks
week_start <- sapply(names(table(all_Weeks)),function(x){grep(x,all_Weeks)[1]}) # We want to add the month when the first day of each week starts
temp <- table(all_Weeks) # to append the month to the week we use a temporary vector
names(temp) <- paste(names(table(all_Weeks)),months_distro[week_start]) # this can only be done because we have the same order
weeks_month_st <- cbind("Week"=names(table(all_Weeks)),"Month (starting day)"=months_distro[week_start])
write.table(weeks_month_st, "07_metadata/weeks2_Month_start.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
head(cbind("Full"=colnames(df),"Genome"=genomes, "Date"=as.character(all_Days), "Month"=all_Months, "Char_Month"=months_distro, "Week"=all_Weeks))

# CLADES & VARIANTS
all_lineages <- find_items(minimetatab[,2],meta, "New_Lineage") # Start with the pango lineages
all_mainvar <- find_items(minimetatab[,2],meta, "Variants") # Next, the main variants (VOCs and 519)
all_clade <- find_items(minimetatab[,2],meta, "Clado Nexstrain") # Next, the main clades
# AGE
all_age <- find_items(minimetatab[,2],meta, "Edad") # Start with the raw age
all_age_vac <- find_items(minimetatab[,2],meta, "Age_vac") # Now with vaccination group
# GENDER
all_sex <- find_items(minimetatab[,2],meta, "Genero")
# PATIENT STATUS
all_status <- find_items(minimetatab[,2],meta, "Estatus del paciente")
# LOCATION
all_region7 <- find_items(minimetatab[,2],meta, "Region_7")
# CT
all_ctrdrp <- round(as.numeric(find_items(minimetatab[,2],meta, "RdRP")))
all_cte <- round(as.numeric(find_items(minimetatab[,2],meta, "E")))
# Additional items
xtra_month <- find_items(minimetatab[,2],meta, "Month")
xtra_mut_Nt <- find_items(minimetatab[,2],meta, "Mutaciones Nucleotidos")
xtra_mut_AA <- find_items(minimetatab[,2],meta, "Mutaciones Aminoacido")

out <- cbind("Full"=colnames(df),"Genome"=genomes, "Date"=as.character(all_Days), "Month"=xtra_month, "Char_Month"=months_distro, "Week"=all_Weeks, "Lineage"=all_lineages, "Variant"=all_mainvar, "Clade"=all_clade, "Age"=all_age, "Age_vac"=all_age_vac, "Sex"=all_sex, "Status"=all_status, "Region"=all_region7, "Ct_RdRp"=all_ctrdrp, "Ct_E"=all_cte)
write.table(out,"07_metadata/metadata_VAMIP_set-Jan2021-Mar2022.tsv",sep='\t', row.names=FALSE, col.names=TRUE)

# BASIC PLOTS
pdf("08_Analyses/Total_genomes-month.pdf")
	barplot(las=2,table(all_Months),border=NA,col="cornflowerblue",yaxt='n',ylim=c(0,2000))
	axis(2,las=1, at=seq(0,2000,100))
dev.off()
pdf("08_Analyses/Total_genomes-week.pdf",width=14)
	barplot(las=2,temp,border=NA,col="cornflowerblue",yaxt='n',ylim=c(0,700))
	axis(2,las=1, at=seq(0,2000,50))
dev.off()
pdf("08_Analyses/Total_genomes-variants.pdf")
	barplot(las=2,table(all_mainvar),border=NA,col="coral1",yaxt='n',ylim=c(0,9000))
	axis(2,las=1, at=seq(0,9000,1000))
dev.off()
pdf("08_Analyses/Total_genomes-clades.pdf",width=14)
	barplot(las=2,table(all_clade),border=NA,col="coral1",yaxt='n',ylim=c(0,9000))
	axis(2,las=1, at=seq(0,9000,1000))
dev.off()
temp <- table(all_age)
temp <- temp[-length(temp)]
temp <- temp[as.character(0:101)]
pdf("08_Analyses/Total_genomes-age.pdf",width=14)
	barplot(las=2,temp,border=NA,col="forestgreen",yaxt='n',ylim=c(0,400))
	axis(2,las=1, at=seq(0,400,50))
dev.off()
pdf("08_Analyses/Total_genomes-age_vac.pdf")
	barplot(las=2,table(all_age_vac),border=NA,col="forestgreen",yaxt='n',ylim=c(0,5000))
	axis(2,las=1, at=seq(0,5000,500))
dev.off()
pdf("08_Analyses/Total_genomes-sex.pdf",width=4)
	barplot(las=2,table(all_sex),border=NA,col="purple",yaxt='n',ylim=c(0,9000))
	axis(2,las=1, at=seq(0,9000,1000))
dev.off()
pdf("08_Analyses/Total_genomes-region7.pdf",width=7)
	barplot(las=2,table(all_region7),border=NA,col="firebrick",yaxt='n',ylim=c(0,4000))
	axis(2,las=1, at=seq(0,4000,500))
dev.off()
pdf("08_Analyses/Total_genomes-status.pdf",width=4)
	barplot(las=2,table(all_status),border=NA,col="gold2",yaxt='n',ylim=c(0,11000))
	axis(2,las=1, at=seq(0,11000,1000))
dev.off()
pdf("08_Analyses/Total_genomes-ctRdRp.pdf",width=14)
	barplot(las=2,table(all_ctrdrp),border=NA,col="aquamarine3",yaxt='n',ylim=c(0,2000))
	axis(2,las=1, at=seq(0,2000,250))
dev.off()
pdf("08_Analyses/Total_genomes-ctE.pdf",width=14)
	barplot(las=2,table(all_cte),border=NA,col="aquamarine3",yaxt='n',ylim=c(0,2750))
	axis(2,las=1, at=seq(0,2750,250))
dev.off()

# COMPARE PREVALENCE
# First, we want to compare a raw view of all totals, with and without the vamips:
compare_all <- cbind("All"=rowSums(group_prevalence(group_map(all_sex),df)$raw), "only_VAMIPs"=rowSums(group_prevalence(group_map(all_sex),df_vamip)$raw),"only_Consensus"=rowSums(group_prevalence(group_map(all_sex),df_novamip)$raw)) # Calculate the total genomes having each mutations, and the number of times they are seen as vamips and as as fixated major mutations
compare_all_total <- cbind(compare_all, compare_all/ncol(df)*100) # also include % of the values above
colnames(compare_all_total) <- c(colnames(compare_all_total)[1:3],paste(colnames(compare_all_total)[1:3],"%")) # And append them
write.table(compare_all_total, "08_Analyses/VAMIP_vs_Consensus.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)
pdf("08_Analyses/Histogram_mutations.pdf") # Summarize with histograms
hist(log10(compare_all_total[,1]),main="Histograma: Genomas con cada mutación (todas)",las=1, col="turquoise3",border="white", xaxt='n', xlab='Cuántos genomas tienen esa mutación (10^log10)', ylab='Mutaciones puntuales', breaks=50,xlim=range(log10(compare_all_total[,1])), ylim=c(0,1200))
temp <- c(5,2,10,20,50,100,200,500,1000,2000,5000,10000,max(compare_all_total[,1]))
axis(1,at=log10(temp),labels=(temp),las=2)
hist(log10(compare_all_total[,2]),main="Histograma: Genomas con cada mutación (sólo VAMIP)",las=1, col="cornflowerblue",border="white", xaxt='n', xlab='Cuántos genomas tienen esa mutación (10^log10)', ylab='Mutaciones puntuales', breaks=50,xlim=range(log10(compare_all_total[,1])), ylim=c(0,1200))
temp <- c(5,2,10,20,50,100,200,500,1000,2000,5000,10000,max(compare_all_total[,1]))
axis(1,at=log10(temp),labels=(temp),las=2)
hist(log10(compare_all_total[,3]),main="Histograma: Genomas con cada mutación (sólo en Consensos)",las=1, col="coral1",border="white", xaxt='n', xlab='Cuántos genomas tienen esa mutación (10^log10)', ylab='Mutaciones puntuales', breaks=50,xlim=range(log10(compare_all_total[,1])), ylim=c(0,1200))
temp <- c(5,2,10,20,50,100,200,500,1000,2000,5000,10000,max(compare_all_total[,1]))
axis(1,at=log10(temp),labels=(temp),las=2)
dev.off()
data <- rev(c("All Mut"=nrow(df),"All VAMIP"=sum(compare_all_total[,2]>0),"Consensus and VAMIP"=sum((compare_all_total[,3]>0)*(compare_all_total[,2]>0)),"All in Consensus"=sum(compare_all_total[,3]>0),"Only VAMIP"=sum(compare_all_total[,2]==compare_all_total[,1]),"Only in Consensus"=sum(compare_all_total[,3]==compare_all_total[,1])))
pdf("08_Analyses/VAMIP_consensus_compare.pdf")
barplot(las=2,data,horiz=T, col=c("cornflowerblue","coral1","coral1","cornflowerblue","coral1","chartreuse3"), border=NA)
dev.off()
# Create a screening of all mutations by month prevalence
g_month <- group_prevalence(group_map(all_Months),df) # create a summary how many items we have per month
g_month_vamip <- group_prevalence(group_map(all_Months),df_vamip) # also a create a subset for only those who are vamips
g_month_novamip <- group_prevalence(group_map(all_Months),df_novamip) # and one for those who were not vamips (major fixated mutations in the consensus)
sub <- head(g_month$rel[names(sort(rowSums(g_month$rel),decreasing=TRUE)),],2000) # get the 2000 most common mutations only
sub[sub==0] <- NA
sub_vamip <- head(g_month_vamip$rel[names(sort(rowSums(g_month_vamip$rel),decreasing=TRUE)),],2000)
sub_vamip[sub_vamip==0] <- NA
sub_novamip <- head(g_month_novamip$rel[names(sort(rowSums(g_month_novamip$rel),decreasing=TRUE)),],2000)
# with these, we can get an overhead view of all the most common mutations
sub_novamip[sub_novamip==0] <- NA
out_heatmap <- pheatmap(sub, cluster_cols=FALSE, cluster_rows=FALSE)
save_pheatmap_pdf(out_heatmap, "08_Analyses/Month_heatmap_2000TopMut.pdf", 9, 230)
out_heatmap <- pheatmap(sub_vamip, cluster_cols=FALSE, cluster_rows=FALSE)
save_pheatmap_pdf(out_heatmap, "08_Analyses/Month_vamip_heatmap_2000TopMut.pdf", 9, 230)
out_heatmap <- pheatmap(sub_novamip, cluster_cols=FALSE, cluster_rows=FALSE)
save_pheatmap_pdf(out_heatmap, "08_Analyses/Month_novamip_heatmap_2000TopMut.pdf", 9, 230)
# Now, use the 2000 most abundant items in sub_novamip to extract them from the novamip set
sub_novamip_same_as_sub_vamip <- g_month_novamip$rel[rownames(sub_vamip),]
sub_novamip_same_as_sub_vamip[sub_novamip_same_as_sub_vamip==0] <- NA
out_heatmap <- pheatmap(sub_novamip_same_as_sub_vamip, cluster_cols=FALSE, cluster_rows=FALSE)
save_pheatmap_pdf(out_heatmap, "08_Analyses/Month_novamip_same_set_as_vamip_heatmap_2000TopMut.pdf", 9, 230)
# Compare those with vamips and novamips:
rownames(sub_novamip_same_as_sub_vamip) <- paste(rownames(sub_novamip_same_as_sub_vamip),"X novamip") # First, rename both sets
rownames(sub_novamip_same_as_sub_vamip) <- paste(sprintf("%04d", 1:2000),rownames(sub_novamip_same_as_sub_vamip))
rownames(sub_vamip) <- paste(rownames(sub_vamip),"VAMIP")
rownames(sub_vamip) <- paste(sprintf("%04d", 1:2000),rownames(sub_vamip))
compare <- rbind(sub_vamip,sub_novamip_same_as_sub_vamip)
compare <- compare[order(rownames(compare)),]
out_heatmap <- pheatmap(compare, cluster_cols=FALSE, cluster_rows=FALSE)
save_pheatmap_pdf(out_heatmap, "08_Analyses/Month_heatmap_compare_2000TopMut.pdf", 9, 460)

freq_mean_month <- group_freq_mean(group_map(all_Months),df) # Create a colection of means (frequency) per  month
freq_mean_month <- freq_mean_month[apply(freq_mean_month, 1, function(x){sum(!is.na(x))})>=5,] # keep only those appearing in 5 months (these may not be contiguous)
out_heatmap <- pheatmap(freq_mean_month, cluster_cols=FALSE, cluster_rows=FALSE) #
save_pheatmap_pdf(out_heatmap, "08_Analyses/freq_mean_month.pdf", 9, 550)
# Repeat but with row clustering
sub <- freq_mean_month
sub[is.na(sub)]=-10
# sub <- sub[,colSums(sub)>0] # Filters columns not having any
dim(sub)
# [1] 4827   15
clust_row <- hclust(dist(sub), method = "mcquitty") # These versions create clusters (WPGMA)
clust_col <- hclust(dist(t(sub)), method = "mcquitty")
out_heatmap <- pheatmap(freq_mean_month, cluster_cols=FALSE, cluster_rows=clust_row) # Repeat but use row clustering
save_pheatmap_pdf(out_heatmap, "08_Analyses/freq_mean_month_crows.pdf", 9, 550)

freq_q50_month <- group_freq_centile(group_map(all_Months), df, 50) # Create a colection of centiles (of the frequency) per  month (no considering NAs). In this case, we use quantile 50.
freq_q50_month <- freq_q50_month[apply(freq_q50_month, 1, function(x){sum(!is.na(x))})>=5,] # keep only those appearing in 5 months (these may not be contiguous)
out_heatmap <- pheatmap(freq_q50_month, cluster_cols=FALSE, cluster_rows=FALSE) #
save_pheatmap_pdf(out_heatmap, "08_Analyses/freq_q50_month.pdf", 9, 550)
































# Same thing, but use depths instead
df2 <- read.table("Alt_variant_TotalDepth.tsv", sep="\t",header=T, skip=0, comment.char='',quote="",fill=F,row.names=1, check.names=FALSE) # Load input table
dim(df2)
# [1]  7340 17927
df2 <- df2[,order(colnames(df2))]
minimetatab <- t(as.data.frame(str_split(colnames(df2), "\\|")))
rownames(minimetatab) <- NULL # row names are removed
genomes <- minimetatab[,2]
xrefnames <- read.table("07_metadata/IMPORT_MLReferenceTable.tsv", header=T, sep ='\t', stringsAsFactors = FALSE, row.names=1, check.names=F) # Load unique mutation information
xrefnames <- xrefnames[xrefnames[,3]!="",1:3] # Remove those with no folio
rownames(xrefnames) <- xrefnames[,"Folio"] # and use folios are keys
xrefnames["newname"] <- paste(xrefnames[,1],xrefnames[,2],sep='-')
temp <- xrefnames[minimetatab[,2],4]
colnames(df2)[!is.na(temp)] <- paste0(colnames(df2)[!is.na(temp)],"||",temp[!is.na(temp)])

# Now select only those flagged by selene as multilineages
multlin <- df2[,grep("\\|ML",colnames(df2))]
multlin <- multlin[keep_mutations,]
dim(multlin)
# [1] 330  20
write.table(multlin,"ML_Total_depth.tsv", quote=FALSE, sep='\t', row.names=TRUE, col.names=NA)

ML_alt <- read.table("ML_alt_freq.tsv", sep="\t",header=T, skip=0, comment.char='',quote="",fill=F,row.names=1, check.names=FALSE) # Load input table
sub <- ML_alt
sub[is.na(sub)]=0
sub <- sub[,colSums(sub)>0] # Filters those not having any
dim(sub)
# [1] 330  20
clust_row <- hclust(dist(sub), method = "mcquitty") # These versions create clusters (WPGMA)
clust_col <- hclust(dist(t(sub)), method = "mcquitty")
# Plot heatmap
out_heatmap <- pheatmap(t(ML_alt), cluster_cols=FALSE, cluster_rows=FALSE)
# save_pheatmap_pdf(out_heatmap, "ML_heatmap_all_mut.pdf", 50, 10)
groups <- locations <- as.numeric(sub(":.*","",rownames(ML_alt)))
groups[groups>0] <- "00 Genomic"
groups[locations<266] <- "01 5-UTR"
groups[locations>=266 & locations<=13467] <- "02 ORF1a"
groups[locations>=13468 & locations<=21555] <- "03 ORF1b" # 13442-13467 have no assignation
groups[locations>=21555 & locations<=25384] <- "04 S"
groups[locations>=25393 & locations<=26220] <- "05 ORF3a"
groups[locations>=26245 & locations<=26472] <- "06 E"
groups[locations>=26523 & locations<=27191] <- "07 M"
groups[locations>=27202 & locations<=27387] <- "08 ORF6"
groups[locations>=27394 & locations<=27759] <- "09 ORF7a"
groups[locations>=27756 & locations<=27887] <- "10 ORF7b"
groups[locations>=27894 & locations<=28259] <- "11 ORF8"
groups[locations>=28274 & locations<=29533] <- "12 N"
groups[locations>=28283 & locations<=28573] <- "13 ORF9b"
groups[locations>=29558 & locations<=29674] <- "14 ORF10"
groups[locations>=29675 & locations<=29903] <- "15 3-UTR"

custom_colours = list(groups = c("00 Genomic"="gray", "01 5-UTR"="cornflowerblue", "02 ORF1a"="coral1","03 ORF1b"="turquoise3", "04 S"="chartreuse2", "05 ORF3a"="purple", "06 E"="firebrick", "07 M"="gold2", "08 ORF6"="hotpink", "09 ORF7a"="forestgreen", "10 ORF7b"="darkblue", "11 ORF8"="darkslategray", "12 N"="darkorange2", "13 ORF9b"="turquoise3", "14 ORF10"="brown2", "15 3-UTR"="mediumpurple"))
row_labels <- data.frame(groups)
rownames(row_labels) <- rownames(ML_alt)
# Plot heatmap
out_heatmap <- pheatmap(t(ML_alt), cluster_cols=FALSE, cluster_rows=FALSE, annotation_colors = custom_colours, border_col="white", annotation_col = row_labels)
save_pheatmap_pdf(out_heatmap, "ML_heatmap_vamip_freq_mut_noClust.pdf", 50, 9)
out_heatmap <- pheatmap(t(ML_alt), cluster_cols=FALSE, cluster_rows=clust_col, annotation_colors = custom_colours, border_col="white", annotation_col = row_labels)
save_pheatmap_pdf(out_heatmap, "ML_heatmap_vamip_freq_mut_clust_genome.pdf", 50, 9)
out_heatmap <- pheatmap(t(ML_alt), cluster_cols=clust_row, cluster_rows=clust_col, annotation_colors = custom_colours, border_col="white", annotation_col = row_labels)
save_pheatmap_pdf(out_heatmap, "ML_heatmap_vamip_freq_mut_clust_both.pdf", 50, 9)

ML_dp <- read.table("ML_Total_depth.tsv", sep="\t",header=T, skip=0, comment.char='',quote="",fill=F,row.names=1, check.names=FALSE) # Load input table
sub <- ML_dp
sub[is.na(sub)]=0
sub <- sub[,colSums(sub)>0] # Filters those not having any
dim(sub)
# [1] 330  20
clust_row <- hclust(dist(sub), method = "mcquitty") # These versions create clusters
clust_col <- hclust(dist(t(sub)), method = "mcquitty")
# Plot heatmap
out_heatmap <- pheatmap(t(ML_dp), cluster_cols=FALSE, cluster_rows=FALSE)
# save_pheatmap_pdf(out_heatmap, "ML_heatmap_all_mut.pdf", 50, 10)
groups <- locations <- as.numeric(sub(":.*","",rownames(ML_dp)))
groups[groups>0] <- "00 Genomic"
groups[locations<=265] <- "01 5-UTR"
groups[locations>265 & locations<=13441] <- "02 ORF1a"
groups[locations>13468 & locations<=21555] <- "03 ORF1b" # 13442-13467 have no assignation
groups[locations>21555 & locations<=25384] <- "04 S"
groups[locations>25393 & locations<=26220] <- "05 ORF3a"
groups[locations>26245 & locations<=26472] <- "06 E"
groups[locations>26523 & locations<=27191] <- "07 M"
groups[locations>27199 & locations<=27387] <- "08 ORF6"
groups[locations>27394 & locations<=27759] <- "09 ORF7a"
groups[locations>27756 & locations<=27887] <- "10 ORF7b"
groups[locations>27894 & locations<=28259] <- "11 ORF8"
groups[locations>28274 & locations<=29533] <- "12 N"
groups[locations>28283 & locations<=28573] <- "13 ORF9b"
groups[locations>29558 & locations<=29674] <- "14 ORF10"
groups[locations>29675 & locations<=29903] <- "15 3-UTR"

custom_colours = list(groups = c("00 Genomic"="gray", "01 5-UTR"="cornflowerblue", "02 ORF1a"="coral1","03 ORF1b"="turquoise3", "04 S"="chartreuse2", "05 ORF3a"="purple", "06 E"="firebrick", "07 M"="gold2", "08 ORF6"="hotpink", "09 ORF7a"="forestgreen", "10 ORF7b"="darkblue", "11 ORF8"="darkslategray", "12 N"="darkorange2", "13 ORF9b"="turquoise3", "14 ORF10"="brown2", "15 3-UTR"="mediumpurple"))
row_labels <- data.frame(groups)
rownames(row_labels) <- rownames(ML_dp)
# Plot heatmap
out_heatmap <- pheatmap(t(log(ML_dp)), cluster_cols=FALSE, cluster_rows=FALSE, annotation_colors = custom_colours, border_col="white", annotation_col = row_labels)
save_pheatmap_pdf(out_heatmap, "ML_heatmap_vamip_tot_depth_noClust.pdf", 50, 9)
out_heatmap <- pheatmap(t(log(ML_dp)), cluster_cols=FALSE, cluster_rows=clust_col, annotation_colors = custom_colours, border_col="white", annotation_col = row_labels)
save_pheatmap_pdf(out_heatmap, "ML_heatmap_vamip_tot_depth_clust_genome.pdf", 50, 9)
out_heatmap <- pheatmap(t(log(ML_dp)), cluster_cols=clust_row, cluster_rows=clust_col, annotation_colors = custom_colours, border_col="white", annotation_col = row_labels)
save_pheatmap_pdf(out_heatmap, "ML_heatmap_vamip_tot_depth_clust_both.pdf", 50, 9)














sub <- df[grep("\\|",rownames(df)),] # First, subset those with the "|" char. These include those from delta or omicron
# sub <- sub[-grep("28247:A->-G\\|Delta:ORF8:DF119-120-",rownames(sub)),]
# sub <- sub[-grep("11282:A->-GTT\\|Omicr:ORF1a:SG3675-3676-",rownames(sub)),]
dim(sub)
# [1]  43 17927
sub <- sub[c(grep("Delta",rownames(sub)),grep("Omicr",rownames(sub))),] # so that we can subset by variant in order
sub2 <- sub
sub2[is.na(sub2)]=0
sub <- sub[,colSums(sub2)>0] # Filters those not having any
dim(sub)
# [1]    43 17747

d1 <- hclust(as.dist(1-cor(sub2, method="pearson")), method="mcquitty")
clust <- cutree(d1, h=max(d1$height/1.5))
write.table(clust,"clusters.tsv")
sub3 <- sub2[,clust>13]
sub4 <- sub[,clust>13]
clust_row <- hclust(dist(sub3), method = "mcquitty") # These versions create clusters
clust_col <- hclust(dist(t(sub3)), method = "mcquitty")
# Plot heatmap for all with only key mutations for delta/omicron
out_heatmap <- pheatmap(t(sub), cluster_cols=FALSE, cluster_rows=FALSE)
save_pheatmap_pdf(out_heatmap, "heatmap_key_mut_Delta_Omicron_all.pdf", 20, 2500)

pos_del <- grep("Delta",rownames(sub))
pos_om <- grep("Omicr",rownames(sub))
delta_omicr <- cbind("Tot_key_delta"=apply(sub,2,function(x){sum(!is.na(x[pos_del]))}), "Tot_key_omicr"=apply(sub,2,function(x){sum(!is.na(x[pos_om]))}))
delta_omicr <- cbind(delta_omicr, "Tot_key_delta-perc"= delta_omicr[,1]/14*100, "Tot_key_omicr-perc"= delta_omicr[,2]/29*100)
write.table(delta_omicr, "Key_mutations-Delta_Omicron.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)
subdeltomi <- sub[,((delta_omicr[,1]>=3)+(delta_omicr[,2]>=4))==2]
dim(subdeltomi)
# [1]  43 226
subdeltomi2 <- subdeltomi
subdeltomi2[is.na(subdeltomi2)]=0
clust_row <- hclust(dist(subdeltomi2), method = "mcquitty") # These versions create clusters
clust_col <- hclust(dist(t(subdeltomi2)), method = "mcquitty")

out_heatmap <- pheatmap(t(subdeltomi), cluster_cols=FALSE, cluster_rows=FALSE)
save_pheatmap_pdf(out_heatmap, "heatmap_key_mut_Delta_Omicron_having_both.pdf",20,75)
out_heatmap <- pheatmap(t(subdeltomi), cluster_cols=FALSE, cluster_rows=clust_col)
save_pheatmap_pdf(out_heatmap, "heatmap_key_mut_Delta_Omicron_having_both_rowclust.pdf",20,75)

###(REINICIA_aqui)








out_heatmap <- pheatmap(t(sub), cluster_cols=clust_row, cluster_rows=clust_col)
save_pheatmap_pdf(out_heatmap, "heatmap_Freq_unique_var_hclust_both_sort_lineage.pdf")
out_heatmap <- pheatmap(t(sub), cluster_cols=FALSE, cluster_rows=clust_col)
save_pheatmap_pdf(out_heatmap, "heatmap_Freq_unique_var_hclust_only_samples_sort_lineage.pdf")

# Next, we'll create a table with % of all delta or omicron mutations
sub_bin <- (sub2>0)*1
del <- grep("Delt",rownames(sub))
omi <- grep("Omic",rownames(sub))
completeness <- cbind(colSums(sub_bin[del,]), colSums(sub_bin[omi,]))
completeness <- cbind(completeness,completeness[,1]*100/length(del),completeness[,2]*100/length(omi))
colnames(completeness) <- c("Delta","Omicron","Delta%","Omicron%")
write.table(completeness, "CritMut_Del_Omi_completeness_per_genome.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)
Genomes_wCritMut_delta <- table(completeness[,1]); names(Genomes_wCritMut_delta)
Genomes_wCritMut_omicron <- table(completeness[,2])
write.table(Genomes_wCritMut_delta, "CritMut_Del_summary.tsv", sep="\t", quote=FALSE, row.names=F, col.names=T)
write.table(Genomes_wCritMut_omicron, "CritMut_Omi_summary.tsv", sep="\t", quote=FALSE, row.names=F, col.names=T)
subcomp <- completeness[(completeness[,1]>=2)+(completeness[,2]>=2)==2,]
inBoth <- sub[,rownames(subcomp)] # subset only those with at least 2 mutations from both delta and omicron
dim(inBoth)
# [1]  45 305
sub2 <- inBoth
sub2[is.na(sub2)]=0
clust_row <- hclust(dist(sub2), method = "complete") # These versions create a complete hclust
clust_col <- hclust(dist(t(sub2)), method = "complete")
save_pheatmap_pdf <- function(x, filename) { # This time, adjust the size of the output pdf file to a smaller size
  pdf(filename,12,50)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
out_heatmap <- pheatmap(t(inBoth), cluster_cols=clust_row, cluster_rows=clust_col)
save_pheatmap_pdf(out_heatmap, "heatmap_Freq_deltacron_hclust_both_sort_lineage.pdf")
out_heatmap <- pheatmap(t(inBoth), cluster_cols=FALSE, cluster_rows=clust_col)
save_pheatmap_pdf(out_heatmap, "heatmap_Freq_deltacron_hclust_only_samples_sort_lineage.pdf")

meta <- read.table("Metadata_all.tsv", header=T, sep ='\t', stringsAsFactors = FALSE, row.names=1, check.names=F, quote="")
names <- sub("\\|.*","",colnames(inBoth))
dates <- meta[names,"Fecha de recoleccion"]
colnames(inBoth) <- paste0(dates,"|",colnames(inBoth))
inBoth <- inBoth[,order(colnames(inBoth))]
out_heatmap <- pheatmap(t(inBoth), cluster_cols=FALSE, cluster_rows=FALSE)
save_pheatmap_pdf(out_heatmap, "heatmap_Freq_deltacron_hclust_bydate.pdf")

# ### Version 2: Total Depths ###
library("pheatmap")
df <- read.table("Alt_variant_TotalDepth.tsv", sep="\t",header=T, skip=0, comment.char='',quote="",fill=F,row.names=1, check.names=FALSE) # Load input table
dim(df)
# [1] 262 144
df2 <- df # Create a copy of the input matrix
df2[is.na(df2)]=0 # to remove NAs so we can cluster items (hclust accepts no non.numeric or NA-including table)
clust_row <- hclust(dist(df2), method = "complete") # These versions create a complete hclust for rows (alleles/mutations; rows)
clust_col <- hclust(dist(t(df2)), method = "complete") # likeswise, this hclust is for the genomes (columns)
# Other options include "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
font = 7# adjust font
# for this, we also need log10 values instead
out_heatmap <- pheatmap(t(log10(df)), cluster_cols=clust_row, cluster_rows=clust_col) # since it plots in the opposite direction, we use a transposed matrix and swap rows columns
save_pheatmap_pdf <- function(x, filename) { # this simple function is used for plotting the heatmap
  pdf(filename,30,30)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(out_heatmap, "heatmap_TotalDepth_hclust_both.pdf") # This clusters both the mutations and genomes

out_heatmap <- pheatmap(t(log10(df)), cluster_cols=FALSE, cluster_rows=clust_col) # Now, repeat but this time, only cluster by genomes (leave mutations in the original order [sorted by position])
save_pheatmap_pdf(out_heatmap, "heatmap_TotalDepth_hclust_only_samples.pdf") # This clusters only genomes, not mutations

# Again, for the most important mutations only
save_pheatmap_pdf <- function(x, filename) { # This time, adjust the size of the output pdf file to a smaller size
  pdf(filename,10,25)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
sub <- df[grep("\\|",rownames(df)),] # First, subset those with the "|" char. These include those from delta or omicron
dim(sub)
# [1]  43 144
sub <- sub[c(grep("Delta",rownames(sub)),grep("Omicr",rownames(sub))),] # so that we can subset by variant in order
sub2 <- sub
sub2[is.na(sub2)]=0
clust_row <- hclust(dist(sub2), method = "complete") # These versions create a complete hclust
clust_col <- hclust(dist(t(sub2)), method = "complete")
out_heatmap <- pheatmap(t(log10(sub)), cluster_cols=clust_row, cluster_rows=clust_col)
save_pheatmap_pdf(out_heatmap, "heatmap_TotalDepth_unique_var_hclust_both_sort_lineage.pdf")
out_heatmap <- pheatmap(t(log10(sub)), cluster_cols=FALSE, cluster_rows=clust_col)
save_pheatmap_pdf(out_heatmap, "heatmap_TotalDepth_unique_var_hclust_only_samples_sort_lineage.pdf")
