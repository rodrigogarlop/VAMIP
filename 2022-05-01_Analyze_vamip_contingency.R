# Started 2022-04-29 by Rodrigo Garcia-Lopez for GAL, iBT, UNAM
# ### FUNCTIONS ###
find_items <- function(queryStrings, dataf, colname){ # IMPORTANT. This is used to avoid partial index matching. Use it to hash a table and use rownames to match a vector. colname may be numeric
	searchVect <- dataf[,colname]
	names(searchVect) <- rownames(dataf)
	out <- searchVect[unlist(queryStrings)]
	return(out)
}
group_map <- function(vect) { # Get a mutation table (columns = genomes; rows = mutations), and a vector of the same rowsize. Then use the vector to define how groups should be splitted and build a list of df with the table observations per group (a map of positions)
	group_names <- levels(as.factor(vect))
	group_positions <- sapply(group_names,function(x){grep(x,vect)})
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
		out_rel[,i] <- out_raw[,i]/nrow(single_grp)
	}
	out <- list(out_raw,out_rel)
	names(out) <- c("raw","rel")
	return(out)
}

# ### LOAD LIBRARIES & INPUT ###
library("pheatmap")
library("stringr")
df <- read.table("Alt_variant_Freq.tsv", sep="\t",header=T, skip=0, comment.char='',quote="",fill=F,row.names=1, check.names=FALSE) # Load input table
# Each column was named as follows: "2021-03-04|L004_202101083013|B.1.1.222(20B_Others)|046|M|Puebla" and they have been sorted before
df<- df[,order(colnames(df))] # In case they were not sorted, we do it here
sum(!is.na(df))
# [1] 1034529
df[df==1e-06] <- NA # First off, remove items that are below the cutoff (this was set when the input table was created to make them 1e-06)
sum(!is.na(df))
# [1] 872199
dim(df)
# [1]	7340 17927
meta <- read.table("07_metadata/Metadata_all.tsv", header=T, sep ='\t',stringsAsFactors = FALSE, check.names=F, fill=T, quote="") # Load metadata with id folio for xref. IDs should be unique
rownames(meta) <- meta[,"ID Folio"]
dim(meta)
# [1] 19024		48
# This may be a larger table but only matching items will be kept

# ### MAIN ###
# Create two copies to compare only vamips (<0.05) and only major variants
df_vamip <- df; df_vamip[df_vamip>0.5] <- NA # this will hold vamips (>0.5 are now NAs)
df_novamip <- df; df_novamip[df_novamip<=0.5] <- NA # this will hold vamips (>0.5 are now NAs)

# PREPARE METADATA
# Process the name to extract basic info and the Folio ID
# names are processed into a table: Date|Folio ID|Pango(clade)|age|gender|State
# Example: "2021-03-04|L004_202101083013|B.1.1.222(20B_Others)|046|M|Puebla"
minimetatab <- t(as.data.frame(str_split(colnames(df), "\\|")))
rownames(minimetatab) <- NULL # row names are removed
# Now, explore basic dataset composition per genome
# Extract dates from the column names (each genome name includes the sample collection date)
# The rest of the metadata needs to be extracted from the external metadata table using the ID Folio for xref
meta <- meta[unique(minimetatab[,2]),]
dim(meta)
# [1] 17927		48
write.table(meta, "07_metadata/Metadata_matched_only.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)

# Some metadata are already found in the genome's name and will be processed directly, others must be extracted from the external meta table
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

# CLADES & VARIANTS
all_mainvar <- find_items(minimetatab[,2],meta, "Variants") # Start with the main variants (VOCs and 519)
all_clade <- find_items(minimetatab[,2],meta, "Clado Nexstrain") # Next, the main variants

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

# COMPARE PREVALENCE
# First, we want to compare a raw view of all totals, with and without the vamips:
compare_all <- cbind("All"=rowSums(group_prevalence(group_map(all_sex),df)$raw), "only_VAMIPs"=rowSums(group_prevalence(group_map(all_sex),df_vamip)$raw),"only_Consensus"=rowSums(group_prevalence(group_map(all_sex),df_novamip)$raw))
# First, by month
g_month <- group_prevalence(group_map(all_Months),df)
# with this, we can get an overhead view of all the most common mutations

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
pdf("08_Analyses/Total_genomes-age.pdf",width=14)
	barplot(las=2,table(all_age),border=NA,col="forestgreen",yaxt='n',ylim=c(0,400))
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

