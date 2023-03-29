# Started 2022-04-29 by Rodrigo Garcia-Lopez for GAL, iBT, UNAM
# ### FUNCTIONS ###
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
save_pheatmap_pdf <- function(x, filename, width, height) { # This time, adjust the size of the output pdf file to a smaller size
  pdf(filename,width,height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
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
sum(!is.na(df)) # we test the total number of items
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




# COMPARE PREVALENCE
# First, we want to compare a raw view of all totals, with and without the vamips:
compare_all <- cbind("All"=rowSums(group_prevalence(group_map(all_sex),df)$raw), "only_VAMIPs"=rowSums(group_prevalence(group_map(all_sex),df_vamip)$raw),"only_Consensus"=rowSums(group_prevalence(group_map(all_sex),df_novamip)$raw))
compare_all_total <- cbind(compare_all, compare_all/ncol(df)*100)
colnames(compare_all_total) <- c(colnames(compare_all_total)[1:3],paste(colnames(compare_all_total)[1:3],"%"))
write.table(compare_all_total, "08_Analyses/VAMIP_vs_Consensus.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)

pdf("08_Analyses/Histogram_mutations.pdf")
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
g_month <- group_prevalence(group_map(all_Months),df)
g_month_vamip <- group_prevalence(group_map(all_Months),df_vamip)
g_month_novamip <- group_prevalence(group_map(all_Months),df_novamip)
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




























# Now, for variants (general)
g_var <- group_prevalence(group_map(all_mainvar),df)
g_var_vamip <- group_prevalence(group_map(all_mainvar),df_vamip)
g_var_novamip <- group_prevalence(group_map(all_mainvar),df_novamip)
sub <- head(g_var$rel[names(sort(rowSums(g_var$rel),decreasing=TRUE)),],2000) # get the 2000 most common mutations only
sub[sub==0] <- NA
sub_vamip <- head(g_var_vamip$rel[names(sort(rowSums(g_var_vamip$rel),decreasing=TRUE)),],2000)
sub_vamip[sub_vamip==0] <- NA
sub_novamip <- head(g_var_novamip$rel[names(sort(rowSums(g_var_novamip$rel),decreasing=TRUE)),],2000)
# with these, we can get an overhead view of all the most common mutations
sub_novamip[sub_novamip==0] <- NA
out_heatmap <- pheatmap(sub, cluster_cols=FALSE, cluster_rows=FALSE)
save_pheatmap_pdf(out_heatmap, "08_Analyses/Variant_heatmap_2000TopMut.pdf", 9, 230)
out_heatmap <- pheatmap(sub_vamip, cluster_cols=FALSE, cluster_rows=FALSE)
save_pheatmap_pdf(out_heatmap, "08_Analyses/Variant_vamip_heatmap_2000TopMut.pdf", 9, 230)
out_heatmap <- pheatmap(sub_novamip, cluster_cols=FALSE, cluster_rows=FALSE)
save_pheatmap_pdf(out_heatmap, "08_Analyses/Variant_novamip_heatmap_2000TopMut.pdf", 9, 230)
# Now, use the 2000 most abundant items in sub_novamip to extract them from the novamip set
sub_novamip_same_as_sub_vamip <- g_var_novamip$rel[rownames(sub_vamip),]
sub_novamip_same_as_sub_vamip[sub_novamip_same_as_sub_vamip==0] <- NA
out_heatmap <- pheatmap(sub_novamip_same_as_sub_vamip, cluster_cols=FALSE, cluster_rows=FALSE)
save_pheatmap_pdf(out_heatmap, "08_Analyses/Variant_novamip_same_set_as_vamip_heatmap_2000TopMut.pdf", 9, 230)
# Compare those with vamips and novamips:
rownames(sub_novamip_same_as_sub_vamip) <- paste(rownames(sub_novamip_same_as_sub_vamip),"X novamip") # First, rename both sets
rownames(sub_novamip_same_as_sub_vamip) <- paste(sprintf("%04d", 1:2000),rownames(sub_novamip_same_as_sub_vamip))
rownames(sub_vamip) <- paste(rownames(sub_vamip),"VAMIP")
rownames(sub_vamip) <- paste(sprintf("%04d", 1:2000),rownames(sub_vamip))
compare <- rbind(sub_vamip,sub_novamip_same_as_sub_vamip)
compare <- compare[order(rownames(compare)),]
out_heatmap <- pheatmap(compare, cluster_cols=FALSE, cluster_rows=FALSE)
save_pheatmap_pdf(out_heatmap, "08_Analyses/Variant_heatmap_compare_2000TopMut.pdf", 9, 460)

# We now want to create a map for tracing each mutations throughout months
longitudinal_map_months <- group_map(all_Months)
test_mutation <- grep("23604:C->A", rownames(df)) # mutation 23604:C->A|Omicr:S:P681H has 7628 observations in the table, not only from omicron but throughout months
test_mut_freqs <- df[test_mutation,]
lapply(longitudinal_map_months,function(x){temp <- test_mut_freqs[x]; temp <- temp[!is.na(temp)]; if(length(temp)==0){temp <- 0}; quantile(temp,seq(0,1,0.1))})
test <- lapply(longitudinal_map_months,function(x){temp <- test_mut_freqs[x]; temp <- temp[!is.na(temp)]; if(length(temp)==0){temp <- 0}; quantile(temp,seq(0,1,0.1))})
test <- lapply(longitudinal_map_months,function(x){temp <- test_mut_freqs[x]; temp <- temp[!is.na(temp)]; if(length(temp)==0){temp <- 0}; temp})
test <- unlist(lapply(longitudinal_map_months,function(x){temp <- test_mut_freqs[x]; temp <- temp[!is.na(temp)]; length(temp)}))
test

plot_single_mutation <- function(vect, map, item_name="Item", comparison_name="Comparison", min_nonzero=2){ # Input items are a frequency vector with NAs (vect), a list containing positions per group (e.g. months, age, etc.; map), names for the actual vector (mutation are expected; item_name) and the category that is explored (comparison_name), and minimum expected subcategories (groups) that should have the item (min_nonzero).
	list_totals <- unlist(lapply(map,function(x){temp <- vect[x]; temp <- temp[!is.na(temp)]; length(temp)}))
	if(sum(list_totals > 0) >= min_nonzero){
		if(max(list_totals)<10){return()} # only consider mutations where the max genomes per category (e.g. any given month) equal at least this
		if(sum(list_totals)<30){return()} # only consider mutations that have at least 30 total items
		list_out <- lapply(map,function(x){temp <- vect[x]; temp <- temp[!is.na(temp)]; if(length(temp)==0){temp <- 0}; temp}) # for each item in the map, extract non-NA values from the vector into a list output
		dir.create(comparison_name,showWarnings = FALSE)
		ref_totals <- unlist(lapply(map, length)) # get the totals per subcategory (e.g. month)
		item_name_new <- sub("\\|", "_", item_name) # Fix characters for file names
		item_name_new <- gsub(":", "_", item_name_new)
		item_name_new <- gsub(" ", "", item_name_new)
		item_name_new <- gsub(",", "n", item_name_new)
		item_name_new <- gsub(">", "", item_name_new)
		item_name_new <- gsub("-", "t", item_name_new)
		pdf(paste0(comparison_name, "/", comparison_name,"_",item_name_new,".pdf"), width=14)
		par(oma = c(1, 1, 1, 10)) # This is just a creative fix to plot the legend outside
# 		par(mar=c(5,4,2,5))
		bp <- barplot(list_totals/ref_totals, las=2, main=paste(item_name, "by", comparison_name), ylab="Frequency", ylim=c(-0.1,1.1), col="coral1", border=NA, yaxt='n')
		boxplot(list_out, add=T, las=2, at=bp, border="navyblue", col=NA, outline = T, yaxt='n')
		axis(2,at=seq(0,1,0.1), las=2)
		axis(4,at=seq(0,1,0.1), labels=seq(0,1,0.1)*max(list_totals), las=2, col="chartreuse4",col.axis="chartreuse4")
		points(bp,list_totals/max(list_totals), col="chartreuse3", pch=17, cex=1.5)
		par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0.5), mar = c(0, 0, 0, 0), new = TRUE)
		plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
		legend("right", legend = c("Mut Freq","Allele Freq","Total with Mut"), xpd = TRUE, horiz = FALSE, inset = c(0,0), bty = "n", col = c("coral1", "navyblue", "chartreuse3"), lty = c(NA,1,NA), cex = 0.8, pch=c(15,NA,17))
		dev.off()
		return(paste0("Created: ",comparison_name,"_",item_name_new,".pdf"))
	}
}

save.image("checkpoint1.Rdata")
load("checkpoint1.Rdata")

for(mut in rownames(df)){
	print(mut)
	plot_single_mutation(df[mut,], longitudinal_map_months, mut, "Month_14_fullname",14)
}

for(mut in rownames(df)){
	print(mut)
	plot_single_mutation(df[mut,], group_map(all_Weeks), mut, "Weeks_4",4)
}

list_ages <- group_map(all_age)
new_names <- sprintf("%03d", (as.numeric(names(list_ages))))
new_names[length(new_names)] <- "u"
# new_names <- paste0("a",new_names)
names(list_ages) <- new_names
list_ages <- list_ages[order(names(list_ages))]
for(mut in rownames(df)){
	print(mut)
	plot_single_mutation(df[mut,], list_ages, mut, "Age",2)
}

for(mut in rownames(df)){
	print(mut)
	plot_single_mutation(df[mut,], group_map(all_age_vac), mut, "Age_vac",1)
}

list_var <- group_map(all_mainvar)
names(list_var)
# [1] "519"     "Alpha"   "Delta"   "Gamma"   "Omicron" "Others"  "Recomb"
list_var <- list_var[c(1,2,4,3,5,7,6)]
for(mut in rownames(df)){
	print(mut)
	plot_single_mutation(df[mut,], list_var, mut, "Variants",1)
}

for(mut in rownames(df)){
	print(mut)
	plot_single_mutation(df[mut,], group_map(all_sex), mut, "Sex",1)
}

for(mut in rownames(df)){
	print(mut)
	plot_single_mutation(df[mut,], group_map(all_region7), mut, "Region_7",1)
}

for(mut in rownames(df)){
	print(mut)
	plot_single_mutation(df[mut,], group_map(all_status), mut, "Status",1)
}

deltas <- group_map(all_mainvar)$Delta
deltas <- deltas[deltas > 100]
deltas <- df[,deltas]
temp <- deltas
temp[is.na(temp)] <- 0
deltas <- deltas[rowSums(temp)>0,]

omicrons <- group_map(all_mainvar)$Omicron
omicrons <- df[,omicrons]
temp <- omicrons
temp[is.na(temp)] <- 0
omicrons <- omicrons[rowSums(temp)>0,]
rm("temp")
rm("df")
save.image("checkpoint2.Rdata")
