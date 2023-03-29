#Started 2022-09-27 by Rodrigo Garcia-Lopez for GAL, iBT, UNAM: We will get percentage of VAMIPS vs non vamips
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

save.image("Genomes-VAMIPsVsSNPs_chkpt1.Rdata")

# test <- apply(df,2,function(x){quantile(x[!is.na(x)],seq(0,1,0.01))})
SNPs <- apply(df,2,function(x){sum(x[!is.na(x)]>=0.5)})
VAMIPs <- apply(df,2,function(x){sum(x[!is.na(x)]<0.5)})
out <- cbind(SNPs, VAMIPs); out <- cbind(out, "Total"=rowSums(out), prop.table(out,1)); colnames(out)[4:5] <-  paste(colnames(out)[4:5], "%") # Create an output table
out <- out[!is.na(out[,5]),] # remove some NAs
# hist(out[,3],breaks=100,xlim=c(0,100))
write.table(out, "2022-09-19_SNPs_vs_VAMIPs_perGenome.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)

hist(out[,3])
abline(v=mean(out[,3])+2*sd(out[,3]))
abline(v=mean(out[,3])+3*sd(out[,3]))
abline(v=mean(out[,3])-2*sd(out[,3]))
abline(v=mean(out[,3])-3*sd(out[,3]))

# again, but repeat with vamips <0.1
df[df<0.1] <- NA
sum(!is.na(df)) # we test the total number of items
# [1] 812811
dim(df)
#  7340 17927

SNPs <- apply(df,2,function(x){sum(x[!is.na(x)]>=0.5)})
VAMIPs <- apply(df,2,function(x){sum(x[!is.na(x)]<0.5)})
out <- out[!is.na(out[,5]),] # remove some NAs
out <- cbind(SNPs, VAMIPs); out <- cbind(out, "Total"=rowSums(out), prop.table(out,1)); colnames(out)[4:5] <-  paste(colnames(out)[4:5], "%") # Create an output table
# hist(out[,3],breaks=100,xlim=c(0,100))
write.table(out, "2022-09-19_SNPs_vs_VAMIPs_perGenome_Freq0.1.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)

# UPDATE 2022-09-28: We now want to get which are the main mutations per genome
# FUNCTIONS
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
# First, load the raw frequency table (% of each mutation seen per position in the genome):
load("Genomes-VAMIPsVsSNPs_chkpt1.Rdata")
sum(!is.na(df)) # we test the total number of items
# [1] 872199
dim(df)
# [1]  7340 17927
# Now, load the metadata to determine lineages
# PREPARE METADATA
meta <- read.table("07_metadata/Metadata_all.tsv", header=T, sep ='\t',stringsAsFactors = FALSE, check.names=F, fill=T, quote="") # Load metadata with id folio for xref. IDs should be unique
rownames(meta) <- meta[,"ID Folio"]
dim(meta)
# [1] 19024		48
# Process the name to extract basic info and the Folio ID
# names are processed into a table: Date|Folio ID|Pango(clade)|age|gender|State
# Example: "2021-03-04|L004_202101083013|B.1.1.222(20B_Others)|046|M|Puebla"
minimetatab <- t(as.data.frame(str_split(colnames(df), "\\|")))
rownames(minimetatab) <- NULL # row names are removed
# This may be a larger table but only matching items will be kept
meta <- meta[unique(minimetatab[,2]),]
dim(meta)
# [1] 17927    48

all_lineages <- find_items(minimetatab[,2],meta, "New_Lineage") # Load the lineage metadata (this is drawn from the "meta" object since the one in the colname may be oudated
lineage_map <- group_map(all_lineages) # Create a map of which columns belong to which lineages
lineage_totals <- sort(unlist(lapply(lineage_map, length)), decreasing=TRUE) # Now, count them
n=20 # Define the least number of genomes for each variant to be considered for the analysis
pass_yes <- names(which(lineage_totals>=n)) # Get which items have at least n genomes
pass_not <- names(which(lineage_totals<n)) # And which ones don't
lineage_map_filt <- lineage_map[pass_yes] # Start by creating a list with those items that pass the filter
temp <- unlist(lineage_map[pass_not[grep("B\\.1\\.631",pass_not)]]);names(temp) <- NULL; lineage_map_filt[["B.1.631"]] <- temp # and append some additional items for general BA.1 subvariants
temp <- unlist(lineage_map[pass_not[grep("BA\\.1",pass_not)]]);names(temp) <- NULL; lineage_map_filt[["BA.1.others"]] <- temp # and append some additional items for general BA.1 subvariants
temp <- unlist(lineage_map[pass_not[grep("BA\\.2",pass_not)]]);names(temp) <- NULL; lineage_map_filt[["BA.2.others"]] <- temp # as well as BA.2 subvariants
temp <- unlist(lineage_map[pass_not[c(grep("B\\.1\\.617\\.1",pass_not),grep("AY",pass_not))]]);names(temp) <- NULL; lineage_map_filt[["AY.others"]] <- temp # and other deltas
compare_variant_fil <- group_prevalence(lineage_map_filt,df) # Now, get the prevalence of each mutation per variant
write.table(compare_variant_fil[[1]], "2022-09-28_Main_mutations_per_variant-filtered.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)
compare_variant_fil_rel <- round(compare_variant_fil[[2]],2)/100
write.table(compare_variant_fil_rel, "2022-09-28_Main_mutations_per_variant-filtered-rel.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)
# multiMatrix <- lapply(rev(seq(50,95,5)), function(x) compare_variant_fil_rel>=x) # Create multiple tables, each with a different cuttoff

multiMatrix <- list() # create a new list, each item will be a matrix with varying numbers of representative mutations per variant (filtered at different permissive levels [percentages] of the genomes that have them per variant)
for (i in rev(seq(.50,.95,.05))) { # We will just be using those in most (>50% samples)
	print(i)
	multiMatrix[[paste0("m",i)]] <- compare_variant_fil_rel>=i # each output in the list will be a logical matrix with various levels of permissiveness
};rm(i)

# MORE FUNCTIONS
variant_completeness <- function(freqs_vect,pos_vect) { # extract only the items in pos_vect positions of vector freqs_vect, then calculate which have not NAs in those positions, the output is a single %
	test_vect <- freqs_vect[pos_vect]
	MutObs <- sum(!is.na(test_vect))
	MutExp <- length(test_vect)
	ObsPerc <- round(MutObs/MutExp*100,2)
	out <- list(MutObs, MutExp, ObsPerc)
	names(out) <- c("MutObs", "MutExp", "ObsPerc")
	return(out)
}

for(mat in names(multiMatrix)){
	out_mat <- t(apply(df, 2, function(genome) {apply(multiMatrix[[mat]], 2, function(var){variant_completeness(genome,var)$ObsPerc})})) # test all variants on each genome and create a table bearing the completeness
	write.table(out_mat, paste0("genome_variant_completeness_", mat, "_confidence.tsv"), sep="\t", quote=FALSE, row.names=T, col.names=NA)
}

# UPDATE 2022-09-29: We will now do the variant analysis (main mutations) in a more general way, using for example Delta, Omicron, etc.
meta[meta[,"Variants"]=="Recomb","Variants"] <- "XB" # First, since all recombinants are XBs, we use that lineage name instead
all_VOCs <- find_items(minimetatab[,2],meta, "Variants") # Load the VOC metadata (this is drawn from the "meta" object since the one in the colname may be oudated
VOC_map <- group_map(all_VOCs) # Create a map of which columns belong to which VOCs
VOC_totals <- sort(unlist(lapply(VOC_map, length)), decreasing=TRUE) # Now, count them
n=20 # Define the least number of genomes for each variant to be considered for the analysis
pass_yes <- names(which(VOC_totals>=n)) # Get which items have at least n genomes
pass_not <- names(which(VOC_totals<n)) # And which ones don't
VOC_map_filt <- VOC_map[pass_yes] # Start by creating a list with those items that pass the filter
compare_variant_fil <- group_prevalence(VOC_map_filt,df) # Now, get the prevalence of each mutation per
write.table(compare_variant_fil[[1]], "2022-09-28_Main_mutations_per_VOC-filtered.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)
compare_variant_fil_rel <- round(compare_variant_fil[[2]],2)/100
write.table(compare_variant_fil_rel, "2022-09-28_Main_mutations_per_VOC-filtered-rel.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)

multiMatrix <- list() # create a new list, each item will be a matrix with varying numbers of representative mutations per variant (filtered at different permissive levels [percentages] of the genomes that have them per variant)
for (i in rev(seq(.50,.95,.05))) { # We will just be using those in most (>50% samples)
	print(i)
	multiMatrix[[paste0("m",i)]] <- compare_variant_fil_rel>=i # each output in the list will be a logical matrix with various levels of permissiveness
};rm(i)

for(mat in names(multiMatrix)){
	out_mat <- t(apply(df, 2, function(genome) {apply(multiMatrix[[mat]], 2, function(var){variant_completeness(genome,var)$ObsPerc})})) # test all variants on each genome and create a table bearing the completeness
	write.table(out_mat, paste0("genome_VOC_completeness_", mat, "_confidence.tsv"), sep="\t", quote=FALSE, row.names=T, col.names=NA)
}
