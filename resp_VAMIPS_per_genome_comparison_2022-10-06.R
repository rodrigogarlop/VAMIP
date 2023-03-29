#Started 2022-10-06 by Rodrigo Garcia-Lopez for GAL, iBT, UNAM: We will get percentage of VAMIPS vs non vamips
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
# There is an item that ends with 0 mutations, we will remove it)
df <- df[apply(df, 1, function(x){sum(x[!is.na(x)])})!=0, apply(df, 2, function(x){sum(x[!is.na(x)])})!=0]
dim(df)
# [1]	7340 17926
save.image("Genomes-VAMIPsVsGenomes_chkpt1.Rdata")

# test <- apply(df,2,function(x){quantile(x[!is.na(x)],seq(0,1,0.01))})
SNPs <- apply(df,2,function(x){sum(x[!is.na(x)]>=0.5)})
VAMIPs <- apply(df,2,function(x){sum(x[!is.na(x)]<0.5)})
out <- cbind(SNPs, VAMIPs); out <- cbind(out, "Total"=rowSums(out), prop.table(out,1)); colnames(out)[4:5] <-  paste(colnames(out)[4:5], "%") # Create an output table
out <- out[!is.na(out[,5]),] # remove some NAs
# hist(out[,3],breaks=100,xlim=c(0,100))
write.table(out, "2022-09-19_SNPs_vs_VAMIPs_perGenome.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)

hist(out[,3],breaks=50,xlim=c(0,100))
abline(v=mean(out[,3])+2*sd(out[,3]))
abline(v=mean(out[,3])+3*sd(out[,3]))
abline(v=mean(out[,3])-2*sd(out[,3]))
abline(v=mean(out[,3])-3*sd(out[,3]))

sum(!is.na(df)) # we test the total number of items
# [1] 872199
# again, but repeat with vamips <0.1
df[df<0.1] <- NA
sum(!is.na(df)) # we test the total number of items
# [1] 812811
dim(df)
#  7340 17926

SNPs <- apply(df,2,function(x){sum(x[!is.na(x)]>=0.5)})
VAMIPs <- apply(df,2,function(x){sum(x[!is.na(x)]<0.5)})
out <- cbind(SNPs, VAMIPs); out <- cbind(out, "Total"=rowSums(out), prop.table(out,1)); colnames(out)[4:5] <-  paste(colnames(out)[4:5], "%") # Create an output table
out <- out[!is.na(out[,5]),] # remove some NAs
# hist(out[,3],breaks=100,xlim=c(0,100))
write.table(out, "2022-09-19_SNPs_vs_VAMIPs_perGenome_Freq0.1.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)

hist(out[,3],breaks=50,xlim=c(0,100))
abline(v=mean(out[,3])+2*sd(out[,3]))
abline(v=mean(out[,3])+3*sd(out[,3]))
abline(v=mean(out[,3])-2*sd(out[,3]))
abline(v=mean(out[,3])-3*sd(out[,3]))

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
load("Genomes-VAMIPsVsGenomes_chkpt1.Rdata")
sum(!is.na(df)) # we test the total number of items
# [1] 872199
dim(df)
# [1]  7340 17926
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
# [1] 17926    48

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
meta[meta[,"New_Lineage"]=="B.1.1.222","Variants"] <- "222" # Also, add variant 222
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

#UPDATE 2022-10-04: We now want to evaluate genome by genome to determine the % of VAMIPs
load("Genomes-VAMIPsVsGenomes_chkpt1.Rdata")
dir.create("genomePos_vs_SNPnVAMIP_contents", showWarnings=FALSE)
df[is.na(df)] <- 0 # We need to remove NAs before anything else
VAMIPs <- SNPs <- df
SNPs[SNPs<0.5] <- 0 # First, get two matrices, one with only SNPs
VAMIPs[VAMIPs>=0.5] <- 0 # and the other one with only VAMIPs
# Now, we need to collate each table where mutations are summed by position. We only want to keep substitutions and indels separated
all_mut <- as.data.frame(cbind("id"=1:nrow(df),"Pos"=sub(":.*","",rownames(df)),"MutType"=rep("_Sub", nrow(df))));rownames(all_mut) <- rownames(df) # we start with the full collection of mutation names, and the type of mutations
# now, fill the indels accordingly to identify the type of mutation
all_mut[grep(">\\+",rownames(df)),3] <- "-Ins"
all_mut[grep(">-",rownames(df)),3] <- "-Del"
SNPs <- rowsum(SNPs, group=paste0(all_mut[,2],all_mut[,3])) # Collate by the postition and type of mutation
dim(SNPs)
# [1]  6762 17926
VAMIPs <- rowsum(VAMIPs, group=paste0(all_mut[,2],all_mut[,3])) # Collate by the postition and type of mutation
dim(VAMIPs)
# [1]  6762 17926
SNPs[SNPs>1] <- 1; VAMIPs[VAMIPs>1] <- 1 # Due to rounding, some items may have values slightly larger than 1 (<1e-5), thus, we will truncate those to 1s

id_Sub <- grep("Sub",rownames(SNPs))
id_Del <- grep("Del",rownames(SNPs))
id_Ins <- grep("Ins",rownames(SNPs))
length(c(id_Sub,id_Del,id_Ins)) # This checkpoint should be the same as the number of rows
# [1] 6762
fill_missing_items <- function(mat,vector){ # Gets a table and a vector with all rows that should be present, outputs an expanded row collection with missing dates for complete calendar in that range. rownames should have date format as %Y-%m-%d
	xtable <- as.data.frame(matrix(0,nrow=length(vector),ncol=ncol(mat)), stringsAsFactors = FALSE) # Create empty vessel for output
	rownames(xtable) <- vector # use the vector as rownames
	colnames(xtable) <- colnames(mat) # and inherit the names
	invisible(sapply(rownames(mat),function(x) {xtable[x,] <<- mat[x,]}))	# append the original values in the corresponding places (write to higher env variable
	return(xtable)
}
strip <- function(string){ # Takes a string object and strips rare characters
	string <- iconv(string,from="UTF-8",to="ASCII//TRANSLIT") # Remove accents and rare chars
# 	string <- tolower(string) # change to lower
	string <- gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", string, perl=TRUE) # Remove multiple spaces
	string <- gsub("\\+","",string) # remove + characters
	string <- gsub("\\|", "_-_",string)
	string <- gsub(" ", "_",string)
	string <- gsub("\\(", "_",string)
	string <- gsub("\\)", "_",string)
	string <- gsub("\\/", "_",string)
	string <- gsub("__", "_",string)
	return(string)
}

rm(df)
save.image("Genomes-VAMIPsVsGenomes_chkpt2.Rdata")
# i <- 4527 # test
load("Genomes-VAMIPsVsGenomes_chkpt2.Rdata")
# for(i in 1:ncol(SNPs)){ # For each genome
for(i in 15001:ncol(SNPs)){ # For each genome
	print(i)
	name <- colnames(SNPs)[i] # This holds the name of the genome
	all_mut <- cbind("SNPs"=SNPs[,i], "VAMIPs"=VAMIPs[,i]); rownames(all_mut) <- rownames(SNPs) # this holds the SNPs and VAMIPs for the genome
	sub <- all_mut[id_Sub,]; rownames(sub) <- sub("_Sub","",rownames(sub)) # The table is split in 3 subtables
	del <- all_mut[id_Del,]; rownames(del) <- sub("-Del","",rownames(del))
	ins <- all_mut[id_Ins,]; rownames(ins) <- sub("-Ins","",rownames(ins))
	all_pos <- sprintf("%05d", 1:29903) # Add a vector for all available positions

	sub <- fill_missing_items(sub,all_pos); colnames(sub) <- paste("Sub",colnames(sub)) # create a full table containing 0s for missing items, then rename the columns
	del <- fill_missing_items(del,all_pos); colnames(del) <- paste("Del",colnames(del))
	ins <- fill_missing_items(ins,all_pos); colnames(ins) <- paste("Ins",colnames(ins))
	data <- cbind(sub, del, ins)
# write.table(data[rowSums(data)>0,], paste0("genomePos_vs_SNPnVAMIP_contents/", str_split(name, "\\|")[[1]][2], ".tsv"), sep="\t", quote=FALSE, row.names=T, col.names=NA)
	write.table(data[rowSums(data)>0,], paste0("genomePos_vs_SNPnVAMIP_contents/", gsub("\\|", "_-_",name), ".tsv"), sep="\t", quote=FALSE, row.names=T, col.names=NA)
	pdf(paste0("genomePos_vs_SNPnVAMIP_contents/", gsub("\\|", "_-_",name), ".pdf"),width=18, height=10)
	cols <- c("turquoise","blue","pink","violetred","chartreuse3","forestgreen")
	ltys <- c(1,1,2,2,4,4)
	par(oma = c(1.5, 0, 1, 6)) # This is just a creative fix to plot the legend outside
	par(mar=c(5.1, 5.1, 4.1, 2.1))
	matplot(las=1, data[,ncol(data):1],type='l', col=rev(cols),lty=rev(ltys),lwd=, ylim=c(0,1), xlim=c(1,29903), ylab="Mutation frequency per position", main=name, xaxt='n', yaxt='n',frame.plot=FALSE, xlab="Genomic Position (Nt)", cex.lab=1.5, cex.main=1.5)
	axis(1, at=seq(1,31000,2500), labels=seq(0,31000,2500),cex=1.5)
	axis(2, las=1, at=seq(0,1,.1),cex.axis=1.2)
	Genome <- diff(c(1, 266, 13468, 21556, 21563, 25385, 25393, 26221, 26245, 26473, 26523, 27192, 27202, 27388, 27394, 27756, 27888, 27894, 28260, 28274, 29534, 29558, 29675, 29903)) # This vector holds each interval (all ORFs). Endings were adjusted with +1 for calculations
	cols_bar = c("01 5-UTR"="cornflowerblue", "02 ORF1a"="coral1","03 ORF1b"="turquoise3", "Sep1"=NA, "04 S"="chartreuse2", "Sep2"=NA, "05 ORF3a"="purple", "Sep3"=NA, "06 E"="firebrick", "Sep4"=NA, "07 M"="gold2", "Sep5"=NA, "08 ORF6"="hotpink", "Sep6"=NA, "09 ORF7a"="forestgreen", "10 ORF7b"="darkblue", "Sep7"=NA, "11 ORF8"="darkslategray", "Sep8"=NA, "12 N"="darkorange2", "Sep9"=NA, "14 ORF10"="brown2", "15 3-UTR"="mediumpurple")
	sub_cols <- cols_bar[!is.na(cols_bar)]
	par(oma = c(6, 0, 43, 6)) # Adjust for plotting the genes
	barplot(las=1,cbind(Genome), horiz = TRUE, beside = FALSE, col=cols_bar, border=NA, xaxt='n',add=T)
	par(fig = c(0, 1, 0, 1), oma = c(1.5, 0, 1, 1.5), mar = c(8.1, 0, 7.1, 0), new = TRUE)
	plot(0,0,type = "n", bty = "n", xaxt = "n", yaxt = "n",xlab="", ylab="")
	legend("topright",legend=c("Sub Fix","Sub VAMIP","Del Fix","Del VAMIP","Ins Fix","Ins VAMIP"), lty=ltys, col=cols, title="Lines", lwd=3, bg="white")
	legend("right",legend=names(sub_cols), pch=15, col=sub_cols, pt.cex=1.5, title="Genomic positions", bg="white")
	dev.off()
}

#UPDATE 2022-10-06: We now want to evaluate the whole set of mutations
load("Genomes-VAMIPsVsGenomes_chkpt1.Rdata")
mut_VOC <- read.table("genome_VOC_completeness/2022-09-28_Main_mutations_per_VOC-filtered-rel.tsv", sep="\t",header=T, skip=0, comment.char='',quote="",fill=F,row.names=1, check.names=FALSE) # Load input table
mut_Lin <- read.table("genome_variant_completeness/2022-09-28_Main_mutations_per_variant-filtered-rel.tsv", sep="\t",header=T, skip=0, comment.char='',quote="",fill=F,row.names=1, check.names=FALSE) # Load input table
# apply(mut_VOC,2,function(x){sum(x>=0.7)}) # Get the total expected mutations per VOC at this confidence level
#   Delta Omicron     519  Others   Gamma   Alpha     222      XB
#      36      51      21       4      36      33      10      34
# apply(mut_Lineage,2,function(x){sum(x>=0.7)}) # Same for
#       AY.20   B.1.1.519      BA.1.1       AY.26      AY.100     BA.1.15
#          39          21          56          35          41          55
#        AY.3     B.1.1.7        BA.1      P.1.17      AY.113      AY.103
#          39          33          52          39          41          38
#   B.1.1.222         P.1          XB       AY.44     B.1.429       AY.39
#          10          36          34          39          22          44
#     B.1.621       B.1.2     B.1.427       AY.62   B.1.617.2       AY.25
#          32          14          18          38          24          37
#      AY.122     B.1.243     B.1.634     AY.25.1        C.37         B.1
#          38           7          34          39          30           4
#       AY.43       A.2.5     B.1.632     B.1.635     BA.1.17       B.1.1
#          39          30          46          22          62           8
#   BA.1.17.2   B.1.1.432     B.1.526      P.1.15     B.1.631 BA.1.others
#          63          14          23          39          34          45
# BA.2.others   AY.others
#          74          25
main_mut_map_VOC <- apply(mut_VOC,2,function(x){out=x[x>=0.7];out[order(names(out))]}) # Extract a map of mutations per variant (at least 70% of all genomes from each variant should have them). This holds how common that mutation is in the variant. The last part just orders by position
main_mut_map_Lin <- apply(mut_Lin,2,function(x){out=x[x>=0.7];out[order(names(out))]})
df[is.na(df)] <- 0 # We need to remove NAs before any subsetting
dfVOC <- df[sort(unique(unlist(lapply(main_mut_map_VOC,names)))),] # For processing, we will create two subsets of the main mutation table, one for VOCs (163 unique items for any VOCs at 70%)
# dim(dfVOC)
# [1]   163 17926
dfLin <- df[sort(unique(unlist(lapply(main_mut_map_Lin,names)))),] # and one for all variants with >=20 genomes in the set (unique mutations seen in at least 70% in any variant)
rownames(dfVOC) <- sub("\\|.*","",rownames(dfVOC)) # Remove tags for delta and omicron as we are expanding the VOCs we are observing
rownames(dfLin) <- sub("\\|.*","",rownames(dfLin))
nt2AA <- read.table("MutNtxAA_complete.tsv", sep="\t",header=T, skip=0, comment.char='',quote="",fill=F,row.names=1, check.names=FALSE) # Load nt to AA name table
nt2AA[nt2AA[,"MutAA"]=="NC","MutAA"] <- paste0(nt2AA[nt2AA[,"MutAA"]=="NC","Gene"],":NC") # add a category to flag UTRs
rownames(dfVOC) <- paste(rownames(dfVOC),nt2AA[rownames(dfVOC),"MutAA"], sep='|') # and append this to the mutation names
rownames(dfLin) <- paste(rownames(dfLin),nt2AA[rownames(dfLin),"MutAA"], sep='|')
# dim(dfLin)
# [1]   438 17926
rm(df)
# now, build some a mask based on the positions, each will hold only the variant's core mutations that were selected

strip <- function(string){ # Takes a string object and strips rare characters
	string <- iconv(string,from="UTF-8",to="ASCII//TRANSLIT") # Remove accents and rare chars
# 	string <- tolower(string) # change to lower
	string <- gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", string, perl=TRUE) # Remove multiple spaces
	string <- gsub("\\+","",string) # remove + characters
	string <- gsub("\\|", "_-_",string)
	string <- gsub(" ", "_",string)
	string <- gsub("\\(", "_",string)
	string <- gsub("\\)", "_",string)
	string <- gsub("\\/", "_",string)
	string <- gsub("__", "_",string)
	return(string)
}
FreqMut_perGenome <- function(i){ # For each genome
	dir.create("Main_mut_Freq-VOC", showWarnings=FALSE)
	dir.create("Main_mut_Freq-Lineage", showWarnings=FALSE)
	print(i) # Start processing the ith genome
	name <- colnames(dfVOC)[i] # This holds the name of the genome
	allFreqsVOC <- dfVOC[i] # Exctact the frequencies of all mutations detected in that genome (it is important to preserve it as a dataframe to keep rownames
	allFreqsLin <- dfLin[i]
	freqVOC <- lapply(lapply(main_mut_map_VOC,names), function(x){round(allFreqsVOC[x,1],2)}) # Now, use the VOC/lineage map to extract the corresponding values
	freqLin <- lapply(lapply(main_mut_map_Lin,names), function(x){round(allFreqsLin[x,1],2)})
	out_tableVOC <- matrix("",ncol=length(main_mut_map_VOC), nrow=nrow(dfVOC)); colnames(out_tableVOC) <- names(main_mut_map_VOC); rownames(out_tableVOC) <- rownames(dfVOC) # and prepare an output tables containing the values in a square matrices
	out_tableLin <- matrix("",ncol=length(main_mut_map_Lin), nrow=nrow(dfLin)); colnames(out_tableLin) <- names(main_mut_map_Lin); rownames(out_tableLin) <- rownames(dfLin)
	for(var in names(main_mut_map_VOC)){out_tableVOC[names(main_mut_map_VOC[[var]]),var] <- freqVOC[[var]]} # fill the matrices with the corresponding items
	out_tableVOC <- out_tableVOC[,c("222","519","Alpha","Gamma","Delta","Omicron","XB")] # Order by VOC appearance and remove others
	for(var in names(main_mut_map_Lin)){out_tableLin[names(main_mut_map_Lin[[var]]),var] <- freqLin[[var]]}
	out_tableLin <- out_tableLin[,sort(colnames(out_tableLin))]
	completenessVOC <- apply(out_tableVOC, 2, function(x){vect=as.numeric(x); vect=vect[!is.na(vect)]; sum(vect>0)/length(vect)}) # Calculate how complete each variant would be
	completenessLin <- apply(out_tableLin, 2, function(x){vect=as.numeric(x); vect=vect[!is.na(vect)]; sum(vect>0)/length(vect)})
	print(round(completenessVOC,2))
	print(round(completenessLin,2))
	if(sum(completenessVOC>=0.80)>=2){write.table(out_tableVOC, paste0("Main_mut_Freq-VOC/", strip(name), ".tsv"), sep="\t", quote=FALSE, row.names=T, col.names=NA)}
	if(sum(completenessLin>=0.80)>=2){write.table(out_tableLin, paste0("Main_mut_Freq-Lineage/", strip(name), ".tsv"), sep="\t", quote=FALSE, row.names=T, col.names=NA)}
}

save.image("Genomes-VAMIPsVsGenomes_chkpt3.Rdata")
load("Genomes-VAMIPsVsGenomes_chkpt3.Rdata")
# i <- 4527 # test
for(i in 1:ncol(SNPs))





















# Now, we need to collate each table where mutations are summed by position. We only want to keep substitutions and indels separated
all_mut <- as.data.frame(cbind("id"=1:nrow(df),"Pos"=sub(":.*","",rownames(df)),"MutType"=rep("_Sub", nrow(df))));rownames(all_mut) <- rownames(df) # we start with the full collection of mutation names, and the type of mutations
# now, fill the indels accordingly to identify the type of mutation
all_mut[grep(">\\+",rownames(df)),3] <- "-Ins"
all_mut[grep(">-",rownames(df)),3] <- "-Del"
SNPs <- rowsum(SNPs, group=paste0(all_mut[,2],all_mut[,3])) # Collate by the postition and type of mutation
dim(SNPs)
# [1]  6762 17926
VAMIPs <- rowsum(VAMIPs, group=paste0(all_mut[,2],all_mut[,3])) # Collate by the postition and type of mutation
dim(VAMIPs)
# [1]  6762 17926
SNPs[SNPs>1] <- 1; VAMIPs[VAMIPs>1] <- 1 # Due to rounding, some items may have values slightly larger than 1 (<1e-5), thus, we will truncate those to 1s

id_Sub <- grep("Sub",rownames(SNPs))
id_Del <- grep("Del",rownames(SNPs))
id_Ins <- grep("Ins",rownames(SNPs))
