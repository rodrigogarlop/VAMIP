# ### Functions ###
find_items <- function(queryStrings, dataf, colname){ # IMPORTANT. This is used to avoid partial index matching. Use it to hash a table and use rownames to match a vector. colname may be numeric
	searchVect <- dataf[,colname]
	names(searchVect) <- rownames(dataf)
	out <- searchVect[unlist(queryStrings)]
	return(out)
}

# ### Preprocess data ###
setwd("/home/rod/Documents/01_Projects/SARS/VAMIP/")
df <- read.table("06_VAMIP_tables/All_mut_NoSingl.tsv", header=F, sep ='\t', stringsAsFactors = FALSE,fill=T)
dim(df)
# [1] 1212142       8
df[!complete.cases(df),] # Some line may have NAs, these should be ignored (this is not ideal but it is like adding an N, basically)
#                            V1    V2 V3 V4 V5       V6 V7 V8
# 346310 L030_202101206551_S199 21973  T  + NA 0.333333  3  3
# 346311                      0    NA       NA       NA NA NA
df <- df[complete.cases(df),]
dim(df)
# [1] 1212140       8
names(df) <- c("Genome", "Position", "In", "Out", "Alt_Freq", "DP", "Ref_rev", "Alt_rev") # Rename ivar columns
# We will now remove the _SXXX part of the name (this was kept until now to avoid repeated items but there should be removed at this point if any are still present)
temp_names <- unique(df[,"Genome"]) # First derreplicate genome names as they are (with the library number)
length(temp_names)
# [1] 18404 # These are the total items that are considered for the analysis
# TEST START
# temp_names[18405]="L003_202101082315_S101"
# test <- tail(temp_names[grep("^S", temp_names)]);test # Some tests were carried out to make sure we only remove the last part
# [1] "S6736_S72" "S6741_S81" "S6737_S73" "S6655_S34" "S6625_S3"  "S6735_S40"
# sub("_S[0-9]*$","", test)
# "S6736" "S6741" "S6737" "S6655" "S6625" "S6735"
# dummy=temp_names[18239];dummy
# [1] "L024_S8000_S72"
# sub("_S[0-9]*$","", dummy)
# [1] "L024_S8000"
# TEST END
temp_names_prefix <- sub("_S[0-9]*$","", temp_names)
dup_for_removal <- names(which(table(temp_names_prefix)>1))
remove_list <- temp_names[as.vector(unlist(sapply(dup_for_removal,function(x){grep(x,temp_names)})))]
row_for_removal <- as.vector(unlist(sapply(remove_list,function(x){grep(x,df[,"Genome"])})))
dim(df[row_for_removal,]) # Total items removed
# [1] 1770    8
write.table(df[row_for_removal,], "07_metadata/RepeatedGenomes_vamips_removed.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)
df <- df[-row_for_removal,]
dim(df)
# [1] 1210466       8
# Now that there are no genomes that have the same name (not repeated but faulty called the same), we can replace the _SXX suffix
df[,"Genome"] <- sub("_S[0-9]*$","", df[,"Genome"])
# write.table(df, "07_metadata/Fixed_All_mut_NoSingl.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)


df[,"Position"] <- sprintf("%05d", df[,"Position"]) # add 0s to fix position sorting
df[,"DP_Alt"] <- round(df[,"DP"]*df[,"Alt_Freq"]) # add the depth of the alternative vars
df[,"DP_Ref"] <- df[,"DP"]-df[,"DP_Alt"] # and the depth of reference (original) vars
df[,"Mutation"] <- paste0(df[,"Position"],":",df[,"In"],"->",df[,"Out"]) # This will hold mutations in a new nomenclature such as 00021:C->T with position:ref->alt
df <- df[order(df[,"Position"]),] # now, sort observations by position
meta <- read.table("07_metadata/Metadata_wID.tsv", header=T, sep ='\t',stringsAsFactors = FALSE, check.names=F, fill=T, quote="") # Load metadata with id folio for xref
dim(meta)
# [1] 19028    46
temp <- names(table(meta[,"ID Folio"])[table(meta[,"ID Folio"])>1]) # Create a dictionary with all posible items (folio ID) there are in the meta table
temp
# [1] "L011_202101008588" "L035_2840" # These two are repeated and will thus be removed
repeated <- rep(TRUE,length(temp)); names(repeated) <- temp # create a vector for a hash where keys are repeated folio IDs
bad <- meta[!is.na(repeated[meta[,"ID Folio"]]),] # extract those that are repeated
bad <- bad[order(bad[,"ID Folio"]),] # sort by the id
# and export it
write.table(bad, "07_metadata/RepeatedFOLIO_metadata.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)
meta <- meta[is.na(repeated[meta[,"ID Folio"]]),] # Remove repeated items
dim(meta)
# [1] 19024    46
rownames(meta) <- meta[,"ID Folio"] # and use the folio id as row identifier
meta[,"Rename"] <- paste0(meta[,"Fecha de recoleccion"],"|", meta[,"ID Folio"],"|",meta[,"Linaje Pangolin"],"(",meta[,"Clado Nexstrain"],"_",meta[,"Variants"],")","|",meta[,"Edad"],"|",meta[,"Genero"],"|",meta[,"Estado"]) # append a new column with a date|name|PANGO(Clade)|age|gender|State
temp <- meta[,"Rename"]; names(temp) <- rownames(meta) # IMPORTANT: create an auxiliary vector to avoid partial matches (df has partial string matching)
# Now, match both sets

df[,"Names"] <- temp[df[,"Genome"]] # add names now
bad <- df[is.na(df[,"Names"]),] # Get those that have no metadata available
nrow(bad) # Mutations with no matching genome
# [1] 94669
nonmatched_Genomes <- unique(bad[,"Genome"])
length(nonmatched_Genomes <- unique(bad[,"Genome"])) # Non-matching genomes (these have no metadata)
# [1] 1187
write.table(bad, "07_metadata/Mutations_withNO_metadata.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)
write.table(nonmatched_Genomes, "07_metadata/Genomes_withNO_metadata.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)
# NOTE TEMPORARY FIX: This was donde to avoid going back to previous steps as the metadata table was missing some lab identifiers from the ID Folio column in batches 40 and 41 (and some other error in batch 27). This is only required once, so comment it afterwards
all_names_nosuf <- meta[,"ID Folio"]
test <- sapply(all_names_nosuf,function(x){grep(x,nonmatched_Genomes)})
recovered_nonmatched <- cbind(names(unlist(test)),nonmatched_Genomes[unlist(test)])
write.table(recovered_nonmatched, "07_metadata/Genomes_withNO_metadata-recovered.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)
rename_vect <- recovered_nonmatched[,2]
names(rename_vect) <- recovered_nonmatched[,1]
meta2 <- meta
meta2[!is.na(rename_vect[meta2[,"ID Folio"]]),"ID Folio"] <- rename_vect
rownames(meta2) <- meta2[,"ID Folio"] # and use the folio id as row identifier
meta2[,"Rename"] <- paste0(meta2[,"Fecha de recoleccion"],"|", meta2[,"ID Folio"],"|",meta2[,"Linaje Pangolin"],"(",meta2[,"Clado Nexstrain"],"_",meta2[,"Variants"],")","|",meta2[,"Edad"],"|",meta2[,"Genero"],"|",meta2[,"Estado"])
temp <- meta2[,"Rename"]; names(temp) <- rownames(meta2)
df[,"Names"] <- temp[df[,"Genome"]] # add names now
bad <- df[is.na(df[,"Names"]),] # Get those that have no metadata available
nrow(bad) # Mutations with no matching genome
# [1] 31520
nonmatched_Genomes <- unique(bad[,"Genome"])
length(nonmatched_Genomes <- unique(bad[,"Genome"])) # Non-matching genomes (these have no metadata)
# [1] 445 # This is after the new recovery step (optional)
write.table(bad, "07_metadata/Mutations_withNO_metadata_2ndchance.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)
write.table(nonmatched_Genomes, "07_metadata/Genomes_withNO_metadata_2ndchance.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)
# NOTE END OF TEMPORARY FIX


df <- df[!is.na(df[,"Names"]),] # Now, filter those that could not be identified
dim(df)
# [1] 1178946      12

df[,"Pango"] <- find_items(df[,"Genome"],meta, "Linaje Pangolin") # add PANGO lineage
df[,"Clado"] <- find_items(df[,"Genome"],meta, "Clado Nexstrain") # add patient status
df[,"Date"] <- find_items(df[,"Genome"],meta, "Fecha de recoleccion") # add date now
df[,"State"] <- find_items(df[,"Genome"],meta, "Estado") # add state
df[,"Status"] <- find_items(df[,"Genome"],meta, "Estatus del paciente") # add patient status
write.table(meta, "07_metadata/Metadata_all.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)
xrefmut <- read.table("07_metadata/UNIQUE_def_mutations_DeltaOmicron.tsv", header=F, sep ='\t', stringsAsFactors = FALSE, row.names=1, check.names=F) # Load unique mutation information
df[,"UniqueMut"] <- find_items(df[,"Mutation"],xrefmut, 1) # add key mutation info
xrefnames <- read.table("07_metadata/IMPORT_All_found_xref_with_names_fixed.tsv", header=F, sep ='\t', stringsAsFactors = FALSE, row.names=1, check.names=F) # Load unique name information based on Folio with Selene's selection ### REPLACE UPDATED FILE HERE
df[,"Names_selene"] <- find_items(df[,"Genome"],xrefnames, 1) # Append original names
df[,"Clade"] <- find_items(df[,"Genome"],xrefnames, 2) # as well as lineages
df[,"Obs"] <- find_items(df[,"Genome"],xrefnames, 3) # and flag those that may be recombinant or should be removed
dim(df)
# [1] 1178946      21
length(grep("ista", df[,"Obs"])) # Search those flagged for removal ("Lista Quitar" or "Lista quitar")
# [1] 318
df <- df[-grep("ista", df[,"Obs"]),] # and remove them
dim(df)
# [1] 1178946      21

# ### Filters ###
df[grep("\\+",df[,"Out"]),"Alt_rev"] <- 1 # This offset is a small fix to avoid loosing all indels (since they always have 0 in the Alt_rev column (only shown in the Fwd strand)
df[grep("\\-",df[,"Out"]),"Alt_rev"] <- 1
# temp <- df[,"Ref_rev"]/df[,"DP_Ref"]; temp[is.na(temp)] <- 0
temp <- df[,"Alt_rev"]/df[,"DP_Alt"]; temp[is.na(temp)] <- 0
temp <- as.logical((temp > 0) * (temp < 1))
sum(!temp) # How many have either all observations on the F or R strand?
# [1] 104405
sum(df[,"Alt_Freq"] < 0.05) # How many are not in at least 0.05 of observations?
# [1] 173116
df[!temp,"Alt_Freq"] <- 0.000001 # Those with no reads in both strands are marked as 0 frequency + an offset
df[df[,"Alt_Freq"] < 0.05,"Alt_Freq"] <- 0.000001 # Those will have there values set to 0 + an offset
sum(df[,"DP"] < 20) # How many have less than 20 reads for depth
# [1] 83950
df[df[,"DP"] < 20,"Alt_Freq"] <- 0.000001 # Same for those with low DP
sum(df[,"Alt_Freq"]==0.000001) # How many were marked as bad?
# [1] 267875
sum(df[,"Alt_Freq"]>0.000001) # How many passed all filters
# [1] 911071

write.table(df, "Alt_Whole_table.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)

# ### Extract tables ###
# First, frequency
genomes <- names(table(df[,"Names"])) # Extract the names based on Selene's table
# genomes <- genomes[order(sub(".*\\|","",genomes))] # Sort them by lineage
genomes <- sort(genomes) # sort by date
mutations <- names(table(df[,"Mutation"])) # and extract the mutations
out <- matrix(NA, ncol=length(genomes), nrow=length(mutations)); colnames(out) <- genomes; rownames(out) <- mutations # Create a void matrix container
for(i in 1:nrow(df)){ # and survey those that are present (the rest will be NAs)
	out[df[i,"Mutation"],df[i,"Names"]] <- df[i,"Alt_Freq"]
# 	print(i)
# 	print(df[i,"Mutation"])
# 	print(df[i,"Genome"])
# 	print(out[df[i,"Mutation"],df[i,"Genome"]])
}
dim(out)
# [1] 46877 17927
# Rename the mutations to include which ones are unique to delta/omicron
temp <- find_items(rownames(out),xrefmut, 1); temp[!is.na(temp)] <- paste0("|",temp[!is.na(temp)]);temp[is.na(temp)] <- ""
rownames(out) <- paste0(mutations,temp)
# Create a mask for removing rare items
mask <- out; mask[is.na(mask)] <- 0; mask[mask<0.03] <- 0; mask <- (mask>0)*1 # Add an additional binary table (presence/absence) for those that were included with an offset before
out <- out[rowSums(mask)>=5,] # and remove those from the output table
dim(out) # How many remain?
# [1]  7340 17927
write.table(out, "Alt_variant_Freq.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)

# Next, depth
genomes2 <- names(table(df[,"Names"])) # Extract the names based on Selene's table
# genomes2 <- genomes2[order(sub(".*\\|","",genomes2))] # Sort them by lineage
genomes2 <- sort(genomes2) # Sort them by date
mutations2 <- names(table(df[,"Mutation"])) # and extract the mutations2
out2 <- matrix(NA, ncol=length(genomes2), nrow=length(mutations2)); colnames(out2) <- genomes2; rownames(out2) <- mutations2 # Create a void matrix container
for(i in 1:nrow(df)){ # and survey those that are present (the rest will be NAs)
	out2[df[i,"Mutation"],df[i,"Names"]] <- df[i,"DP"]
# 	print(i)
# 	print(df[i,"Mutation"])
# 	print(df[i,"Genome"])
# 	print(out2[df[i,"Mutation"],df[i,"Genome"]])
}
dim(out2)
# [1] 46877 17927
# Rename the mutations2 to include which ones are unique to delta/omicron
temp <- find_items(rownames(out2),xrefmut, 1); temp[!is.na(temp)] <- paste0("|",temp[!is.na(temp)]);temp[is.na(temp)] <- ""
rownames(out2) <- paste0(mutations2,temp)
# Create a mask for removing rare items
out2 <- out2[rowSums(mask)>=5,] # and remove those from the out2put table
dim(out2) # How many remain?
# [1]  7340 17927
write.table(out2, "Alt_variant_TotalDepth.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)

# And now, with alt allelic depth
genomes2 <- names(table(df[,"Names"])) # Extract the names based on Selene's table
genomes2 <- genomes2[order(sub(".*\\|","",genomes2))] # Sort them by lineage
mutations2 <- names(table(df[,"Mutation"])) # and extract the mutations2
out2 <- matrix(NA, ncol=length(genomes2), nrow=length(mutations2)); colnames(out2) <- genomes2; rownames(out2) <- mutations2 # Create a void matrix container
for(i in 1:nrow(df)){ # and survey those that are present (the rest will be NAs)
	out2[df[i,"Mutation"],df[i,"Names"]] <- df[i,"DP_Alt"]
# 	print(i)
# 	print(df[i,"Mutation"])
# 	print(df[i,"Genome"])
# 	print(out2[df[i,"Mutation"],df[i,"Genome"]])
}
dim(out2)
# [1] 46877 17927
# Rename the mutations2 to include which ones are unique to delta/omicron
temp <- find_items(rownames(out2),xrefmut, 1); temp[!is.na(temp)] <- paste0("|",temp[!is.na(temp)]);temp[is.na(temp)] <- ""
rownames(out2) <- paste0(mutations2,temp)
# Create a mask for removing rare items
out2 <- out2[rowSums(mask)>=5,] # and remove those from the out2put table
dim(out2) # How many remain?
# [1]  7340 17927
write.table(out2, "Alt_variant_AltAllelicDepth.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)

# ### Allele prevalence table ###
counts <- rowSums(mask[rowSums(mask)>4,]) # Create a table for counting how many samples have each mutation
ncol(mask) # How many total observations are there?
# [1] 17927
counts <- cbind("Total_genomes_w/var"=counts, "Total_genomes_w/var %"=round(counts*100/ncol(mask),2)) # Now create a raw and % table
deltas <- mask[,grep("elta", genomes)] # now extract only deltas
count_delta <- rowSums(deltas[rowSums(deltas)>4,])
ncol(deltas) # How many deltas are there?
# [1] 8740
count_delta <- cbind("Delta_genomes_w/var"=count_delta, "Delta_genomes_w/var %"=round(count_delta*100/ncol(deltas),2)) # Now create a raw and % table
temp <- merge(counts, count_delta, by="row.names", all=TRUE); rownames(temp) <- rownames(counts); temp <- temp[,-1]; counts <- temp
omicrons <- mask[,grep("micro", genomes)] # now extract only omicrons
count_omicron <- rowSums(omicrons[rowSums(omicrons)>4,])
ncol(omicrons) # How many omicrons are there?
# [1] 3585
count_omicron <- cbind("Omicron_genomes_w/var"=count_omicron, "Omicron_genomes_w/var %"=round(count_omicron*100/ncol(omicrons),2)) # Now create a raw and % table
temp <- merge(counts, count_omicron, by="row.names", all=TRUE); rownames(temp) <- rownames(counts); temp <- temp[,-1]; counts <- temp
write.table(counts, "Alt_variants_per_genome_Delt_Omi.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)
save.image("explore_vamips_end_worskspace.RData")
load("explore_vamips_end_worskspace.RData")



