library("pheatmap")
library("stringr")
save_pheatmap_pdf <- function(x, filename, width, height) { # This time, adjust the size of the output pdf file to a smaller size
  pdf(filename,width,height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
# ### Version 1: Allelic frequency ###
df <- read.table("Alt_variant_Freq.tsv", sep="\t",header=T, skip=0, comment.char='',quote="",fill=F,row.names=1, check.names=FALSE) # Load input table
dim(df)
# [1]  7340 17927
df <- df[,order(colnames(df))]
minimetatab <- t(as.data.frame(str_split(colnames(df), "\\|")))
rownames(minimetatab) <- NULL # row names are removed
genomes <- minimetatab[,2]
xrefnames <- read.table("07_metadata/IMPORT_MLReferenceTable.tsv", header=T, sep ='\t', stringsAsFactors = FALSE, row.names=1, check.names=F) # Load unique mutation information
xrefnames <- xrefnames[xrefnames[,3]!="",1:3] # Remove those with no folio
rownames(xrefnames) <- xrefnames[,"Folio"] # and use folios are keys
xrefnames["newname"] <- paste(xrefnames[,1],xrefnames[,2],sep='-')
temp <- xrefnames[minimetatab[,2],4]
colnames(df)[!is.na(temp)] <- paste0(colnames(df)[!is.na(temp)],"||",temp[!is.na(temp)])

# Now select only those flagged by selene as multilineages
keep_genomes <- grep("\\|ML",colnames(df))
multlin <- df[,keep_genomes]
multlin[multlin<0.03] <- NA
keep_mutations <- rowSums((!is.na(multlin))*1)!=0
multlin <- multlin[keep_mutations,]
dim(multlin)
# [1] 330  20
write.table(multlin,"ML_alt_freq.tsv", quote=FALSE, sep='\t', row.names=TRUE, col.names=NA)


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
groups[locations>=266 & locations<=13441] <- "02 ORF1a"
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
