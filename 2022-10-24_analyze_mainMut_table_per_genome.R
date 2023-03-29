indir <- "/home/rod/Documents/01_Projects/SARS/VAMIP/10_Analysis_per_genome/Main_mut_Freq-VOC/"
setwd(indir)
outdir <- "Freq_analysis"
dir.create(outdir,showWarnings=FALSE)

bords=c("hotpink", "chartreuse3", "coral1", "purple", "cornflowerblue", "firebrick", "tan2") # This scheme is for VOCs
cols=c("pink", "darkolivegreen2", "lightsalmon", "slateblue1", "darkslategray3", "indianred1", "khaki1") # This scheme is for VOCs
file_list <- list.files(pattern="\\.tsv")
# infile <- "2021-06-08_-_L034_202101118198_-_B.1.1.7_20I_501Y.V1_Alpha_-_035_-_M_-_State_of_Mexico.tsv" # for testing
for(file in file_list){
	print(file)
	infile <- file
	df <- read.table(infile, sep="\t",header=T, skip=0, comment.char='',quote="",fill=F,row.names=1, check.names=FALSE) # Load input table
	nam <- gsub("\\.tsv", "", gsub("_-_", "|", infile))
	ignore <- c("00241", "03037", "14408", "23403") # These are always present in VOCs, they add nothing of value to comparisons
	ignore <- sapply(ignore,function(x){grep(x,rownames(df))}) # We get which rows have them
	df <- df[-ignore,] # and we remove them
	comp <- df[nrow(df),] # store completeness
	exp <- df[nrow(df)-2,] # as well as the expected
	obs <- df[nrow(df)-1,] # and observed mutations per variant
	df <- df[1:(nrow(df)-3),] # ignore observed and total mutations, keep only actual mutation frequency
	varList <- apply(df, 2, function(x){x[!is.na(x)]})
	pdf(paste0(outdir, "/",sub("\\.tsv", "\\.pdf", infile)))
	boxplot(las=1, varList, boxwex=as.numeric(comp), main=nam, cex.main=0.8, ylab="Key Mutations' Frequency", col=cols, border=bords, ylim=c(0,1), outline=FALSE)
# 	boxplot(las=1,lapply(varList, function(x) quantile(x, seq(0,1,0.1))), boxwex=as.numeric(comp), main=nam)
	stripchart(varList, vertical = TRUE, method = "jitter",pch = 2, add = TRUE, col = "gray30", cex=0.5, jitter = 0.1)
	mtext("Box width depicts variant completeness")
	dev.off()
}

# library(vioplot)
# vioplot(varList, col=cols,boxwex=as.numeric(comp))


# Repeat for variants:
indir <- "/home/rod/Documents/01_Projects/SARS/VAMIP/10_Analysis_per_genome/Main_mut_Freq-Lineage/"
setwd(indir)
outdir <- "Freq_analysis"
dir.create(outdir,showWarnings=FALSE)
colScheme <- read.table("/home/rod/Documents/01_Projects/SARS/VAMIP/00_docs/Variant_boxplot_cols.tsv", sep="\t",header=T, skip=0, comment.char='',quote="",fill=F,row.names=1, check.names=FALSE) # Load input table
# bords=c("hotpink", "chartreuse3", "coral1", "purple", "cornflowerblue", "firebrick", "tan2") # This scheme is for VOCs
# cols=c("pink", "darkolivegreen2", "lightsalmon", "slateblue1", "darkslategray3", "indianred1", "khaki1") #
file_list <- list.files(pattern="\\.tsv")
# infile <- "2021-06-08_-_L034_202101118198_-_B.1.1.7_20I_501Y.V1_Alpha_-_035_-_M_-_State_of_Mexico.tsv" # for testing
for(file in file_list){
	print(file)
	infile <- file
	df <- read.table(infile, sep="\t",header=T, skip=0, comment.char='',quote="",fill=F,row.names=1, check.names=FALSE) # Load input table
	df <- df[,rownames(colScheme)] # sort by the color scheme
	nam <- gsub("\\.tsv", "", gsub("_-_", "|", infile))
	ignore <- c("00241", "03037", "14408", "23403") # These are always present in VOCs, they add nothing of value to comparisons
	ignore <- sapply(ignore,function(x){grep(x,rownames(df))}) # We get which rows have them
	df <- df[-ignore,] # and we remove them
	comp <- df[nrow(df),] # store completeness
	exp <- df[nrow(df)-2,] # as well as the expected
	obs <- df[nrow(df)-1,] # and observed mutations per variant
	df <- df[1:(nrow(df)-3),] # ignore observed and total mutations, keep only actual mutation frequency
	varList <- apply(df, 2, function(x){x[!is.na(x)]})
	pdf(paste0(outdir, "/",sub("\\.tsv", "\\.pdf", infile)), height=4, width=10)
	par(oma=c(1.2,1,0,0))
	boxplot(las=2, varList, boxwex=as.numeric(comp), main=nam, cex.main=0.8, ylab="Key Mutations' Frequency", col=colScheme[,1], border=colScheme[,2], ylim=c(0,1), outline=FALSE)
# 	boxplot(las=1,lapply(varList, function(x) quantile(x, seq(0,1,0.1))), boxwex=as.numeric(comp), main=nam)
	stripchart(varList, vertical = TRUE, method = "jitter",pch = 2, add = TRUE, col = "gray30", cex=0.5, jitter = 0.1)
	mtext("Box width depicts variant completeness")
	dev.off()
}
