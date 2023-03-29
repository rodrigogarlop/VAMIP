# Started: 2022-07-15
# by Rodrigo García-López for Carlos Arias's Viromics Group at IBt, UNAM, Cuernavaca, Mexico as part of the CoViGen-Mex SARS-CoV-2 survillance in Mexico.
# Under GNU GPLv3 license
# Disclamer: This script was written for use with specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.

# The script is used with the full collection of ivar tables where the Wuhan-Hu-1 genome (MN908947.3) gff was included for SARS-CoV-2 ORF identification so that nt could be linked to AA mutations. For the input, multiple ivar tables (one per genome) were contcatenated.

# ### DEFINE FUNCTIONS ###
gene_by_pos <- function(int){ # From a genomic position, return the actual gene and the actual starting position within the gene
# 	st <- c(0, 265, 13468, 21563, 25393, 26245, 26523, 27199, 27394, 27756, 27894, 28274, 29558, 100000)
	st <- c(1, 266, 13468, 21556, 21563, 25385, 25393, 26221, 26245, 26473, 26523, 27192, 27202, 27388, 27394, 27756, 27888, 27894, 28260, 28274, 29534, 29558, 29675, 100000) # This vector holds each interval (all ORFs). Endings were adjusted with +1 for calculations
	genes <- c("5pUTR", "ORF1a", "ORF1b", "NC", "S", "NC", "ORF3a", "NC", "E", "NC", "M", "NC", "ORF6", "NC", "ORF7a", "ORF7b", "NC", "ORF8", "NC", "N", "NC", "ORF10", "3pUTR") # This one holds the corresponding gene per interval, including genomic outer regions (as non-coding).
	x <- findInterval(int,st) # Extract the bin it corresponds to
	start <- st[x] #  This will hold the exact position of the starting gene
	end <- round(st[x+1])-1 # And this, the end of the gene
	gene <- genes[x] # This holds the name of the gene
# 	p_st <- g_st <- "NC"
# print(gene)
# 	if(gene!="NC") {
		g_st <- int-start+1 # Get the starting nt position within the gene
		p_st <- trunc((g_st+2)/3) # Get the starting AA position within the protein
# 	}
	out <- list(start,end,gene,g_st,p_st)
	names(out) <- c("NtSt","NtEnd","Gene","GSt","PSt")
	return(out)
}
fix_indels <- function(inmat){ # For indels, we can remap the missing genes if any
	del <- grep("^-", inmat[,"Nt_alt"]) # Search indels by their -/+ sign
	ins <- grep("^\\+", inmat[,"Nt_alt"])
	indels <- sort(c(ins,del)) # get a single vector for both
	indels_found <- gene_by_pos(as.numeric(inmat[indels,"PosNt"])) # Use a custom function to extract the gene and initial positions (nt and AA)
	inmat[indels,"Gene"] <- indels_found$Gene # Append gene
	inmat[indels,"MutAA"] <- paste0(inmat[indels,"Gene"],":",nchar(inmat[indels,"Nt_alt"])-1,";",indels_found$PSt) # add the mutation (a placeholder is used for ins and del type)
	inmat[indels,"GenePosNt"] <- indels_found$GSt
	inmat[indels,"GenePosAA"] <- indels_found$PSt
	residue <- ((nchar(inmat[indels,"Nt_alt"])-1)%%3) # Calculate modulo for detecting frameshifts
	inmat[indels,"MutType"] <- ifelse(residue==0,"_NFS","_FS") # determine if frameshifts are present
	inmat[del,"MutType"] <- paste0("Del",inmat[del,"MutType"]) # and append either del or ins accodingly
	inmat[ins,"MutType"] <- paste0("Ins",inmat[ins,"MutType"])
	inmat[del,"MutAA"] <- sub(";","Del",inmat[del,"MutAA"]) # and the same for the mutations
	inmat[ins,"MutAA"] <- sub(";","Ins",inmat[ins,"MutAA"])
	ignore <- which(as.logical((inmat[,"Gene"]=="5pUTR")+(inmat[,"Gene"]=="3pUTR")+(inmat[,"Gene"]=="NC"))) # determine which are in the UTR ends to remove any coding info
	inmat[ignore,"MutType"] <- "NC"
	inmat[ignore,"MutAA"] <- "NC"
	inmat[ignore,"GenePosNt"] <- "NC"
	inmat[ignore,"GenePosAA"] <- "NC"
	return(inmat)
}
# Gene map (reference only, this wasn't actually used)
# groups[locations<=265] <- "01 5-UTR"
# groups[locations>=266 & locations<=13467] <- "02 ORF1a"
# groups[locations>=13468 & locations<=21555] <- "03 ORF1b" # 13442-13467 have no assignation
# groups[locations>=21563 & locations<=25384] <- "04 S"
# groups[locations>=25393 & locations<=26220] <- "05 ORF3a"
# groups[locations>=26245 & locations<=26472] <- "06 E"
# groups[locations>=26523 & locations<=27191] <- "07 M"
# groups[locations>=27202 & locations<=27387] <- "08 ORF6"
# groups[locations>=27394 & locations<=27759] <- "09 ORF7a"
# groups[locations>=27756 & locations<=27887] <- "10 ORF7b"
# groups[locations>=27894 & locations<=28259] <- "11 ORF8"
# groups[locations>=28274 & locations<=29533] <- "12 N"
# groups[locations>=28283 & locations<=28573] <- "13 ORF9b"
# groups[locations>=29558 & locations<=29674] <- "14 ORF10"
# groups[locations>=29675 & locations<=29903] <- "15 3-UTR"

# ### MAIN ###
df <- read.table(gzfile("single_ivar_table.tsv.gz"), sep='\t', stringsAsFactors = FALSE, header = TRUE, comment.char = "", check.names = FALSE) # This is a concatenated file of all ivar tables in the set, they must have AA info
dim(df) # Each row is a separate point mutation, even if it is a VAMIP (not part of the consensus)
# [1] 1497685      20
df <- df[df[,"REGION"]!="REGION",] # Since multiple tables were included, we should first remove the headers
dim(df)
# [1] 1477784      20
df[is.na(df[,"GFF_FEATURE"]),"GFF_FEATURE"] <- "NC" # Change missing values to NC
# we should ignore those with NA in AA or in UTR regions (as these are not translated)
MutNt <- apply(df, 1, function(x){paste0(sprintf("%05d",as.numeric(x["POS"])),":",x["REF"],"->",x["ALT"])})
# Build a map of initial positions of each gene
# genomeMap <- c(0, 266, 13468, 21563, 25393, 26245, 26523, 27199, 27394, 27756, 27894, 28274, 28283, 29558, 29675,100000)+1 # Original positions
genomeMap <- c(100000, 266, 13468, 21563, 25393, 26245, 26523, 27202, 27394, 27756, 27894, 28274, 28283, 29558, 100000,100000)-1 # This was modified to consider those that are non-coding
names(genomeMap) <- c("5pUTR","ORF1a","ORF1b","S","ORF3a","E","M","ORF6","ORF7a","ORF7b","ORF8","N","ORF9b","ORF10","3pUTR", "NC")
init <- as.numeric(df[,"POS"]) - as.numeric(genomeMap[df[,"GFF_FEATURE"]])
# Build an output table with items describing the mutations (in Nt and AA)
out <- cbind("MutNt"=MutNt, "PosNt"=df[,"POS"], "Nt_ref"=df[,"REF"], "Nt_alt"=df[,"ALT"], "Gene"=df[,"GFF_FEATURE"], "GenePosNt"=init, "GenePosAA"=ceiling(init/3), "AA_ref"=df[,"REF_AA"], "AA_alt"=df[,"ALT_AA"])
# Add a new category for the standardized name for AA mutations
MutAA <- apply(out, 1, function(x){paste0(x["Gene"],":",x["AA_ref"],x["GenePosAA"],x["AA_alt"])})
out <- cbind(out, MutAA) # and append it
test <- as.numeric(out[,"GenePosNt"])
out[test<0,"MutAA"] <- "NC" # Label those non coding for AA
out[test<0,"GenePosNt"] <- "NC"
out[test<0,"GenePosAA"] <- "NC"
out[test<0,"AA_ref"] <- "NC"
out[test<0,"AA_alt"] <- "NC"
# head(out)
# Finally, dereplicate items:
vect <- out[,"MutNt"] # Get the total items
allMut <- split(seq_along(vect), vect) # and split them in a list of positions
out <- out[unlist(lapply(allMut,function(x){x[1]})),]
# head(out,100)
out <- out[order(as.numeric(out[,"PosNt"])),]
out <- as.data.frame(out)
out[,"MutType"] <- "NonSyn"
out[out[,"AA_ref"]==out[,"AA_alt"],"MutType"] <- "Syn"
out[out[,"AA_ref"]=="NC","MutType"] <- "NC"
# table(out["MutType"])
# MutType
# NC NonSyn    Syn 
# 13106  52252  16570 
dim(out)
# [1] 81928    11
out2 <- fix_indels(out) # Append information on indels
write.table(out2,"MutNtxAA_complete.tsv", sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)
