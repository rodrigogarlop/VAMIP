# Based on /home/rod/Documents/02_Collaborations/Minireview_variantes/Epi_2022-02-28/COVIGEN/preFilter_raw_metatable.R
# # # # # # # First, load the complete table
# # # # # # library("readxl")
# # # # # # # library("xtable")
# # # # # # xlsx <- "01_data/2022_02_28_DB_MexCov2.xlsx"
# # # # # # df <- as.data.frame(read_excel(xlsx, sheet = 1))
# # # # # # # df <- read.table("01_data/2022_02_07_Metadata_BDMexCov2_FINAL.tsv",header=T, sep='\t', skip=0, comment.char='',fill=T, check.names=FALSE, stringsAsFactors = FALSE)
# # # # # # # df <- read.table(unz("01_data/17Enero2022_Metadata_BDMexCov2_FINAL.zip", "17Enero2022_Metadata_BDMexCov2_FINAL.tsv"), header=T, quote="\"", sep="\t", check.names=FALSE, fill=TRUE)
# # # # # # print("Starting dimensions:")
# # # # # # dim(df)
# # # # # # # [1] 54488    27
# # # # # # # Keep only the columns we want
# # # # # # #  [1] "Accession ID"           "Nombre"                 "Fecha de recolección"
# # # # # # #  [4] "Fecha de envío"         "Estado"                 "Municipio"
# # # # # # #  [7] "Genero"                 "Edad"                   "Estatus del paciente"
# # # # # # # [10] "Linaje Pangolin"        "Clado Nexstrain"        "ID Folio"
# # # # # # # [13] "RdRP"                   "Mutaciones Nucleótidos" "Mutaciones Aminoácido"
# # # # # # # [16] "Deleciones Nucleótidos" "Deleciones Aminoácido"
# # # # # #
# # # # # # df <- df[,c(1,2,5,6,7,8,12,13,14,26,27,20,23,28:31)]
# # # # # # # dim(df)
# # # # # # # [1] 54488    17
#START of alt version
df <- read.table("2022_04_27_DB_MexCov2_v1_xtras.tsv", header=T, sep ='\t', stringsAsFactors = FALSE, fill=TRUE,quote="\"",check.name=FALSE)
dim(df)
# [1] 58737    24
# Keep only the columns we want
#End of alt version
# df <- df[,c(1,2,3,4,5,6,7,8,9,10,13,14,15,18,19,13,14,15,20,21,22,23)]
dim(df)
# [1] 58727    17

# Remove rare punctuations or accents
strip <- function(string){ # Takes a string object and strips rare characters
	string <- iconv(string,from="UTF-8",to="ASCII//TRANSLIT") # Remove accents and rare chars
# 	string <- tolower(string) # change to lower
	string <- gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", string, perl=TRUE) # Remove multiple spaces
	return(string)
}
names(df) <- strip(names(df))
df[,"Estado"] <- strip(df[,"Estado"])
df[,"Municipio"] <- strip(df[,"Municipio"])
df[,"Edad"] <- strip(df[,"Edad"])
# Fix the Ambulatory column
df[df[,"Estatus del paciente"]=="Asymptomatic","Estatus del paciente"] <- "Ambulatory"
df[df[,"Estatus del paciente"]=="Asymptomatic and Ambulatory","Estatus del paciente"] <- "Ambulatory"
df[df[,"Estatus del paciente"]=="Outpatient","Estatus del paciente"] <- "Ambulatory"
df[df[,"Estatus del paciente"]=="Mild","Estatus del paciente"] <- "Ambulatory"
df[df[,"Estatus del paciente"]=="Symptomatic and Ambulatory","Estatus del paciente"] <- "Ambulatory"
df[df[,"Estatus del paciente"]=="Symptomatic and Hospitalized","Estatus del paciente"] <- "Hospitalized"
df[df[,"Estatus del paciente"]=="Released","Estatus del paciente"] <- "Hospitalized"
df[df[,"Estatus del paciente"]=="Fatal","Estatus del paciente"] <- "Deceased"
df[df[,"Estatus del paciente"]=="Live","Estatus del paciente"] <- "unknown"
df[df[,"Estatus del paciente"]=="Moderate","Estatus del paciente"] <- "unknown"
df[df[,"Estatus del paciente"]=="Symptomatic","Estatus del paciente"] <- "unknown"
# # # # Add a new column to emulate the old input table for recycling the script Filter_Delta.R
# # # df["FLAG"] <- df[,"Linaje Pangolin"]=="B.1.617.2"
# # # df[grep("^AY\\.",df[,"Linaje Pangolin"]),"FLAG"] <- TRUE
# # # df[df[,"FLAG"],"ID"] <-  df[df[,"FLAG"],"Accession ID"]
# # # df[df[,"FLAG"],"New_lineage"] <-  df[df[,"FLAG"],"Linaje Pangolin"]
# # # df[is.na(df)] <- ""
# # # df <- df[,-which(names(df)=="FLAG")]
print("Ending dimensions:")
dim(df)
# [1] 58727    24
write.table(df, "07_metadata/preFiltered.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
