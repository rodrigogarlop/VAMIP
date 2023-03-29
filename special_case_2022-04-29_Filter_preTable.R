df <- read.table("07_metadata/preFiltered.tsv", header=T, quote="\"", sep="\t", check.names=FALSE)
# df <- read.table("01_data/Xai_tabla_filtrada_Delta_fixed5.tsv",header=T, sep='\t', skip=0, comment.char='',fill=FALSE, check.names=FALSE, stringsAsFactors = FALSE)
print("Starting total items:")
dim(df)
# [1] 58727    17
# Now, test how many were different
df[,"ID"] <- TRUE # Just to keep an old format intact (nor really required)
df[,"New_Lineage"] <- df[,"Linaje Pangolin"]
df[!df[,"ID"],"New_Lineage"] <- df[!df[,"ID"],"Linaje Pangolin"] # New lineage was originally bound to Blanca's table. This is just legacy for namesake
df[,"Delta"] <- grepl("^AY\\.",df[,"New_Lineage"]) # Flag those that are delta variants
df[df[,"New_Lineage"]=="B.1.617.2","Delta"] <- TRUE # And the original Maharashtra variant
df[,"Omicron"] <- grepl("^BA\\.",df[,"New_Lineage"]) # Flag those that are omicron variants
df[df[,"New_Lineage"]=="B.1.1.529","Omicron"] <- TRUE # And the original Maharashtra variant
df[,"Recomb"] <- grepl("^X.",df[,"New_Lineage"]) # Flag those that are known recombinants

# Change new column to categories, switch the labels and change to characters
df[,"Delta"] <- as.factor(df[,"Delta"]);levels(df[,"Delta"]) <- c("Others","Delta"); df[,"Delta"] <- as.character(df[,"Delta"])
df[,"Omicron"] <- as.factor(df[,"Omicron"]);levels(df[,"Omicron"]) <- c("Others","Omicron"); df[,"Omicron"] <- as.character(df[,"Omicron"])
df[,"Recomb"] <- as.factor(df[,"Recomb"]);levels(df[,"Recomb"]) <- c("Others","Recomb"); df[,"Recomb"] <- as.character(df[,"Recomb"])
sum(df[,"Delta"]=="Delta") # Total Delta
# [1] 25026
# Create subcategories
df[,"Delta_sub"] <- df[,"Delta"]
sort(table(df[df[,"Delta_sub"]=="Delta","New_Lineage"]),decreasing=T)
#      AY.20      AY.26     AY.100       AY.3     AY.113     AY.103      AY.44
#      11413       5490       3083       1479        945        572        304
#      AY.62  B.1.617.2     AY.122    AY.25.1      AY.39      AY.25      AY.43
#        250        229        215        155        155        153        107
#      AY.13      AY.47     AY.118      AY.75   AY.119.2       AY.4     AY.119
#         83         51         31         24         20         20         19
#      AY.15   AY.122.4      AY.38      AY.42      AY.74     AY.127     AY.117
#         16         15         14         14         13         12         11
#      AY.14      AY.54   AY.119.1     AY.114     AY.121       AY.2    AY.39.2
#         11         11          8          7          7          7          7
#     AY.3.1      AY.53       AY.5    AY.98.1   AY.116.1     AY.101      AY.33
#          6          6          5          5          4          3          3
#    AY.99.2     AY.105     AY.126      AY.35     AY.4.2     AY.4.9      AY.46
#          3          2          2          2          2          2          2
#    AY.46.4      AY.48      AY.52   AY.120.1   AY.121.1     AY.124 AY.124.1.1
#          2          2          2          1          1          1          1
#     AY.125      AY.27     AY.3.3      AY.32      AY.34    AY.34.1    AY.34.2
#          1          1          1          1          1          1          1
#    AY.39.1    AY.4.10    AY.43.8       AY.6      AY.64      AY.66      AY.67
#          1          1          1          1          1          1          1
#      AY.73      AY.77     AY.9.2      AY.98
#          1          1          1          1
# I'd rather add them manually
df[df[,"Delta_sub"]=="Delta","Delta_sub"] <- "Other Delta" # First, change the name
df[df[,"New_Lineage"]=="AY.20","Delta_sub"] <- "AY.20"
df[df[,"New_Lineage"]=="AY.26","Delta_sub"] <- "AY.26"
df[df[,"New_Lineage"]=="AY.3","Delta_sub"] <- "AY.3"
df[df[,"New_Lineage"]=="AY.113","Delta_sub"] <- "AY.113"
df[df[,"New_Lineage"]=="AY.103","Delta_sub"] <- "AY.103"
df[df[,"New_Lineage"]=="AY.100","Delta_sub"] <- "AY.100"
df[df[,"New_Lineage"]=="AY.4","Delta_sub"] <- "AY.44"
df[df[,"New_Lineage"]=="AY.62","Delta_sub"] <- "AY.62"
df[df[,"New_Lineage"]=="B.1.617.2","Delta_sub"] <- "B.1.617.2"
df[df[,"New_Lineage"]=="AY.122","Delta_sub"] <- "AY.122"
df[df[,"New_Lineage"]=="AY.25.1","Delta_sub"] <- "AY.25"
df[df[,"New_Lineage"]=="AY.39","Delta_sub"] <- "AY.39"
df[df[,"New_Lineage"]=="AY.25","Delta_sub"] <- "AY.25" # This is joined with AY.25.1
df[df[,"New_Lineage"]=="AY.43","Delta_sub"] <- "AY.43"
table(df[,"Delta_sub"])
#      AY.100      AY.103      AY.113      AY.122       AY.20       AY.25
#        3083         572         945         215       11413         153
#     AY.25.1       AY.26        AY.3       AY.39       AY.43       AY.44
#         155        5490        1479         155         107          20
#       AY.62   B.1.617.2 Other Delta      Others
#         250         229         760       33701

sum(df[,"Omicron"]=="Omicron") # Total Omicron
# [1] 13316
df[,"Omicron_sub"] <- df[,"Omicron"]
sort(table(df[df[,"Omicron_sub"]=="Omicron","New_Lineage"]),decreasing=T)
#    BA.1.1   BA.1.15      BA.1   BA.1.17      BA.2 BA.1.17.2    BA.2.9 BA.1.1.15
#      8059      2350      1990       181       149       117        94        58
#   BA.1.18  BA.1.1.1    BA.2.3   BA.1.13 BA.1.1.14 BA.1.1.16   BA.1.14    BA.1.5
#        56        41        36        35        26        25        19        17
#  BA.1.1.2  BA.1.1.8 BA.1.1.10 BA.1.15.1 BA.1.1.11    BA.2.1   BA.1.19    BA.1.4
#        11         8         7         6         4         4         3         3
#    BA.1.7   BA.2.12  BA.1.1.6   BA.1.20   BA.1.16    BA.1.2    BA.1.6    BA.1.9
#         3         3         2         2         1         1         1         1
#   BA.2.10 BA.2.12.1  BA.2.3.2
#         1         1         1
df[grep("^BA.1",df[,"New_Lineage"]), "Omicron_sub"] <- "BA.1"
df[grep("^BA.2",df[,"New_Lineage"]), "Omicron_sub"] <- "BA.2"
table(df[,"Omicron_sub"])
#   BA.1   BA.2 Others
#  13027    289  45411
df[,"Variants"] <- df[,"Delta"] # make one more column for general VOCs and VOIs (start with this one as it has the most items
df[df[,"New_Lineage"]=="B.1.1.7","Variants"] <- "Alpha"
df[grep("^Q",df[,"New_Lineage"]), "Variants"] <- "Alpha"
df[grep("^P.1",df[,"New_Lineage"]), "Variants"] <- "Gamma"
df[df[,"New_Lineage"]=="B.1.1.519","Variants"] <- "519"
df[df[,"Omicron"]=="Omicron","Variants"] <- "Omicron"
df[df[,"Recomb"]=="Recomb","Variants"] <- "Recomb"

# Now, we'll be adding regions
edos <- read.table("07_metadata/EstadoRegion_v4.tsv",header=T, sep='\t', skip=0, comment.char='',fill=FALSE, check.names=FALSE, stringsAsFactors = FALSE, row.names=1)
df <- cbind(df, edos[df[,"Estado"],]); rownames(df) <- NULL # Append region groupings
# Fix dates
df[,"Fecha de recoleccion"] <-  as.Date(df[,"Fecha de recoleccion"],"%d/%m/%Y")
df[,"Fecha de envio"] <-  as.Date(df[,"Fecha de envio"],"%d/%m/%Y")
# Remove all those newer than the last 16 nov date
print("Last date in the whole table:")
test_date <- max(as.Date(df[,"Fecha de recoleccion"]));test_date
# [1] "2022-04-13"
# Fix gender
df[,"Genero"] <- sub("[Ff]emale","F",df[,"Genero"])
df[,"Genero"] <- sub("[Mm]ale","M",df[,"Genero"])
df[,"Genero"] <- sub("[Uu]nknown","u",df[,"Genero"])
# Fix age
df[,"Edad"] <- sub("[Uu]nknown","u",df[,"Edad"])
df[,"Edad"] <- sub("Sin informacion","u",df[,"Edad"])
# DEPRECATED: START
# # # # Fix identifiers
# # # df[,"ID Folio"] <- sub("Sin información","u",df[,"ID Folio"])
# # # df[,"SINAVE ID"] <- sub("Sin información","u",df[,"SINAVE ID"])
# # # df[,"SINOLAVE ID"] <- sub("Sin información","u",df[,"SINOLAVE ID"])
# DEPRECATED: END
# Fix type of patient
df[,"Estatus del paciente"] <- sub("[Aa]mbulatory","Amb",df[,"Estatus del paciente"])
df[,"Estatus del paciente"] <- sub("[Hh]ospitalized","Hosp",df[,"Estatus del paciente"])
df[,"Estatus del paciente"] <- sub("[Dd]eceased","Dec",df[,"Estatus del paciente"])
df[,"Estatus del paciente"] <- sub("[Uu]nknown","Unk",df[,"Estatus del paciente"])
# Append the month
df[,"Month"] <- format(as.Date(df[,"Fecha de recoleccion"]),"%Y-%m(%b)")
pdf("07_metadata/2022-04-28_Weekly_genomes.pdf")
	barplot(las=2,table(df[,"Month"]), border=F, col="cornflowerblue")
dev.off()
# Append the Year
df[,"Year"] <- format(as.Date(df[,"Fecha de recoleccion"]),"%Y")
# Append week
	# To do this, first create a calendar (week 1 starts on sunday of the first complete week in the year and reaches up to 52)
	week <- data.frame("week"=c(paste0("20W",sprintf('%0.2d', rep(1:52,each=7))),paste0("21W",sprintf('%0.2d', rep(1:52,each=7))),paste0("22W",sprintf('%0.2d', rep(1:52,each=7)))))
	# Now rename the rows to use them as a hash (dictionary)
	rownames(week) <- seq(as.Date("2020/01/05"), as.Date("2022/12/31"), by="day")[1:nrow(week)]
	# write.table(week, "week_calendar.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
	# Now used the newly created dictionary to append the corresponding week
	df[,"Week"] <- week[as.character(df[,"Fecha de recoleccion"]),1]
# Next, append age categories
	# First create dictionary
	age1 <- data.frame("Age"=as.character(rep(paste0(sprintf('%0.2d',1:14),":",sapply(seq(0,130,10),function(x) paste0(x,"-",x+9))),each=10)), stringsAsFactors = FALSE)
	rownames(age1) <- 0:(nrow(age1)-1) # and use age1s as rownames
	age1[as.numeric(rownames(age1))>=100,1] <- "11:100+" # This was later added to collate all 100 an over
	df[,"Age_range_by10"] <- age1[as.character(df[,"Edad"]),1] # Now append it to the main dataframe
	# Repeat with smaller bins
	age2 <- data.frame("Age"=rep(paste0(sprintf('%0.2d',1:28),":",sapply(seq(0,135,5),function(x) paste0(x,"-",x+4))),each=5), stringsAsFactors = FALSE)
	rownames(age2) <- 0:(nrow(age2)-1) # and use age2s as rownames
	age2[as.numeric(rownames(age2))>=100,1] <- "21:100+" # This was later added to collate all 100 an over
	df[,"Age_range_by5"] <- age2[as.character(df[,"Edad"]),1] # Now append it to the main dataframe
	# Repeat with a new classification
	age3 <- data.frame("age3"=c(rep("1:0-12", 13), rep("2:13-25",13),rep("3:26-45",20), rep("4:46-60",15), rep("5:61-74",14), rep("6:75+",55)), stringsAsFactors = FALSE)
	rownames(age3) <- 0:(nrow(age3)-1) # and use age3s as rownames
	df[,"Age_range_manual"] <- age3[as.character(df[,"Edad"]),1] # Now append it to the main dataframe
	# Added a new age scheme matching vaccination
	age4 <- data.frame("age3"=c(rep("1:0-17", 18), rep("2:18-29",12),rep("3:30-39",10), rep("4:40-49",10), rep("5:50-59",10), rep("6:60+",100)), stringsAsFactors = FALSE)
	rownames(age4) <- 0:(nrow(age4)-1) # and use age3s as rownames
	df[,"Age_vac"] <- age4[as.character(df[,"Edad"]),1] #
pdf("07_metadata/2022-04-28_age5_genomes.pdf")
	barplot(las=2,table(df[,"Age_range_by10"]),border=F, col="coral1")
	barplot(las=2,table(df[,"Age_range_by5"]),border=F, col="coral1")
	barplot(las=2,table(df[,"Age_range_manual"]),border=F, col="coral1")
	barplot(las=2,table(df[,"Age_vac"]),border=F, col="coral1")
dev.off()
	fix_age <- function(vect){ # Ages may be missing. Fix all and convert them to character.
		vect[vect=="u"]=NA # change "u"s to NAs to prevent warnings
		vect <- sprintf('%0.3d',as.numeric(vect))
		vect[vect=="NA"]="u" # Revert NAs (they are now text) to "u"s
		return(vect)
	}
	df[,"Edad"] <- fix_age(df[,"Edad"])
	Date <- as.Date(df[,"Fecha de recoleccion"], format= "%Y-%m-%d") # Create a vector with all dates
	df[Date > as.Date("2021-02-18", "%Y-%m-%d"),"Last_Vaccinated"] <- "A:60+"
	df[Date > as.Date("2021-05-06", "%Y-%m-%d"),"Last_Vaccinated"] <- "B:50-59"
	df[Date > as.Date("2021-06-09", "%Y-%m-%d"),"Last_Vaccinated"] <- "C:40-49"
	df[Date > as.Date("2021-07-06", "%Y-%m-%d"),"Last_Vaccinated"] <- "D:30-39"
	df[Date > as.Date("2021-08-02", "%Y-%m-%d"),"Last_Vaccinated"] <- "E:18-29"
	df[Date > as.Date("2021-11-29", "%Y-%m-%d"),"Last_Vaccinated"] <- "F:15-17"
print("Resulting total items:")
dim(df)
# [1] 58727    46
write.table(df, "07_metadata/AllVariants.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

print("Output total items:")
nrow(df)
subset <- df[df["ID Folio"]!="",] # Filter only those having a folio ID
write.table(subset, "07_metadata/Metadata_wID.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
