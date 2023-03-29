# Started: 2021-11-03
# by Rodrigo García-López for Carlos Arias's Viromics Group at IBt, UNAM, Cuernavaca, Mexico as part of the CoViGen-Mex SARS-CoV-2 survillance in Mexico.
# Under GNU GPLv3 license
# Disclamer: This script was written for use with specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.
# This script was tested with R version 4.1.2 (2021-11-01) -- "Bird Hippie"

# This script is intended for subsetting an input table based on the initial and end dates (spanish format dd/mm/YYYY), and exctract only matching dates (both start and ending dates included).
# The input a 4-column table containing the following items: Accession ID	Collection date	Lineage	Location
# Only the Collection date is considered for processing (it should be in position 2). This was provided with SPA format: dd/mm/YYYY
# Both target dates should be written in universal date format: YYYY-mm-dd to avoid confusion.

# Execute in bash as follows:
# Rscript extract_dates.R <input_table> <start_date> <end_date> <prefix_output>

# Test in  R:
# intable="01_GISAID_Tables/2021-10-28_All_items_for_areaplots.tsv"
# intable="01_GISAID_Tables/2021-10-28_First_year_for_areaplots.tsv"
# intable="01_GISAID_Tables/2021-10-28_First_year_for_areaplots.tsv"
# date_s="2020-01-01"
# date_f="2020-05-01"
# prefix="02_Split_dates/2021-10-28_First_year"

# Test in bash:
# Rscript extract_dates.R 01_GISAID_Tables/2021-10-28_All_items_for_areaplots.tsv 2020-01-01 2020-03-01 02_Split_dates/First_year

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) { # at least, 4 arguments should be included
  stop("A minimum of 4 arguments are mandatory: Rscript extract_dates.R <input_table> <start_date> <end_date> <prefix_output>", call.=FALSE)
}
intable <- as.character(args[1]) # Get a string handle for the input table
date_s <- as.character(args[2]) # Get a string handle for the starting date
date_f <- as.character(args[3]) # Get a string handle for the ending date
prefix <- as.character(args[4]) # Get a string handle for the output (it can contain a complete path but requires a seed name for the output file)
library(areaplot)
### Functions ###
printable_dates <- function(vector){ # Gets a vector with ordered dates (fomat must have %d at the end), get the position where each month stats or at day 15 and adjust them accordingly for plotting
	dates <- unique(sort(c(grep("01$|15$",vector),length(vector),1)))
	if(length(dates)>20){dates <- unique(sort(c(grep("01$",vector),length(vector),1)))}
	if((dates[length(dates)]-dates[length(dates)-1])<10){dates <- dates[-(length(dates)-1)]} # Fix dates showing
	if((dates[2]-dates[1])<10){dates <- dates[-2]}
	return(dates)
}
define_plot_scheme <- function(){
	lty <- c("solid","dashed","dotted","dotdash","longdash","twodash")
	col <- c('chartreuse3', 'cornflowerblue', 'darkgoldenrod1', 'peachpuff3','mediumorchid2', 'turquoise3', 'wheat4','slategray2',"black","coral1","aquamarine2","blue2","violetred2","palegreen3","purple3","magenta1","limegreen","darkorange2","darkgray")
	lwd <- c(1.3,2)
	pars <- expand.grid(col = col, lty = lty, lwd = lwd, stringsAsFactors = FALSE) # This will create all
	return(pars)
}
plot_dates_lines <- function(mat,scheme,name,total){ # Gets an input matrix with rownames containing continuous date. This returns a plot for the current graphical device. The second parameter should be a 3-column table containing color names at col1, linetype at col2, and lwd at col3 as produced by define_plot_scheme(). A string defines the name. Total contains a number of total items for printing in the X axis
	if(length(grep("Mexico",colnames(mat))>1)){mat=t(fill_missing_items(t(mat),alllocations))} # Fix missing States if present
	oripar <- par(no.readonly=TRUE) # save default params
	cex.X=0.53
	dates <- printable_dates(rownames(mat)) # Determine where each month starts
	par(oma = c(1, 1, 1, 4)) # This is just a creative fix to plot the legend outside
	with(scheme[1:ncol(mat),],matplot(mat[,1:ncol(mat)], type = 'l', col = col, lty = lty, lwd = lwd, main = paste0(name," by date"), ylab = NA, xlab = NA, las = 2, xaxt = 'n'))
	mtext(paste("(n = ",total," total observations)"))
	axis(1, las=2, at = dates, labels = rownames(mat)[dates], cex.axis = 0.75)
	par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0.5), mar = c(0, 0, 0, 0), new = TRUE)
	plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
	legend("right", legend = colnames(mat), xpd = TRUE, horiz = FALSE, inset = c(0,0), bty = "n", col = scheme$col[1:ncol(mat)], lty = scheme$lty[1:ncol(mat)], lwd = scheme$lwd[1:ncol(mat)], cex = cex.X)
	par(oripar) # reset parameters
	return(name)
}
fill_missing_items <- function(mat,vector){ # Gets a table and a vector with all rows that should be present, outputs an expanded row collection with missing dates for complete calendar in that range. rownames should have date format as %Y-%m-%d
	xtable <- as.data.frame(matrix(0,nrow=length(vector),ncol=ncol(mat)), stringsAsFactors = FALSE) # Create empty vessel for output
	rownames(xtable) <- vector # use the vector as rownames
	colnames(xtable) <- colnames(mat) # and inherit the names
	invisible(sapply(rownames(mat),function(x) {xtable[x,] <<- mat[x,]}))	# append the original values in the corresponding places (write to higher env variable
	return(xtable)
}
alllocations <- c("Baja California", "Baja California Sur", "Sonora", "Sinaloa", "Chihuahua", "Durango", "Coahuila", "Nuevo Leon", "Tamaulipas", "Zacatecas", "San Luis Potosi", "Aguascalientes", "Guanajuato", "Queretaro", "Hidalgo", "State of Mexico", "Mexico City", "Puebla", "Tlaxcala", "Morelos", "Nayarit", "Jalisco", "Colima", "Michoacan", "Veracruz", "Guerrero", "Oaxaca", "Tabasco", "Chiapas", "Campeche", "Yucatan", "Quintana Roo")
rolling_N_avg <- function(mat,int){ # Input should be a table with continuous data at rows and the desired interval for the mean. If an even number is provided, the average will be placed one position to the right of as there is no single middle number
	before <- after <- trunc(int/2) # initialize with same range before and after each position
	if(int%%2==0){after <- after-1} # If even, shift the upper half of the range by 1 position (the mean will be calculated for the next value next to the middle as there is no exact number in it)
	roll <- sapply((before+1):(nrow(mat)-after),function(y) {apply(mat,2, function(x) mean(as.numeric(x[(y-before):(y+after)])))})
	colnames(roll) <- row.names(mat)[(before+1):(nrow(mat)-after)]
	roll <- t(roll)
	return(roll)
}
date_vs_X <- function(mat,string,int) { # Input should have at least a "Date" column with format as %Y-%m-%d, a target vector of dates that should be present (passed as a string vector) and an integer for the days in the rolling average
	dateX <- twoway_table(mat[,c("Collection.date",string)])
	dates <- seq(as.Date(rownames(dateX)[1]),as.Date(rownames(dateX)[nrow(dateX)]),by="day") # create the complete date range
	dateX <- fill_missing_items(dateX,dates) # use the predicted missing days to get the whole date spectrum (adding 0s when required)
	dateX <- rolling_N_avg(dateX,int) # smoothen with a N-day average
	return(dateX)
}
twoway_table <- function(mat){ # Cross two variables from a 2 column matrix (var1, var2), return a table of var2 x var1 with n columns depending on total var2 items
	xtable <- xtabs(rep(1,nrow(mat))~mat[,1]+mat[,2], data=mat)
	return(xtable)
}
plot_dates_areas <- function(mat,scheme,name,total){ # Gets an input matrix with rownames containing continuous date. This returns a plot for the current graphical device. The second parameter should be a 3-column table containing color names at col1, linetype at col2, and lwd at col3 as produced by define_plot_scheme(). A string defines the name. Total contains a number of total items for printing in the X axis
	if(length(grep("Mexico",colnames(mat))>1)){mat=t(fill_missing_items(t(mat),alllocations))} # Fix missing States if present
	oripar <- par(no.readonly=TRUE) # save default params
	cex.X <- 0.53
	dates <- printable_dates(rownames(mat)) # Determine where each month starts
	par(oma = c(1, 1, 1, 4)) # This is just a creative fix to plot the legend outside
	cols <- rev(scheme[1:ncol(mat),1])
	areaplot(mat[,ncol(mat):1],col=cols,ylim=c(0,100),border=NA, las=2, main=paste0(name," by date"), ylab=NA, xaxt='n', xlab=NA)
	mtext(paste("(n = ",total," total observations)"))
	axis(1, las=2, at = dates, labels = rownames(mat)[dates], cex.axis = 0.75)
	par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0.5), mar = c(0, 0, 0, 0), new = TRUE)
	plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
	legend("right", legend = colnames(mat), xpd = TRUE, horiz = FALSE, inset = c(0,0), bty = "n", col = rev(cols), pch=15,cex = cex.X)
	par(oripar) # reset parameters
	return(name)
}
append_region <- function(mat){
	states <- data.frame("States"=sort(unique(mat[,4])), stringsAsFactors = FALSE)
	rownames(states) <- states[,1] # use ages as rownames
	# Now, manually edit the table based on Blanca's classification
	states["Baja California",1] <- "1:NW"
	states["Baja California Sur",1] <- "1:NW"
	states["Sonora",1] <- "1:NW"
	states["Sinaloa",1] <- "1:NW"
	states["Chihuahua",1] <- "1:NW"
	states["Durango",1] <- "1:NW"
	states["Coahuila",1] <- "2:NE"
	states["Nuevo Leon",1] <- "2:NE"
	states["Tamaulipas",1] <- "2:NE"
	states["Zacatecas",1] <- "3:CN"
	states["San Luis Potosi",1] <- "3:CN"
	states["Aguascalientes",1] <- "3:CN"
	states["Guanajuato",1] <- "3:CN"
	states["Queretaro",1] <- "3:CN"
	states["Hidalgo",1] <- "4:CS"
	states["State of Mexico",1] <- "4:CS"
	states["Mexico City",1] <- "4:CS"
	states["Puebla",1] <- "4:CS"
	states["Tlaxcala",1] <- "4:CS"
	states["Morelos",1] <- "4:CS"
	states["Nayarit",1] <- "5:W"
	states["Jalisco",1] <- "5:W"
	states["Colima",1] <- "5:W"
	states["Michoacan",1] <- "5:W"
	states["Veracruz",1] <- "6:S"
	states["Guerrero",1] <- "6:S"
	states["Oaxaca",1] <- "6:S"
	states["Tabasco",1] <- "6:S"
	states["Chiapas",1] <- "6:S"
	states["Campeche",1] <- "7:SE"
	states["Yucatan",1] <- "7:SE"
	states["Quintana Roo",1] <- "7:SE"
	states["aaaa",1] <- "8:SE"
	mat[,5] <- states[mat[,4],1] # Now append it to the main dataframe
	colnames(mat)[5]="Region"
	return(mat)
}

### MAIN ###
df <- read.table(intable, header=T, sep ='\t', stringsAsFactors = FALSE) # Load table (header is expected)
df[,2] <- as.Date(df[,2], "%d/%m/%Y")
# Edos <- xtabs(rep(1,nrow(df))~df[,2]+df[,4], data=df)
xtra <- append_region(df)
pdf(paste(prefix,"sub.pdf",sep="_"),width=14)
	Edos <- date_vs_X(xtra,"Location",1)
	scheme <- define_plot_scheme()
	pb <- barplot(las=1,rowSums(Edos),col="darkgray",border=NA,xaxt='n',main="Total sampled genomes per date")
	axis(1, las=2, at=pb[printable_dates(rownames(Edos))], labels=rownames(Edos)[printable_dates(rownames(Edos))], cex.axis = 0.75)
	mtext(paste("(n = ",nrow(xtra)," total observations)"))
	Edos <- date_vs_X(xtra,"Location",7)
	plot_dates_lines(Edos, scheme, "States",nrow(xtra))
	plot_dates_areas(prop.table(Edos,1)*100, scheme, "States",nrow(xtra))
	Regions <- date_vs_X(xtra,"Region",7)
	plot_dates_lines(Regions, scheme, "Regions",nrow(xtra))
	plot_dates_areas(prop.table(Regions,1)*100, scheme, "Regions",nrow(xtra))
dev.off()

subtime <- df # Start with a copy
subtime <- subset(df, Collection.date >= date_s & Collection.date <= date_f)
write.table(subtime, paste(prefix,"from",date_s,"to",date_f,"sub.tsv",sep="_"), sep="\t", quote=FALSE, row.names=F, col.names=T)
xtra <- append_region(subtime)
pdf(paste(prefix,"from",date_s,"to",date_f,"sub.pdf",sep="_"),width=14)
	Edos <- date_vs_X(xtra,"Location",1)
	scheme <- define_plot_scheme()
	pb <- barplot(las=1,rowSums(Edos),col="darkgray",border=NA,xaxt='n',main="Total sampled genomes per date")
	axis(1, las=2, at=pb[printable_dates(rownames(Edos))], labels=rownames(Edos)[printable_dates(rownames(Edos))], cex.axis = 0.75)
	mtext(paste("(n = ",nrow(xtra)," total observations)"))
	Edos <- date_vs_X(xtra,"Location",7)
	plot_dates_lines(Edos, scheme, "States",nrow(xtra))
	plot_dates_areas(prop.table(Edos,1)*100, scheme, "States",nrow(xtra))
	Regions <- date_vs_X(xtra,"Region",7)
	plot_dates_lines(Regions, scheme, "Regions",nrow(xtra))
	plot_dates_areas(prop.table(Regions,1)*100, scheme, "Regions",nrow(xtra))
dev.off()
print("End of execution")
