library("ggplot2")
library("ggpubr")
library(dplyr)
library(plyr)

# Test in bash:
# Rscript plot.lineages_v8.R 01_GISAID_Tables/2021-10-28_First_year_for_areaplots.tsv EstadoRegion_v3.txt 03_Lineages/test
# Test in R:
intable="01_GISAID_Tables/2021-10-28_First_year_for_areaplots.tsv"
regions_table="EstadoRegion_v3.txt"
prefix="03_Lineages/test"

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) { # at least, 4 arguments should be included
  stop("A minimum of 3 arguments are mandatory: Rscript extract_dates.R <input_table> <regions_table> <prefix_output>", call.=FALSE)
}
intable <- as.character(args[1]) # Get a string handle for the input table
regions_table <- as.character(args[2]) # Get a string handle for the starting date
prefix <- as.character(args[3]) # Get a string handle for the output (it can contain a complete path but requires a seed name for the output file)

#Tiene dos entradas.
#Archivo separado por tabulador con cuatro columnas: 1) Accession ID, 2)Collection date, 3) Lineage y 4) Location
#linCMX= read.table(file= 'gisaid_hcov-19_1Nov20-18Jun21_Metadata.txt', header=T, sep ='\t' )
#linCMX= read.table(file= 'gisaid_hcov-19_1Dic21-18Jun21_Metadata.txt', header=T, sep ='\t' )
linCMX= read.table(file= intable, header=T, sep ='\t', stringsAsFactors = FALSE )

#Archivo de dos columnas separadas por tabulador: 1)Estado y 2) Region
EstReg= read.table(file= regions_table, header=T, sep ='\t', stringsAsFactors = FALSE )
sort(unique(linCMX$Lineage))

linajes <- c(
             "B.1.1.7",    #Alfa
             "Q.1",
             "Q.3",
             "Q.8",
             
             "B.1.1.519",  #Mexicana
             
             "B.1.427",    #Ex Epsilon
             "B.1.429",    #Ex Epsilon
             
             "B.1.628",
             "B.1.631",
             
             "B.1.351",    #Beta
             
             "B.1.617.2",  #Delta
             
             "AY.2",
             "AY.3",
             "AY.4", 
             "AY.5", 
             
             "AY.7.1",
             "AY.9", 
             
             "AY.10",
             "AY.11",
            
             "AY.12",
             "AY.13",
             "AY.14",
             "AY.15",
             
             "AY.17",
             "AY.19",
             "AY.20",
             "AY.21",
             "AY.23",
             "AY.24",
             "AY.25",
             "AY.26",
             
             "AY.32",
             "AY.33",
             "AY.34",
             "AY.35",
             "AY.37",
             "AY.38",
             
             "P.1",        #Gamma
             "P.1.1", 
             "P.1.2",
             "P.1.3",
             "P.1.4",
             "P.1.7",
             "P.1.8",
             "P.1.9", 
             "P.1.10",
             "P.1.10.2",
             "P.1.12",
             
             "C.37",       # Lambda
             "C.37.1",
             
             "B.1.621",    # Mu
             "B.1.621.1",
             
             ##################
             
             "B.1.526",   #Ex Iota
             
             "B.1.617.1"  #Ex Kappa
             
             
             ##################
             
             #"P.2",
             #"P.3"
             
             )

#Poner la fecha en formato y otra variable que sea por mes
linCMX$Collection.date<- as.Date(linCMX$Collection.date, "%d/%m/%Y")
linCMX$M<-as.factor(format(as.Date(linCMX$Collection.date), "%Y-%m"))
#Poner como "Otros" a todos los regristros que no estan en esta lista
Cambiar<-!(linCMX$Lineage %in% linajes)
linCMX$Lineage[Cambiar]<-"Otros"

#linCMX2 <-linCMX[(linCMX$Lineage != "Otros"),]
#linCMX2 <- na.omit(linCMX2)
#linCMX <- linCMX2
nameLinage<-c(
  "B.1.1.7",    #Alfa
  "Q.1",
  "Q.3",
  "Q.8",
  
  "B.1.1.519",  #Mexicana
  
  "B.1.427",    #Ex Epsilon
  "B.1.429",    #Ex Epsilon
  
  "B.1.628",
  "B.1.631",
  
  "B.1.351",    #Beta
  
  "B.1.617.2",  #Delta
  
  "AY.2",
  "AY.3",
  "AY.4", 
  "AY.5", 
  
  "AY.7.1",
  "AY.9", 
  
  "AY.10",
  "AY.11",
  
  "AY.12",
  "AY.13",
  "AY.14",
  "AY.15",
  
  "AY.17",
  "AY.19",
  "AY.20",
  "AY.21",
  "AY.23",
  "AY.24",
  "AY.25",
  "AY.26",
  
  "AY.32",
  "AY.33",
  "AY.34",
  "AY.35",
  "AY.37",
  "AY.38",
  
  "P.1",        #Gamma
  "P.1.1", 
  "P.1.2",
  "P.1.3",
  "P.1.4",
  "P.1.7",
  "P.1.8",
  "P.1.9", 
  "P.1.10",
  "P.1.10.2",
  "P.1.12",
  
  "C.37",       # Lambda
  "C.37.1",
  
  "B.1.621",    # Mu
  "B.1.621.1",
  
  ##################
  
  "B.1.526",   #Ex Iota
  
  "B.1.617.1",  #Ex Kappa
  
  
  ##################
  
  #"P.2",
  #"P.3",
  "Otros"
  
)

# test <- linCMX
# test <- na.omit(test) 
# 
# test$Lineage[test$Lineage%in% c('B.1.1.7', 'Q.1', 'Q.3')] <- 'Alfa'
# test$Lineage[test$Lineage %in% c('B.1.351', 'B.1.351.1')] <- 'Beta'
# test$Lineage[test$Lineage %in% c("P.1",
#                                  "P.1.1", 
#                                  "P.1.2",
#                                  "P.1.4",
#                                  "P.1.7",
#                                  "P.1.8",
#                                  "P.1.9")] <- 'Gamma'
# test$Lineage[test$Lineage %in% c("B.1.617.2",
#                                  "AY.2", 
#                                  "AY.3", 
#                                  "AY.3.1",
#                                  "AY.4", 
#                                  "AY.5",
#                                  "AY.6",
#                                  "AY.7.1",
#                                  "AY.7.2",
#                                  "AY.9", 
#                                  "AY.10", 
#                                  "AY.11", 
#                                  "AY.12",
#                                  "AY.13",
#                                  "AY.14",
#                                  "AY.15",
#                                  "AY.16",
#                                  "AY.17",
#                                  "AY.18",
#                                  "AY.19",
#                                  "AY.20",
#                                  "AY.21",
#                                  "AY.22",
#                                  "AY.23",
#                                  "AY.24",
#                                  "AY.25")] <- 'Delta'
# test$Lineage[test$Lineage == 'B.1.526'] <- 'Iota'
# test$Lineage[test$Lineage == 'B.1.617.1'] <- 'Kappa'
# test$Lineage[test$Lineage %in% c("C.37",
#                                  "C.37.1")] <- 'Lambda'
# test$Lineage[test$Lineage %in% c("B.1.621",
#                                  "B.1.621.1")] <- 'Mu'



#Cambiar Mexico City a español
Cambiar<-(linCMX$Location %in% c("Mexico City","Distrito Federal"))
linCMX$Location[Cambiar]<-"Ciudad de Mexico"
#Cambiar State of Mexico a español
Cambiar<-(linCMX$Location %in% c("State of Mexico"))
linCMX$Location[Cambiar]<-"Estado de Mexico"
#Poner a cada registro a que región pertenece
linCMX$Region<-linCMX$Location
for (i in 1:32){
  linCMX[linCMX$Location==EstReg[i,1],]$Region<-EstReg[i,2] 
}
FechaInicio<-min(linCMX[,2])
Fecha<-min(linCMX[,2])
# linCMX<-linCMX[(linCMX$Collection.date>FechaInicio & linCMX$Collection.date<Fecha),]
sort(unique(linCMX$Lineage))

#Graficas por regiones
#Titulos<-c("Region NW (BCN, BCS, DUR, CHH, SON, SIN)","Region NE (COA, NLE,TAM)","Region WE (COL, JAL,MIC, NAY)","Region CN (AGU, GUA, QUE, SLP, ZAC)","Region CS (CMX, HID, MEX,MOR, PUE, TLA)","Region SE (CHP, GRO, OAX, TAB, VER)","Region S (CAM, ROO, YUC)")
Titulos<-c("Region NW (BCN, BCS, DUR, CHH, SON, SIN)","Region NE (COA, NLE,TAM)","Region WE (COL, JAL,MIC, NAY)","Region CN (AGU, GUA, QUE, SLP, ZAC)","Region CS (CMX, HID, MEX,MOR, PUE, TLA)", "Region SE (CAM, ROO, YUC)", "Region S (CHP, GRO, OAX, TAB, VER)")
Regiones<-c("NW","NE","WE","CN","CS","SE","S")
#Fecha final por region
FechaInicio<-c(min(linCMX[,2]),min(linCMX[,2]),min(linCMX[,2]),min(linCMX[,2]),min(linCMX[,2]),min(linCMX[,2]),min(linCMX[,2]))
Fecha<-c(min(linCMX[,2]),min(linCMX[,2]),min(linCMX[,2]),min(linCMX[,2]),min(linCMX[,2]),min(linCMX[,2]),min(linCMX[,2]))
#Fecha<-c('2021-08-07','2021-08-09','2021-08-09','2021-08-09','2021-08-19','2021-08-19','2021-08-10')
#Fecha<-c('2021-08-02','2021-08-09','2021-08-09','2021-08-09','2021-08-16','2021-08-16','2021-08-09')
Imag<- list() 
for (i in 1:length(Regiones)){
  linCMX$F<-1
  #Obtener registros de la recion NW
  Variant<-linCMX[linCMX$Region==Regiones[i],]
  Variant$Fac<-factor(Variant$Lineage,levels=nameLinage)
  #Sacar cuantos registros hay por mes en la region
  Ttm<-ddply(Variant, .(M), summarize,  sum(F))
  print (Regiones[i])
  print (Ttm)
  GraphVariant<-ggplot(Variant, aes(Collection.date, after_stat(count),fill = Fac)) +
    geom_density(position = "fill",adjust=2.5) +  #adjust se usa para hacer mas grande el promedio en los valores relativo, quita ruido y picos extraños
    scale_fill_manual(breaks=nameLinage,
                      values=c(
                               "B.1.1.7"="darkorange",
                               "Q.1"	=	 "darkorange1",
                               "Q.3"	=	 "darkorange1",
                               "Q.8"	=	 "darkorange1",
                               
                               "B.1.1.519"="indianred1",
                               
                               "B.1.427"="darkolivegreen2",
                               "B.1.429"="darkolivegreen3",
                               
                               "B.1.619"="aquamarine1",
                               "B.1.628"="azure2",
                               "B.1.631"="azure4",
                               
                               "B.1.351"="purple",
                               
                               "B.1.617.2"="gold3",
                               "AY.2"="lightsalmon1",
                               "AY.3"="lightpink4",
                               "AY.4"="peru",
                               "AY.5"="lightgoldenrod4",
                               
                               "AY.7.1"	=	 "khaki2",
                               "AY.9"	=	 "khaki3",
                               
                               "AY.10"	=	 "lightgoldenrod2",
                               "AY.11"	=	 "lightgoldenrod2",
                               
                               "AY.12"	=	 "lightgoldenrod4",
                               "AY.13"	=	 "lemonchiffon",
                               "AY.14"	=	 "lemonchiffon",
                               "AY.15"	=	 "lemonchiffon",
                               "AY.17"	=	 "lemonchiffon2",
                               "AY.19"	=	 "lemonchiffon2",
                               
                               "AY.20"	=	 "sienna3",
                               
                               "AY.21"="lemonchiffon2",
                               "AY.23"="lemonchiffon2",
                               "AY.24"="lemonchiffon3",
                               "AY.25"="lemonchiffon3",
                               "AY.26"="orangered4",
                               
                               "AY.32"	=	 "lightyellow2",
                               "AY.33"	=	 "lightyellow2",
                               "AY.34"	=	 "lightyellow3",
                               "AY.35"	=	 "lightyellow3",
                               "AY.37"	=	 "lightyellow4",
                               "AY.38"	=	 "lightyellow4",
                               
                               "P.1"="green3",
                               "P.1.1"="green4",
                               "P.1.2"="lightseagreen",
                               "P.1.3"="lightseagreen",
                               "P.1.4"="lightseagreen",
                               "P.1.7"="mediumseagreen",
                               "P.1.8"="mediumseagreen",
                               "P.1.9"="mediumspringgreen",
                               "P.1.10"="mediumspringgreen",
                               "P.1.10.2"="mediumspringgreen",
                               "P.1.12"="greenyellow",
                               "C.37"="blueviolet",
                               "C.37.1"="blueviolet",
                               
                               "B.1.621"="azure4",
                               "B.1.621.1"="azure4",
                               
                               ##################
                               
                               "B.1.526"="steelblue1",
                               
                               "B.1.617.1"="cyan2",
                               
                               ##################
                               
                               "Otros"="magenta"),
                      labels=c(
                         "Alfa (B.1.1.7)",
                         "Alfa (Q.1)",
                         "Alfa (Q.3)",
                         "Alfa (Q.8)",
                         
                                "B.1.1.519**",
                         
                                "B.1.427*",
                                "B.1.429*",
                         
                                "B.1.628",
                                "B.1.631",
                         
                                "Beta (B.1.351)",
                         
                                "Delta (B.1.617.2)",
                                "Delta (AY.2)",
                                "Delta (AY.3)",
                                "Delta (AY.4)",
                                "Delta (AY.5)",
                         "Delta (AY.7.1)",
                                "Delta (AY.9)",
                         "Delta (AY.10)",
                         "Delta (AY.11)",
                                "Delta (AY.12)", 
                                "Delta (AY.13)", 
                                "Delta (AY.14)", 
                                "Delta (AY.15)",
                         "Delta (AY.17)",
                         "Delta (AY.19)",
                                
                                "Delta (AY.20)",
                         "Delta (AY.21)",
                                "Delta (AY.23)",
                                "Delta (AY.24)",
                                "Delta (AY.25)",
                        "Delta (AY.26)",
                        "Delta (AY.32)",
                        "Delta (AY.33)",
                        "Delta (AY.34)",
                        "Delta (AY.35)",
                        "Delta (AY.37)",
                        "Delta (AY.38)",
                        #        
                                "Gamma (P.1)",
                                "Gamma (P.1.1)", 
                                "Gamma (P.1.2)", 
                        "Gamma (P.1.3)", 
                                "Gamma (P.1.4)", 
                                "Gamma (P.1.7)", 
                                "Gamma (P.1.8)", 
                                "Gamma (P.1.9)",
                        "Gamma (P.1.10)",
                        "Gamma (P.1.10.2)",
                        "Gamma (P.1.12)",
                        
                        #        
                                
                                "Lambda (C.37)",
                                "Lambda (C.37.1)",
                        #        
                                "Mu (B.1.621)",
                                "Mu (B.1.621.1)",
                        #        
                        "B.1.526 ***",
                        "B.1.617.1 ***",
                        
                        "Otros"),drop = FALSE )                       
  #scale_fill_manual(breaks=c("B.1.1.7","B.1.1.519","B.1.351","B.1.351.1","B.1.427","B.1.429","B.1.526","B.1.526.2","B.1.617.1","B.1.617.2","C.37","P.1","P.1.1","P.2","P.3","Otros"),labels=c("B.1.1.7","B.1.1.519","B.1.351","B.1.351.1","B.1.427","B.1.429","B.1.526","B.1.526.2","B.1.617.1","B.1.617.2","C.37","P.1","P.1.1","P.2","P.3","Otros") )
  GraphVariant<- ggpar(GraphVariant,font.label = list(size=8, color = "grey27"),xlab="",ylab="Abundancia relativa")+
    theme_bw() + scale_y_continuous(expand = c(0, 0))+
    scale_x_date(date_breaks = "1 month",date_labels="%b/%Y",date_minor_breaks="1 month", limits = as.Date(c(FechaInicio[i],Fecha[i])),expand = c(0,0))+
    theme(panel.border = element_blank(),axis.title=element_text(size=8.5),axis.text= element_text(size=8,margin=margin(c(-1,0,0,-5))))+
    theme( plot.title = element_text(hjust = 0.5, size=9,face = "bold", margin=margin(c(5,0,1,0)))) +
    ggtitle(Titulos[i])+
    theme(panel.grid.major = element_blank()) +
    theme(
      legend.title = element_blank(),
      legend.position="right", 
      legend.margin=margin(),
      legend.box="vertical",
      legend.text = element_text(size = 6.0),
      #legend.box.margin=margin(0,10,-5,10),
      legend.key.size = unit(0.4, "cm"))+ guides(fill=guide_legend(ncol=3))+
    labs(x = paste("n =",nrow(Variant),sep=" ")) +
    theme(plot.margin = unit(c(0,0,0,0.4), "lines")) 
  Imag[[i]]<-GraphVariant
  #tiff(paste(paste("Test/Regiones/", "12Jul2021_VOC_VOI_", sep=""),Regiones[i],".tiff",sep=""), width=17.5, height=12, units="cm", res=300,compression ="lzw")
  pdf(paste(paste(prefix, "_VOC_VOI_", sep=""),Regiones[i],".pdf",sep=""), width=18, height=12, units="cm", res=300)
  print (GraphVariant)
  dev.off()
}
aux <- Imag[[7]]
leg<-get_legend(aux + theme(legend.text = element_text(size = 9)) + guides(fill=guide_legend(ncol=3)))  
for (i in 1:length(Regiones)){
  Imag[[i]]<-Imag[[i]] + theme(legend.position = "none")
}
#tiff("Test/Regiones/12Jul2021_VOC_VOI_Regiones.tiff", width=17.5, height=21, units="cm", res=300,compression ="lzw")
pdf(paste0(prefix,"_VOC_VOI_Regiones.pdf", width=24, height=32, units="cm", res=300)
ggarrange(Imag[[1]],Imag[[2]],Imag[[3]],Imag[[4]],Imag[[5]],Imag[[6]],Imag[[7]],leg,ncol=2,nrow=4,common.legend = FALSE,widths=c(2,2),labels = c("A","B","C","D","E","F","G"),font.label = list(size = 9),hjust=-1.2)
dev.off()

# ##Grafica a nivel nacional
FechaInicio<-min(linCMX[,2])
Fecha<-min(linCMX[,2])
linCMX$F<-1
#Obtener todos los registros
Variant<-linCMX
Variant$Fac<-factor(Variant$Lineage,levels=nameLinage)
#Sacar cuantos registros hay por mes a nivel nacional
Ttm<-ddply(Variant, .(M), summarize,  sum(F))
print (Ttm)
GraphVariant<-ggplot(Variant, aes(Collection.date, after_stat(count),fill = Fac)) +
  geom_density(position = "fill",adjust=2.5) +
  scale_fill_manual(breaks=nameLinage,
                    values=c(
                      "B.1.1.7"="darkorange",
                      "Q.1"	=	 "darkorange1",
                      "Q.3"	=	 "darkorange1",
                      "Q.8"	=	 "darkorange1",
                      
                      "B.1.1.519"="indianred1",
                      
                      "B.1.427"="darkolivegreen2",
                      "B.1.429"="darkolivegreen3",
                      
                      "B.1.628"="azure2",
                      "B.1.631"="azure4",
                      
                      "B.1.351"="purple",
                      
                      "B.1.617.2"="gold3",
                      "AY.2"="lightsalmon1",
                      "AY.3"="lightpink4",
                      "AY.4"="peru",
                      "AY.5"="lightgoldenrod4",
                      
                      "AY.7.1"	=	 "khaki2",
                      "AY.9"	=	 "khaki3",
                      
                      "AY.10"	=	 "lightgoldenrod2",
                      "AY.11"	=	 "lightgoldenrod2",
                      
                      "AY.12"	=	 "lightgoldenrod4",
                      "AY.13"	=	 "lemonchiffon",
                      "AY.14"	=	 "lemonchiffon",
                      "AY.15"	=	 "lemonchiffon",
                      
                      "AY.17"	=	 "lemonchiffon2",
                      "AY.19"	=	 "lemonchiffon2",
                      
                      "AY.20"	=	 "sienna3",
                      
                      "AY.21"="lemonchiffon2",
                      "AY.23"="lemonchiffon2",
                      "AY.24"="lemonchiffon3",
                      "AY.25"="lemonchiffon3",
                      "AY.26"="orangered4",
                      
                      "AY.32"	=	 "lightyellow2",
                      "AY.33"	=	 "lightyellow2",
                      "AY.34"	=	 "lightyellow3",
                      "AY.35"	=	 "lightyellow3",
                      "AY.37"	=	 "lightyellow4",
                      "AY.38"	=	 "lightyellow4",
                      
                      "P.1"="green3",
                      "P.1.1"="green4",
                      "P.1.2"="lightseagreen",
                      "P.1.3"="lightseagreen",
                      "P.1.4"="lightseagreen",
                      "P.1.7"="mediumseagreen",
                      "P.1.8"="mediumseagreen",
                      "P.1.9"="mediumspringgreen",
                      "P.1.10"="mediumspringgreen",
                      "P.1.10.2"="mediumspringgreen",
                      "P.1.12"="greenyellow",
                      
                      "C.37"="blueviolet",
                      "C.37.1"="blueviolet",
                      
                      "B.1.621"="azure4",
                      "B.1.621.1"="azure4",
                      
                      ##################
                      
                      "B.1.526"="steelblue1",
                      
                      "B.1.617.1"="cyan2",
                      
                      ##################
                      
                      "Otros"="magenta"),
                    labels=c(
                      "Alfa (B.1.1.7)",
                      "Alfa (Q.1)",
                      "Alfa (Q.3)",
                      "Alfa (Q.8)",
                      
                      "B.1.1.519**",
                      
                      "B.1.427*",
                      "B.1.429*",
                      
                      "B.1.628",
                      "B.1.631",
                      
                      "Beta (B.1.351)",
                      
                      "Delta (B.1.617.2)",
                      "Delta (AY.2)",
                      "Delta (AY.3)",
                      "Delta (AY.4)",
                      "Delta (AY.5)",
                      "Delta (AY.7.1)",
                      "Delta (AY.9)",
                      "Delta (AY.10)",
                      "Delta (AY.11)",
                      "Delta (AY.12)", 
                      "Delta (AY.13)", 
                      "Delta (AY.14)", 
                      "Delta (AY.15)",
                      "Delta (AY.17)",
                      "Delta (AY.19)",
                      "Delta (AY.20)",
                      "Delta (AY.21)",
                      "Delta (AY.23)",
                      "Delta (AY.24)",
                      "Delta (AY.25)",
                      "Delta (AY.26)",
                      
                      "Delta (AY.32)",
                      "Delta (AY.33)",
                      "Delta (AY.34)",
                      "Delta (AY.35)",
                      "Delta (AY.37)",
                      "Delta (AY.38)",
                      #        
                      "Gamma (P.1)",
                      "Gamma (P.1.1)", 
                      "Gamma (P.1.2)",
                      "Gamma (P.1.3)",
                      "Gamma (P.1.4)", 
                      "Gamma (P.1.7)", 
                      "Gamma (P.1.8)", 
                      "Gamma (P.1.9)",
                      "Gamma (P.1.10)",
                      "Gamma (P.1.10.2)",
                      "Gamma (P.1.12)",
                      
                      #        
                      
                      "Lambda (C.37)",
                      "Lambda (C.37.1)",
                      #        
                      "Mu (B.1.621)",
                      "Mu (B.1.621.1)",
                      #        
                      "B.1.526 ***",
                      "B.1.617.1 ***",
                      
                      "Otros"),drop = FALSE )  
#scale_fill_manual(breaks=c("B.1.1.7","B.1.1.519","B.1.351","B.1.351.1","B.1.427","B.1.429","B.1.526","B.1.526.2","B.1.617.1","B.1.617.2","C.37","P.1","P.1.1","P.2","P.3","Otros"),labels=c("B.1.1.7","B.1.1.519","B.1.351","B.1.351.1","B.1.427","B.1.429","B.1.526","B.1.526.2","B.1.617.1","B.1.617.2","C.37","P.1","P.1.1","P.2","P.3","Otros") )
GraphVariant<- ggpar(GraphVariant,font.label = list(size=14, color = "grey27"),xlab="",ylab="Abundancia relativa")+
  theme_bw() + scale_y_continuous(expand = c(0, 0))+
  scale_x_date(date_breaks = "1 month",date_labels="%b/%Y",date_minor_breaks="1 month", limits = as.Date(c(FechaInicio,Fecha)),expand = c(0,0))+
  theme(panel.border = element_blank(),axis.title=element_text(size=14),axis.text= element_text(size=12,margin=margin(c(0,0,0,-5))))+
  theme( plot.title = element_text(hjust = 0.5, size=8,face = "bold", margin=margin(c(0,0,4,0)))) +
  ggtitle("Variantes en la Republica Mexicana")+
  theme(panel.grid.major = element_blank()) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.box.margin=margin(-10,-5,-5,-10),
    legend.key.size = unit(0.4, "cm"))+ guides(fill=guide_legend(ncol=2))+
  labs(x = paste("n =",nrow(Variant),sep=" "))+
  theme(plot.margin = unit(c(0,0,0,0), "lines"))
#tiff("Test/12Jul2021_VOC_VOI_Total.tiff", width=17, height=11, units="cm", res=300,compression ="lzw")
pdf(paste0(prefix, "_VOC_VOI_Total.jpg"), width=32, height=18, units="cm", res = 300)
GraphVariant
dev.off()


# # # #Graficas por estado
# uniq(linCMX[linCMX$,])
# Location<-unique(linCMX$Location)
# Fecha<-min(linCMX[,2])
# for (i in 1:32){
# 	#Si hay pocas muestras en un estado en la fecha de corte, se debe correr solo para ese estado.
# 	#No se corre el ciclo sin solo se pone i<-3 o 4 o el número de estado que se quiera corregir.
# 	#i<-10
# 	TL<-linCMX[linCMX$Location==Location[i],]
# 	Ttm<-ddply(TL, .(M), summarize,  sum(F))
# 	print (Location[i])
# 	print (Ttm)
# 	n<-sum(Ttm$..1)
# 	L<-ggplot(TL, aes(Collection.date, after_stat(count),fill = Lineage)) +
# 		geom_density(position = "fill",adjust=1.5) +
# 	  scale_fill_manual(breaks=nameLinage,
# 	                    values=c(
# 	                      "B.1.1.7"="darkorange",
# 	                      
# 	                      "B.1.1.519"="indianred1",
# 	                      
# 	                      "B.1.427"="darkolivegreen2",
# 	                      "B.1.429"="darkolivegreen3",
# 	                      
# 	                      "B.1.619"="aquamarine1",
# 	                      "B.1.628"="azure2",
# 	                      "B.1.631"="azure4",
# 	                      
# 	                      "B.1.351"="indianred3",
# 	                      
# 	                      "B.1.617.2"="gold3",
# 	                      "AY.2"="gold",
# 	                      "AY.3"="gold2",
# 	                      "AY.4"="goldenrod2",
# 	                      "AY.5"="goldenrod3",
# 	                      
# 	                      "AY.9"="khaki3",
# 	                      
# 	                      "AY.12"="lightgoldenrod4",
# 	                      "AY.13"="lemonchiffon",
# 	                      "AY.14"="lemonchiffon",
# 	                      "AY.15"="lemonchiffon",
# 	                      
# 	                      "AY.20"="sienna3",
# 	                      
# 	                      "AY.23"="lemonchiffon3",
# 	                      "AY.24"="lemonchiffon3",
# 	                      "AY.25"="lemonchiffon3",
# 	                      "AY.26"="indianred4",
# 	                      "AY.27"="indianred3",
# 	                      "AY.32"="indianred3",
# 	                      
# 	                      "P.1"="green3",
# 	                      "P.1.1"="green4",
# 	                      "P.1.2"="lightseagreen",
# 	                      "P.1.4"="lightseagreen",
# 	                      "P.1.7"="lightseagreen",
# 	                      "P.1.8"="lightseagreen",
# 	                      "P.1.9"="lightseagreen",
# 	                      "P.1.10"="mediumseagreen",
# 	                      "P.1.10.2"="mediumseagreen",
# 	                      
# 	                      "C.37"="blueviolet",
# 	                      "C.37.1"="blueviolet",
# 	                      
# 	                      "B.1.621"="azure4",
# 	                      "B.1.621.1"="azure4",
# 	                      
# 	                      ##################
# 	                      
# 	                      "B.1.526"="steelblue1",
# 	                      
# 	                      "B.1.617.1"="cyan2",
# 	                      
# 	                      ##################
# 	                      
# 	                      "Otros"="magenta"),
# 	                    labels=c(
# 	                      "Alfa (B.1.1.7)",
# 	                      
# 	                      "B.1.1.519**",
# 	                      
# 	                      "B.1.427*",
# 	                      "B.1.429*",
# 	                      
# 	                      "B.1.619",
# 	                      "B.1.628",
# 	                      "B.1.631",
# 	                      
# 	                      "Beta (B.1.351)",
# 	                      
# 	                      "Delta (B.1.617.2)",
# 	                      "Delta (AY.2)",
# 	                      "Delta (AY.3)",
# 	                      "Delta (AY.4)",
# 	                      "Delta (AY.5)",
# 	                      
# 	                      "Delta (AY.9)",
# 	                      
# 	                      "Delta (AY.12)", 
# 	                      "Delta (AY.13)", 
# 	                      "Delta (AY.14)", 
# 	                      "Delta (AY.15)",
# 	                      
# 	                      "Delta (AY.20)",
# 	                      
# 	                      "Delta (AY.23)",
# 	                      "Delta (AY.24)",
# 	                      "Delta (AY.25)",
# 	                      "Delta (AY.26)",
# 	                      "Delta (AY.27)",
# 	                      "Delta (AY.32)",
# 	                      #        
# 	                      "Gamma (P.1)",
# 	                      "Gamma (P.1.1)", 
# 	                      "Gamma (P.1.2)", 
# 	                      "Gamma (P.1.4)", 
# 	                      "Gamma (P.1.7)", 
# 	                      "Gamma (P.1.8)", 
# 	                      "Gamma (P.1.9)",
# 	                      "Gamma (P.1.10)",
# 	                      "Gamma (P.1.10.2)",
# 	                      
# 	                      #        
# 	                      
# 	                      "Lambda (C.37)",
# 	                      "Lambda (C.37.1)",
# 	                      #        
# 	                      "Mu (B.1.621)",
# 	                      "Mu (B.1.621.1)",
# 	                      #        
# 	                      "B.1.526 ***",
# 	                      "B.1.617.1 ***",
# 	                      
# 	                      "Otros"),drop = FALSE ) 
# 	#scale_fill_manual(breaks=c("B.1.1.7","B.1.1.519","B.1.351","B.1.351.1","B.1.427","B.1.429","B.1.526","B.1.526.2","B.1.617.1","B.1.617.2","C.37","P.1","P.1.1","P.2","P.3","Otros"),values=c("B.1.1.7"="orange1","B.1.1.519"="indianred1","B.1.351"="gold3","B.1.351.1"="gold2","B.1.427"="darkolivegreen2","B.1.429"="darkolivegreen3","B.1.526"="green3","B.1.526.2"="mediumseagreen","B.1.617.1"="mediumturquoise","B.1.617.2"="steelblue1","C.37"="dodgerblue","P.1"="slateblue1","P.1.1"="mediumpurple1","P.2"="orchid1","Otros"="hotpink1"),labels=c("B.1.1.7","B.1.1.519","B.1.351","B.1.351.1","B.1.427","B.1.429","B.1.526","B.1.526.2","B.1.617.1","B.1.617.2","C.37","P.1","P.1.1","P.2","P.3","Otros"),drop = FALSE )
# 	L<- ggpar(L,font.label = list(size=9, color = "grey27"),xlab="",ylab="Abundancia relativa")+
# 		theme_bw() + scale_y_continuous(expand = c(0, 0))+
# 		scale_x_date(date_breaks = "1 month",date_labels="%d/%b/%y",date_minor_breaks="1 month", limits = as.Date(c(NA,Fecha)),expand = c(0,0))+
# 		theme(panel.border = element_blank(),axis.text= element_text(margin=margin(c(-10,0,0,-10))))+
# 		theme( plot.title = element_text(hjust = 0.5, size=11.5,face = "bold", margin=margin(c(0,0,10,0)))) +
# 		ggtitle(Location[i])+
# 		theme(panel.grid.major = element_blank()) +
# 		theme(
# 			legend.title = element_text(size = 9),
# 			legend.text = element_text(size = 8),
# 			legend.box.margin=margin(-10,-5,-5,-10),
# 			legend.key.size = unit(0.4, "cm"))+
# 		labs(x = paste("n =",nrow(TL),sep=" "))
# 	tiff(paste(paste("Test/Estados/",paste(Location[i], Fecha, sep="_"), sep=""),".tiff",sep=""), width=27.5, height=22, units="cm", res=300,compression ="lzw")
# 		print (L)
# 	dev.off()
# }
