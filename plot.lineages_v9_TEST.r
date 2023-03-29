library("ggplot2")
library("ggpubr")
library(dplyr)
library(plyr)
setwd("C:/Users/beyon/Dropbox/Difusion-CoViGEn/")

#Tiene dos entradas.
#Archivo separado por tabulador con cuatro columnas: 1) Accession ID, 2)Collection date, 3) Lineage y 4) Location
#linCMX= read.table(file= 'gisaid_hcov-19_1Nov20-18Jun21_Metadata.txt', header=T, sep ='\t' )
#linCMX= read.table(file= 'gisaid_hcov-19_1Dic21-18Jun21_Metadata.txt', header=T, sep ='\t' )
linCMX= read.table(file= 'Metadata/2021-10-28_All_batches_for_areaplots.txt', header=T, sep ='\t', stringsAsFactors = FALSE )

#Archivo de dos columnas separadas por tabulador: 1)Estado y 2) Region
EstReg= read.table(file= 'Metadata/EstadoRegion_v3.txt', header=T, sep ='\t', stringsAsFactors = FALSE )

#Poner la fecha en formato y otra variable que sea por mes
linCMX$Collection.date<- as.Date(linCMX$Collection.date, "%d/%m/%Y")
linCMX$M<-as.factor(format(as.Date(linCMX$Collection.date), "%Y-%m"))

#Cambiar Mexico City a español
Cambiar<-(linCMX$Location %in% c("Mexico City","Distrito Federal"))
linCMX$Location[Cambiar]<-"Ciudad de Mexico"
#Cambiar State of Mexico a español
Cambiar<-(linCMX$Location %in% c("State of Mexico"))
linCMX$Location[Cambiar]<-"Estado de Mexico"
#Poner a cada registro a que región pertenece
linCMX$Region<-linCMX$Location
for (i in 1:32){
  if(EstReg[i,1] %in% linCMX$Location){
    linCMX[linCMX$Location==EstReg[i,1],]$Region<-EstReg[i,2] 
  }
}
FechaInicio<-'2021-02-02'
Fecha<-'2021-09-27'
linCMX<-linCMX[(linCMX$Collection.date>FechaInicio & linCMX$Collection.date<Fecha),]

unique(sort(linCMX$Location))

############
# test <- linCMX
# test <- na.omit(test)
sort(unique(linCMX$Lineage))

# nameLinage <- c("A.1", "A.3", "A.5", 
#                 "B.1", "B.1.1", "B.1.1.133", 
#                 "B.1.1.222", "B.1.1.244", "B.1.1.322", 
#                 "B.1.1.329", "B.1.1.344", "B.1.1.432", 
#                 "B.1.1.512", "B.1.1.517", "B.1.1.519", 
#                 "B.1.1.526", "B.1.111", "B.1.153", "B.1.160", 
#                 "B.1.189", "B.1.2", "B.1.221", "B.1.232", 
#                 "B.1.239", "B.1.241", "B.1.243", "B.1.245", 
#                 "B.1.289", "B.1.320", "B.1.36.39", 
#                 "B.1.366", "B.1.396", "B.1.397", "B.1.399", 
#                 "B.1.400", "B.1.415", "B.1.499", "B.1.503",
#                 "B.1.551", "B.1.558", "B.1.564", "B.1.578", "B.1.609")

nameLinage <- c("A.2.5",
"AZ.3",
"B.1",
"B.1.1",
"B.1.1.121",
"B.1.1.159",
"B.1.1.189",
"B.1.1.207",
"B.1.1.220",
"B.1.1.221",

"B.1.1.222",

"B.1.1.28",
"B.1.1.285",
"B.1.1.291",

"B.1.1.307",

"B.1.1.316",
"B.1.1.318",
"B.1.1.322",
"B.1.1.329",
"B.1.1.362",

"B.1.1.432",

"B.1.1.517",
"B.1.1.519",

"B.1.119",
"B.1.126",
"B.1.127",

"B.1.189",
"B.1.2",
"B.1.232",
"B.1.237",
"B.1.239",
"B.1.241",
"B.1.243",
"B.1.243.1",
"B.1.243.2",

"B.1.396",
"B.1.397",
"B.1.400",
"B.1.404",
"B.1.415.1",
"B.1.427",
"B.1.429",
"B.1.438.1",
"B.1.459",
"B.1.499",
"B.1.526",
"B.1.551",
"B.1.558",
"B.1.561",
"B.1.596",
"B.1.599",
"B.1.609",
"B.1.617.1",
"B.1.623",
"B.1.627",
"B.1.628",
"B.1.631",
"B.1.632",
"B.1.634",
"B.1.635",
"B.1.636",
"B.1.637",
"C.23",
"D.2",
"None",
"P.2",

"B.1.1.7",
"Q.3",

"B.1.351",

"B.1.617.2",
"AY.13",
"AY.14",
"AY.15",
"AY.17",
"AY.2",
"AY.20",
"AY.21",
"AY.25",
"AY.26",
"AY.3",
"AY.37",
"AY.4",

"P.1",
"P.1.10",
"P.1.10.2",
"P.1.12",
"P.1.2",
"P.1.3",
"P.1.7",
"P.1.9",

"C.37",
"C.37.1",

"B.1.621",
"B.1.621.1")

unique(sort(linCMX$Lineage))
unique(linCMX$Lineage)

# #Graficas por regiones
Titulos<-c("Region NW (BCN, BCS, DUR, CHH, SON, SIN)","Region NE (COA, NLE,TAM)","Region WE (COL, JAL,MIC, NAY)","Region CN (AGU, GUA, QUE, SLP, ZAC)","Region CS (CMX, HID, MEX,MOR, PUE, TLA)","Region SE (CHP, GRO, OAX, TAB, VER)","Region S (CAM, ROO, YUC)")
Titulos<-c("Region NW (BCN, BCS, DUR, CHH, SON, SIN)","Region NE (COA, NLE,TAM)","Region WE (COL, JAL,MIC, NAY)","Region CN (AGU, GUA, QUE, SLP, ZAC)","Region CS (CMX, HID, MEX,MOR, PUE, TLA)", "Region SE (CAM, ROO, YUC)", "Region S (CHP, GRO, OAX, TAB, VER)")
Regiones<-c("NW","NE","WE","CN","CS","SE","S")
#Fecha final por region
FechaInicio<-c('2021-02-03','2021-02-03','2021-02-03','2021-02-03','2021-02-03','2021-02-03','2021-02-03')
Fecha<-c('2021-09-29','2021-09-29','2021-09-29','2021-09-29','2021-09-29','2021-09-29','2021-09-29')
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
    geom_density(position = "fill",adjust=2.0) +  #adjust se usa para hacer mas grande el promedio en los valores relativo, quita ruido y picos extraños
# scale_fill_manual(breaks=
#                     nameLinage,
#                   values=c(
#                            "A.1"	=	 "white",
#                            "A.3"	=	 "white",
#                            "A.5"	=	 "white",
#                            
#                            "B.1"	=	 "pink",
#                            "B.1.1"	=	 "gray21",
#                            "B.1.1.133"	=	 "gray27",
#                            "B.1.1.222"	=	 "firebrick3",
#                            "B.1.1.244"	=	 "gray61",
#                            "B.1.1.322"	=	 "pink2",
#                            "B.1.1.329"	=	 "pink2",
#                            "B.1.1.344"	=	 "lightsteelblue3",
#                            "B.1.1.432"	=	 "skyblue",
#                            
#                            "B.1.1.512"	=	 "aliceblue",
#                            "B.1.1.517"	=	 "aliceblue",
#                            "B.1.1.519"	=	 "indianred1",
#                            "B.1.1.526"	=	 "aliceblue",
#                            
#                            "B.1.111"	=	 "steelblue",
#                            "B.1.153"	=	 "steelblue",
#                            "B.1.160"	=	 "steelblue",
#                            "B.1.189"	=	 "steelblue",
#                            "B.1.2"	=	 "pink3",
#                            "B.1.221"	=	 "pink3",
#                            "B.1.232"	=	 "pink3",
#                            "B.1.239"	=	 "pink3",
#                            "B.1.241"	=	 "pink3",
#                            "B.1.243"	=	 "pink3",
#                            "B.1.245"	=	 "pink4",
#                            "B.1.289"	=	 "pink4",
#                            "B.1.320"	=	 "pink4",
#                            "B.1.36.39"	=	 "darkseagreen",
#                            "B.1.366"	=	 "darkseagreen",
#                            "B.1.396"	=	 "darkseagreen",
#                            "B.1.397"	=	 "darkseagreen",
#                            "B.1.399"	=	 "darkseagreen",
#                            "B.1.400"	=	 "darkseagreen1",
#                            "B.1.415"	=	 "darkseagreen1",
#                            "B.1.499"	=	 "darkseagreen2",
#                            "B.1.503"	=	 "darkseagreen2",
#                            "B.1.551"	=	 "darkslateblue",
#                            "B.1.558"	=	 "darkslateblue",
#                            "B.1.564"	=	 "darkslateblue",
#                            "B.1.578"	=	 "lavenderblush",
#                            "B.1.609"	=	 "darkslategray"
#                            
#                   ),
#                   labels=nameLinage
#                   ,drop = FALSE )
  scale_fill_manual(breaks=
                      nameLinage,
                    values=c(
                             "A.2.5"	=	 "white",
                             "AZ.3"	=	 "white",
                             "B.1"	=	 "pink",
                             "B.1.1"	=	 "gray21",
                             "B.1.1.10" =	 "gray22",
                             "B.1.1.121"	=	 "gray23",
                             
                             "B.1.1.159"	=	 "gray27",
                             
                             "B.1.1.189"	=	 "gray28",
                             
                             "B.1.1.207"	=	 "gray61",
                             "B.1.1.221"	=	 "gray61",
                             "B.1.1.222"	=	 "firebrick3",
                             "B.1.1.28"	=	 "gray61",
                             "B.1.1.285"	=	 "gray61",
                             "B.1.1.291"	=	 "gray61",
                             
                             "B.1.1.307"	=	 "pink2",
                             "B.1.1.316"	=	 "pink2",
                             "B.1.1.318"	=	 "pink2",
                             "B.1.1.322"	=	 "pink2",
                             "B.1.1.329"	=	 "pink2",
                             "B.1.1.362"	=	 "lightsteelblue3",
                             "B.1.1.432"	=	 "skyblue",
                             
                             "B.1.1.517"	=	 "aliceblue",
                             "B.1.1.519"	=	 "indianred1",
                             
                             "B.1.119"	=	 "steelblue",
                             "B.1.126"	=	 "steelblue",
                             "B.1.127"	=	 "steelblue",
                             "B.1.189"	=	 "steelblue",
                             "B.1.2"	=	 "pink3",
                             "B.1.232"	=	 "pink3",
                             "B.1.237"	=	 "pink3",
                             "B.1.239"	=	 "pink3",
                             "B.1.241"	=	 "pink3",
                             "B.1.243"	=	 "pink3",
                             "B.1.243.1"	=	 "pink3",
                             "B.1.243.2"	=	 "pink3",
                             
                             "B.1.396"	=	 "darkseagreen",
                             "B.1.397"	=	 "darkseagreen",
                             "B.1.400"	=	 "darkseagreen1",
                             "B.1.404"	=	 "darkseagreen1",
                             "B.1.415.1"	=	 "darkseagreen1",
                             "B.1.427"	=	 "darkolivegreen2",
                             "B.1.429"	=	 "darkolivegreen3",
                             "B.1.438.1"	=	 "darkseagreen2",
                             "B.1.459"	=	 "darkseagreen2",
                             "B.1.499"	=	 "darkseagreen2",
                             
                             "B.1.526"	=	 "steelblue1",
                             
                             "B.1.551"	=	 "darkslateblue",
                             "B.1.558"	=	 "darkslateblue",
                             "B.1.561"	=	 "darkslateblue",
                             
                             "B.1.596"	=	 "lavenderblush",
                             "B.1.599"	=	 "lavenderblush",
                             "B.1.609"	=	 "darkslategray",
                             
                             "B.1.617.1"	=	 "cyan2",
                             
                             "B.1.623"	=	 "darkslategray",
                             "B.1.627"	=	 "darkslategray",
                             "B.1.628"	=	 "azure2",
                             "B.1.631"	=	 "darkslategray",
                             "B.1.632"	=	 "darkslategray",
                             "B.1.634"	=	 "darkslategray",
                             "B.1.635"	=	 "darkslategray",
                             "B.1.636"	=	 "darkslategray",
                             "B.1.637"	=	 "darkslategray",
                             "C.23"	=	 "darkslategray3",
                             "D.2"	=	 "darkslategray3",
                             "None"	=	 "black",
                             "P.2"	=	 "darkslategray4",
                             
                             "B.1.1.7"	=	 "darkorange",
                             "Q.3"	=	 "darkorange1",
                             
                             "B.1.351"	=	 "lightyellow3",
                             
                             "B.1.617.2"	=	 "gold3",
                             
                             "AY.13"	=	 "lemonchiffon",
                             "AY.14"	=	 "lemonchiffon",
                             "AY.15"	=	 "lemonchiffon",
                             "AY.17"	=	 "lemonchiffon",
                             
                             "AY.2"	=	 "gold",
                             
                             "AY.20"	=	 "sienna3",
                             "AY.21"	=	 "lemonchiffon3",
                             "AY.25"	=	 "lemonchiffon3",
                             "AY.26"	=	 "indianred4",
                             
                             "AY.3"	=	 "gold2",
                             "AY.37"	=	 "lightyellow3",
                             "AY.4"	=	 "goldenrod2",
                             
                             "P.1"	=	 "green3",
                             "P.1.10"	=	 "mediumseagreen",
                             "P.1.10.2"	=	 "mediumseagreen",
                             "P.1.12"	=	 "mediumseagreen",
                             "P.1.2"	=	 "lightseagreen",
                             "P.1.3"	=	 "lightseagreen",
                             "P.1.7"	=	 "lightseagreen",
                             "P.1.9"	=	 "lightseagreen",
                             
                             "C.37"	=	 "blueviolet",
                             "C.37.1"	=	 "blueviolet",
                             
                             "B.1.621"	=	 "azure4",
                             "B.1.621.1"	=	 "azure4"
                             
                    ),
                    labels=c(
                             "A.2.5",
                             "AZ.3",
                             "B.1",
                             "B.1.1",
                             "B.1.1.121",
                             "B.1.1.159",
                             "B.1.1.189",
                             "B.1.1.207",
                             "B.1.1.220",
                             "B.1.1.221",
                             
                             "B.1.1.222",
                             
                             "B.1.1.28",
                             "B.1.1.285",
                             "B.1.1.291",
                             
                             "B.1.1.307",
                             
                             "B.1.1.316",
                             "B.1.1.318",
                             "B.1.1.322",
                             "B.1.1.329",
                             "B.1.1.362",
                             
                             "B.1.1.432",
                             
                             "B.1.1.517",
                             "B.1.1.519",
                             
                             "B.1.119",
                             "B.1.126",
                             "B.1.127",
                             
                             "B.1.189",
                             "B.1.2",
                             "B.1.232",
                             "B.1.237",
                             "B.1.239",
                             "B.1.241",
                             "B.1.243",
                             "B.1.243.1",
                             "B.1.243.2",
                             
                             "B.1.396",
                             "B.1.397",
                             "B.1.400",
                             "B.1.404",
                             "B.1.415.1",
                             "B.1.427",
                             "B.1.429",
                             "B.1.438.1",
                             "B.1.459",
                             "B.1.499",
                             "B.1.526",
                             "B.1.551",
                             "B.1.558",
                             "B.1.561",
                             "B.1.596",
                             "B.1.599",
                             "B.1.609",
                             "B.1.617.1",
                             "B.1.623",
                             "B.1.627",
                             "B.1.628",
                             "B.1.631",
                             "B.1.632",
                             "B.1.634",
                             "B.1.635",
                             "B.1.636",
                             "B.1.637",
                             "C.23",
                             "D.2",
                             "None",
                             "P.2",
                             
                             "Alfa (B.1.1.7)",
                             "Alfa (Q.3)",
                             
                             "Beta (B.1.351)",
                             
                             "Delta (B.1.617.2)",
                             "Delta (AY.13)",
                             "Delta (AY.14)",
                             "Delta (AY.15)",
                             "Delta (AY.17)",
                             "Delta (AY.2)",
                             "Delta (AY.20)",
                             "Delta (AY.21)",
                             "Delta (AY.25)",
                             "Delta (AY.26)",
                             "Delta (AY.3)",
                             "Delta (AY.37)",
                             "Delta (AY.4)",
                             
                             "Gamma (P.1)",
                             "Gamma (P.1.10)",
                             "Gamma (P.1.10.2)",
                             "Gamma (P.1.12)",
                             "Gamma (P.1.2)",
                             "Gamma (P.1.3)",
                             "Gamma (P.1.7)",
                             "Gamma (P.1.9)",
                             
                             "Lambda (C.37)",
                             "Lambda (C.37.1)",
                             
                             "Mu (B.1.621)",
                             "Mu (B.1.621.1)"
                    )
                    ,drop = FALSE )
  #scale_fill_manual(breaks=c("B.1.1.7","B.1.1.519","B.1.351","B.1.351.1","B.1.427","B.1.429","B.1.526","B.1.526.2","B.1.617.1","B.1.617.2","C.37","P.1","P.1.1","P.2","P.3","Otros"),labels=c("B.1.1.7","B.1.1.519","B.1.351","B.1.351.1","B.1.427","B.1.429","B.1.526","B.1.526.2","B.1.617.1","B.1.617.2","C.37","P.1","P.1.1","P.2","P.3","Otros") )
  GraphVariant<- ggpar(GraphVariant,font.label = list(size=8, color = "grey27"),xlab="",ylab="Abundancia relativa")+
    theme_bw() + scale_y_continuous(expand = c(0, 0))+
    scale_x_date(date_breaks = "2 month",date_labels="%b/%y",date_minor_breaks="2 month", limits = as.Date(c(FechaInicio[i],Fecha[i])),expand = c(0,0))+
    theme(panel.border = element_blank(),axis.title=element_text(size=8),axis.text= element_text(size=8,margin=margin(c(-1,0,0,-5))))+
    theme( plot.title = element_text(hjust = 0.5, size=8,face = "bold", margin=margin(c(5,0,1,0)))) +
    ggtitle(Titulos[i])+
    theme(panel.grid.major = element_blank()) +
    theme(
      legend.title = element_blank(),
      legend.position="right",
      legend.margin=margin(),
      legend.box="vertical",
      legend.text = element_text(size = 6),
      #legend.box.margin=margin(0,10,-5,10),
      legend.key.size = unit(0.2, "cm"))+ guides(fill=guide_legend(ncol=3))+
    labs(x = paste("n =",nrow(Variant),sep=" ")) +
    theme(plot.margin = unit(c(0,0,0,0.4), "lines"))
  Imag[[i]]<-GraphVariant
  #tiff(paste(paste("Test/Regiones/", "test2_VOC_VOI_", sep=""),Regiones[i],".tiff",sep=""), width=60, height=35, units="cm", res=300,compression ="lzw")
  jpeg(paste(paste("Test/Regiones/", "test2_VOC_VOI_", sep=""),Regiones[i],".jpeg",sep=""), width=20, height=24, units="cm", res=300)
  print (GraphVariant)
  dev.off()
}
leg<-get_legend(Imag[[7]] + theme(legend.text = element_text(size = 6)) + guides(fill=guide_legend(ncol=6)))
for (i in 1:length(Regiones)){
  Imag[[i]]<-Imag[[i]] + theme(legend.position = "none")
}
#tiff("Test/Regiones/test2_VOC_VOI_Regiones.tiff", width=60, height=80, units="cm", res=300,compression ="lzw")
jpeg("Test/Regiones/test2_VOC_VOI_Regiones.jpeg", width=20, height=24, units="cm", res=300)
ggarrange(Imag[[1]],Imag[[2]],Imag[[3]],Imag[[4]],Imag[[5]],Imag[[6]],Imag[[7]], leg, ncol=2,nrow=4,common.legend = FALSE,widths=c(2,2),labels = c("A","B","C","D","E","F","G"),font.label = list(size = 9),hjust=-1.2)
dev.off()
#
# # # ##Grafica a nivel nacional
FechaInicio<-'2021-02-03'
Fecha<-'2021-09-29'
linCMX$F<-1
#Obtener todos los registros
Variant<-linCMX
Variant$Fac<-factor(Variant$Lineage,levels=nameLinage)
#Sacar cuantos registros hay por mes a nivel nacional
Ttm<-ddply(Variant, .(M), summarize,  sum(F))
print (Ttm)
GraphVariant<-ggplot(Variant, aes(Collection.date, after_stat(count),fill = Fac)) +
  geom_density(position = "fill",adjust=2.0) +
  # scale_fill_manual(breaks=
  #                     nameLinage,
  #                   values=c(
  #                            "A.1"	=	 "white",
  #                            "A.3"	=	 "white",
  #                            "A.5"	=	 "white",
  #                            
  #                            "B.1"	=	 "pink",
  #                            "B.1.1"	=	 "gray21",
  #                            "B.1.1.133"	=	 "gray27",
  #                            "B.1.1.222"	=	 "firebrick3",
#                            "B.1.1.244"	=	 "gray61",
#                            "B.1.1.322"	=	 "pink2",
#                            "B.1.1.329"	=	 "pink2",
#                            "B.1.1.344"	=	 "lightsteelblue3",
#                            "B.1.1.432"	=	 "skyblue",
#                            
#                            "B.1.1.512"	=	 "aliceblue",
#                            "B.1.1.517"	=	 "aliceblue",
#                            "B.1.1.519"	=	 "indianred1",
#                            "B.1.1.526"	=	 "aliceblue",
#                            
#                            "B.1.111"	=	 "steelblue",
#                            "B.1.153"	=	 "steelblue",
#                            "B.1.160"	=	 "steelblue",
#                            "B.1.189"	=	 "steelblue",
#                            "B.1.2"	=	 "pink3",
#                            "B.1.221"	=	 "pink3",
#                            "B.1.232"	=	 "pink3",
#                            "B.1.239"	=	 "pink3",
#                            "B.1.241"	=	 "pink3",
#                            "B.1.243"	=	 "pink3",
#                            "B.1.245"	=	 "pink4",
#                            "B.1.289"	=	 "pink4",
#                            "B.1.320"	=	 "pink4",
#                            "B.1.36.39"	=	 "darkseagreen",
#                            "B.1.366"	=	 "darkseagreen",
#                            "B.1.396"	=	 "darkseagreen",
#                            "B.1.397"	=	 "darkseagreen",
#                            "B.1.399"	=	 "darkseagreen",
#                            "B.1.400"	=	 "darkseagreen1",
#                            "B.1.415"	=	 "darkseagreen1",
#                            "B.1.499"	=	 "darkseagreen2",
#                            "B.1.503"	=	 "darkseagreen2",
#                            "B.1.551"	=	 "darkslateblue",
#                            "B.1.558"	=	 "darkslateblue",
#                            "B.1.564"	=	 "darkslateblue",
#                            "B.1.578"	=	 "lavenderblush",
#                            "B.1.609"	=	 "darkslategray"
#                            
#                   ),
#                   labels=nameLinage
#                   ,drop = FALSE )
scale_fill_manual(breaks=
                    nameLinage,
                  values=c(
                    "A.2.5"	=	 "white",
                    "AZ.3"	=	 "white",
                    "B.1"	=	 "pink",
                    "B.1.1"	=	 "gray21",
                    "B.1.1.10" =	 "gray22",
                    "B.1.1.121"	=	 "gray23",
                    
                    "B.1.1.159"	=	 "gray27",
                    
                    "B.1.1.189"	=	 "gray28",
                    
                    "B.1.1.207"	=	 "gray61",
                    "B.1.1.221"	=	 "gray61",
                    "B.1.1.222"	=	 "firebrick3",
                    "B.1.1.28"	=	 "gray61",
                    "B.1.1.285"	=	 "gray61",
                    "B.1.1.291"	=	 "gray61",
                    
                    "B.1.1.307"	=	 "pink2",
                    "B.1.1.316"	=	 "pink2",
                    "B.1.1.318"	=	 "pink2",
                    "B.1.1.322"	=	 "pink2",
                    "B.1.1.329"	=	 "pink2",
                    "B.1.1.362"	=	 "lightsteelblue3",
                    "B.1.1.432"	=	 "skyblue",
                    
                    "B.1.1.517"	=	 "aliceblue",
                    "B.1.1.519"	=	 "indianred1",
                    
                    "B.1.119"	=	 "steelblue",
                    "B.1.126"	=	 "steelblue",
                    "B.1.127"	=	 "steelblue",
                    "B.1.189"	=	 "steelblue",
                    "B.1.2"	=	 "pink3",
                    "B.1.232"	=	 "pink3",
                    "B.1.237"	=	 "pink3",
                    "B.1.239"	=	 "pink3",
                    "B.1.241"	=	 "pink3",
                    "B.1.243"	=	 "pink3",
                    "B.1.243.1"	=	 "pink3",
                    "B.1.243.2"	=	 "pink3",
                    
                    "B.1.396"	=	 "darkseagreen",
                    "B.1.397"	=	 "darkseagreen",
                    "B.1.400"	=	 "darkseagreen1",
                    "B.1.404"	=	 "darkseagreen1",
                    "B.1.415.1"	=	 "darkseagreen1",
                    "B.1.427"	=	 "darkolivegreen2",
                    "B.1.429"	=	 "darkolivegreen3",
                    "B.1.438.1"	=	 "darkseagreen2",
                    "B.1.459"	=	 "darkseagreen2",
                    "B.1.499"	=	 "darkseagreen2",
                    
                    "B.1.526"	=	 "steelblue1",
                    
                    "B.1.551"	=	 "darkslateblue",
                    "B.1.558"	=	 "darkslateblue",
                    "B.1.561"	=	 "darkslateblue",
                    
                    "B.1.596"	=	 "lavenderblush",
                    "B.1.599"	=	 "lavenderblush",
                    "B.1.609"	=	 "darkslategray",
                    
                    "B.1.617.1"	=	 "cyan2",
                    
                    "B.1.623"	=	 "darkslategray",
                    "B.1.627"	=	 "darkslategray",
                    "B.1.628"	=	 "azure2",
                    "B.1.631"	=	 "darkslategray",
                    "B.1.632"	=	 "darkslategray",
                    "B.1.634"	=	 "darkslategray",
                    "B.1.635"	=	 "darkslategray",
                    "B.1.636"	=	 "darkslategray",
                    "B.1.637"	=	 "darkslategray",
                    "C.23"	=	 "darkslategray3",
                    "D.2"	=	 "darkslategray3",
                    "None"	=	 "black",
                    "P.2"	=	 "darkslategray4",
                    
                    "B.1.1.7"	=	 "darkorange",
                    "Q.3"	=	 "darkorange1",
                    
                    "B.1.351"	=	 "lightyellow3",
                    
                    "B.1.617.2"	=	 "gold3",
                    
                    "AY.13"	=	 "lemonchiffon",
                    "AY.14"	=	 "lemonchiffon",
                    "AY.15"	=	 "lemonchiffon",
                    "AY.17"	=	 "lemonchiffon",
                    
                    "AY.2"	=	 "gold",
                    
                    "AY.20"	=	 "sienna3",
                    "AY.21"	=	 "lemonchiffon3",
                    "AY.25"	=	 "lemonchiffon3",
                    "AY.26"	=	 "indianred4",
                    
                    "AY.3"	=	 "gold2",
                    "AY.37"	=	 "lightyellow3",
                    "AY.4"	=	 "goldenrod2",
                    
                    "P.1"	=	 "green3",
                    "P.1.10"	=	 "mediumseagreen",
                    "P.1.10.2"	=	 "mediumseagreen",
                    "P.1.12"	=	 "mediumseagreen",
                    "P.1.2"	=	 "lightseagreen",
                    "P.1.3"	=	 "lightseagreen",
                    "P.1.7"	=	 "lightseagreen",
                    "P.1.9"	=	 "lightseagreen",
                    
                    "C.37"	=	 "blueviolet",
                    "C.37.1"	=	 "blueviolet",
                    
                    "B.1.621"	=	 "azure4",
                    "B.1.621.1"	=	 "azure4"
                    
                  ),
                  labels=c(
                    "A.2.5",
                    "AZ.3",
                    "B.1",
                    "B.1.1",
                    "B.1.1.121",
                    "B.1.1.159",
                    "B.1.1.189",
                    "B.1.1.207",
                    "B.1.1.220",
                    "B.1.1.221",
                    
                    "B.1.1.222",
                    
                    "B.1.1.28",
                    "B.1.1.285",
                    "B.1.1.291",
                    
                    "B.1.1.307",
                    
                    "B.1.1.316",
                    "B.1.1.318",
                    "B.1.1.322",
                    "B.1.1.329",
                    "B.1.1.362",
                    
                    "B.1.1.432",
                    
                    "B.1.1.517",
                    "B.1.1.519",
                    
                    "B.1.119",
                    "B.1.126",
                    "B.1.127",
                    
                    "B.1.189",
                    "B.1.2",
                    "B.1.232",
                    "B.1.237",
                    "B.1.239",
                    "B.1.241",
                    "B.1.243",
                    "B.1.243.1",
                    "B.1.243.2",
                    
                    "B.1.396",
                    "B.1.397",
                    "B.1.400",
                    "B.1.404",
                    "B.1.415.1",
                    "B.1.427",
                    "B.1.429",
                    "B.1.438.1",
                    "B.1.459",
                    "B.1.499",
                    "B.1.526",
                    "B.1.551",
                    "B.1.558",
                    "B.1.561",
                    "B.1.596",
                    "B.1.599",
                    "B.1.609",
                    "B.1.617.1",
                    "B.1.623",
                    "B.1.627",
                    "B.1.628",
                    "B.1.631",
                    "B.1.632",
                    "B.1.634",
                    "B.1.635",
                    "B.1.636",
                    "B.1.637",
                    "C.23",
                    "D.2",
                    "None",
                    "P.2",
                    
                    "Alfa (B.1.1.7)",
                    "Alfa (Q.3)",
                    
                    "Beta (B.1.351)",
                    
                    "Delta (B.1.617.2)",
                    "Delta (AY.13)",
                    "Delta (AY.14)",
                    "Delta (AY.15)",
                    "Delta (AY.17)",
                    "Delta (AY.2)",
                    "Delta (AY.20)",
                    "Delta (AY.21)",
                    "Delta (AY.25)",
                    "Delta (AY.26)",
                    "Delta (AY.3)",
                    "Delta (AY.37)",
                    "Delta (AY.4)",
                    
                    "Gamma (P.1)",
                    "Gamma (P.1.10)",
                    "Gamma (P.1.10.2)",
                    "Gamma (P.1.12)",
                    "Gamma (P.1.2)",
                    "Gamma (P.1.3)",
                    "Gamma (P.1.7)",
                    "Gamma (P.1.9)",
                    
                    "Lambda (C.37)",
                    "Lambda (C.37.1)",
                    
                    "Mu (B.1.621)",
                    "Mu (B.1.621.1)"
                  )
                  ,drop = FALSE )
#scale_fill_manual(breaks=c("B.1.1.7","B.1.1.519","B.1.351","B.1.351.1","B.1.427","B.1.429","B.1.526","B.1.526.2","B.1.617.1","B.1.617.2","C.37","P.1","P.1.1","P.2","P.3","Otros"),labels=c("B.1.1.7","B.1.1.519","B.1.351","B.1.351.1","B.1.427","B.1.429","B.1.526","B.1.526.2","B.1.617.1","B.1.617.2","C.37","P.1","P.1.1","P.2","P.3","Otros") )
GraphVariant<- ggpar(GraphVariant,font.label = list(size=8, color = "grey27"),xlab="",ylab="Abundancia relativa")+
  theme_bw() + scale_y_continuous(expand = c(0, 0))+
  scale_x_date(date_breaks = "2 month",date_labels="%b/%y",date_minor_breaks="2 month", limits = as.Date(c(FechaInicio,Fecha)),expand = c(0,0))+
  theme(panel.border = element_blank(),axis.title=element_text(size=8),axis.text= element_text(size=8,margin=margin(c(0,0,0,-5))))+
  theme( plot.title = element_text(hjust = 0.5, size=8,face = "bold", margin=margin(c(0,0,4,0)))) +
  ggtitle("Variantes en la Republica Mexicana")+
  theme(panel.grid.major = element_blank()) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.box.margin=margin(-10,-5,-5,-10),
    legend.key.size = unit(0.4, "cm"))+ guides(fill=guide_legend(ncol=3))+
    
  labs(x = paste("n =",nrow(Variant),sep=" "))+
  theme(plot.margin = unit(c(0,0,0,0), "lines"))
#tiff("Test/Regiones/test2_VOC_VOI_Total.tiff", width=60, height=28, units="cm", res=300,compression ="lzw")
jpeg("Test/Regiones/test2_VOC_VOI_Total.jpg", width=30, height=20, units="cm", res = 300)
GraphVariant
dev.off()