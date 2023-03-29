library("ggplot2")
library("ggpubr")
library(dplyr)
library(plyr)
setwd("C:/Users/beyon/Dropbox/Difusion-CoViGEn/")

#Tiene dos entradas.
#Archivo separado por tabulador con cuatro columnas: 1) Accession ID, 2)Collection date, 3) Lineage y 4) Location
#linCMX= read.table(file= 'gisaid_hcov-19_1Nov20-18Jun21_Metadata.txt', header=T, sep ='\t' )
#linCMX= read.table(file= 'gisaid_hcov-19_1Dic21-18Jun21_Metadata.txt', header=T, sep ='\t' )
linCMX= read.table(file= 'Metadata/gisaid_hcov-19_13Septiembre2021_Metadata.txt', header=T, sep ='\t', stringsAsFactors = FALSE )

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
  linCMX[linCMX$Location==EstReg[i,1],]$Region<-EstReg[i,2] 
}
FechaInicio<-'2020-01-01'
Fecha<-'2021-09-14'
linCMX<-linCMX[(linCMX$Collection.date>FechaInicio & linCMX$Collection.date<Fecha),]

test <- linCMX
test <- na.omit(test) 

test$Lineage[test$Lineage%in% c('B.1.1.7', 'Q.1', 'Q.3')] <- 'Alfa'
test$Lineage[test$Lineage %in% c('B.1.351', 'B.1.351.1')] <- 'Beta'
test$Lineage[test$Lineage %in% c("P.1",
                                 "P.1.1", 
                                 "P.1.2",
                                 "P.1.4",
                                 "P.1.7",
                                 "P.1.8",
                                 "P.1.9")] <- 'Gamma'
test$Lineage[test$Lineage %in% c("B.1.617.2",
                                 "AY.2", 
                                 "AY.3", 
                                 "AY.3.1",
                                 "AY.4", 
                                 "AY.5",
                                 "AY.6",
                                 "AY.7.1",
                                 "AY.7.2",
                                 "AY.9", 
                                 "AY.10", 
                                 "AY.11", 
                                 "AY.12",
                                 "AY.13",
                                 "AY.14",
                                 "AY.15",
                                 "AY.16",
                                 "AY.17",
                                 "AY.18",
                                 "AY.19",
                                 "AY.20",
                                 "AY.21",
                                 "AY.22",
                                 "AY.23",
                                 "AY.24",
                                 "AY.25")] <- 'Delta'
test$Lineage[test$Lineage == 'B.1.526'] <- 'Iota'
test$Lineage[test$Lineage == 'B.1.617.1'] <- 'Kappa'
test$Lineage[test$Lineage %in% c("C.37",
                                 "C.37.1")] <- 'Lambda'
test$Lineage[test$Lineage %in% c("B.1.621",
                                 "B.1.621.1")] <- 'Mu'

b <- count(test, vars = "Lineage")
b$freq <- (b$freq / nrow(test)) * 100
c <- filter(b, freq >= 0.5)

nameLinage <- sort(c$Lineage)

Cambiar<-!(test$Lineage %in% nameLinage)
test$Lineage[Cambiar]<-"Otros"

nameLinage <- sort(unique(test$Lineage))
####################################
FechaInicio<-'2020-01-01'
Fecha<-'2021-08-31'
test$F<-1
#Obtener todos los registros
Variant<-test
Variant$Fac<-factor(Variant$Lineage,levels=nameLinage)
#Sacar cuantos registros hay por mes a nivel nacional
Ttm<-ddply(Variant, .(M), summarize,  sum(F))
print (Ttm)
GraphVariant<-ggplot(Variant, aes(Collection.date, after_stat(count),fill = Fac)) +
  geom_density(position = "fill",adjust=2.5) +
  scale_fill_manual(breaks=
    nameLinage,
    values=c("Alfa" = "darkorange",
             "Beta" = "indianred3",
             "Gamma" = "green3",
             "Delta" = "gold3",
             "Iota" = "steelblue1",
             "Kappa" = "cyan2",
             "Lambda" = "blueviolet",
             "Mu" = "azure4",
             "A" = "ivory1",
             "B.1.1" = "blue",
             "None" = "black",
             "B.1.243" = "bisque1",
             "B.1.609" = "darkkhaki",
             "B.1.561" = "darkslategray4",
             "B.1.427" = "darkolivegreen2",
             "B.1.429" = "darkolivegreen3",
             "B.1.2" = "mediumpurple1",
             "B.1.619" = "aquamarine1",
             "B.1.620" = "aquamarine2",
             "B.1.628" = "cornflowerblue",
             "B.1.1.318" = "lightseagreen",
             "B.1.1.519" = "indianred1",
             "B.1.1.222" = "firebrick3",
             "B.1" = "lightpink",
             
             
             "A.1" = "white",
             "A.2" = "white",
             "A.2.3" = "white",
             "A.2.5" = "white",
             "A.2.5.1" = "white",
             "A.2.5.2" = "white",
             "A.3" = "white",
             "A.5" = "white",
             "AZ.3" = "lightpink4",
             "B" = "tan",
             "B.1.1.10" = "gray21",
             "B.1.1.101" = "gray21",
             "B.1.1.117" = "gray21",
             "B.1.1.128" = "gray21",
             "B.1.1.129" = "gray21",
             "B.1.1.133" = "gray21",
             "B.1.1.161" = "gray21",
             "B.1.1.182" = "gray21",
             "B.1.1.198" = "gray21",
             "B.1.1.205" = "gray61",
             "B.1.1.207" = "gray61",
             "B.1.1.220" = "gray61",
             "B.1.1.231" = "gray61",
             "B.1.1.236" = "gray61",
             "B.1.1.244" = "gray61",
             "B.1.1.274" = "gray61",
             "B.1.1.28" = "gray61",
             "B.1.1.285" = "gray61",
             "B.1.1.307" = "pink2",
             "B.1.1.310" = "pink2",
             "B.1.1.311" = "pink2",
             "B.1.1.316" = "pink2",
             "B.1.1.317" = "pink2",
             "B.1.1.322" = "pink2",
             "B.1.1.329" = "pink2",
             "B.1.1.33" = "lightsteelblue1",
             "B.1.1.334" = "lightsteelblue1",
             "B.1.1.344" = "lightsteelblue1",
             "B.1.1.348" = "lightsteelblue1",
             "B.1.1.362" = "lightsteelblue1",
             "B.1.1.398" = "lightsteelblue1",
             "B.1.1.343" = "skyblue",
             "B.1.1.416" = "skyblue",
             "B.1.1.418" = "skyblue",
             "B.1.1.422" = "skyblue",
             "B.1.1.432" = "skyblue",
             "B.1.1.434" = "skyblue",
             "B.1.1.47" = "skyblue",
             "B.1.1.512" = "aliceblue",
             "B.1.1.517" = "aliceblue",
             "B.1.1.518" = "aliceblue",
             "B.1.1.526" = "aliceblue",
             "B.1.1.70" = "aliceblue",
             "B.1.1.71" = "aliceblue",
             "B.1.1.8" = "aliceblue",
             "B.1.1.93" = "aliceblue",
             "B.1.105" = "steelblue",
             "B.1.111" = "steelblue",
             "B.1.119" = "steelblue",
             "B.1.126" = "steelblue",
             "B.1.127" = "steelblue",
             "B.1.128" = "steelblue",
             "B.1.160.32" = "steelblue",
             "B.1.165" = "steelblue",
             "B.1.177" = "steelblue",
             "B.1.177.53" = "steelblue",
             "B.1.179" = "steelblue",
             "B.1.189" = "steelblue",
             "B.1.201" = "pink3",
             "B.1.206" = "pink3",
             "B.1.208" = "pink3",
             "B.1.214.3" = "pink3",
             "B.1.221" = "pink3",
             "B.1.229" = "pink3",
             "B.1.232" = "pink3",
             "B.1.234" = "pink3",
             "B.1.236" = "pink3",
             "B.1.239" = "pink3",
             "B.1.240" = "pink3",
             "B.1.241" = "pink3",
             "B.1.243" = "pink3",
             "B.1.243.1" = "pink3",
             "B.1.245" = "pink3",
             "B.1.258.21" = "pink3",
             "B.1.267" = "pink3",
             "B.1.319" = "pink4",
             "B.1.320" = "pink4",
             "B.1.324" = "pink4",
             "B.1.333" = "pink4",
             "B.1.346" = "pink4",
             "B.1.349" = "pink4",
             "B.1.36.10" = "pink4",
             "B.1.36.31" = "pink4",
             "B.1.361" = "pink4",
             "B.1.366" = "pink4",
             "B.1.369" = "pink4",
             "B.1.375" = "pink4",
             "B.1.393" = "pink4",
             "B.1.395" = "pink4",
             "B.1.396" = "pink4",
             "B.1.397" = "pink4",
             "B.1.398" = "pink4",
             "B.1.399" = "pink4",
             "B.1.400" = "darkseagreen",
             "B.1.400.1" = "darkseagreen",
             "B.1.402" = "darkseagreen",
             "B.1.404" = "darkseagreen",
             "B.1.405" = "darkseagreen",
             "B.1.409" = "darkseagreen",
             "B.1.411" = "darkseagreen",
             "B.1.415" = "darkseagreen",
             "B.1.415.1" = "darkseagreen",
             "B.1.438.1" = "darkseagreen",
             "B.1.446" = "darkseagreen",
             "B.1.451" = "darkseagreen",
             "B.1.459" = "darkseagreen",
             "B.1.499" = "darkseagreen",
             "B.1.503" = "darkseagreen1",
             "B.1.527" = "darkseagreen1",
             "B.1.540" = "darkseagreen1",
             "B.1.551" = "darkseagreen1",
             "B.1.556" = "darkseagreen1",
             "B.1.558" = "darkseagreen1",
             "B.1.561" = "darkseagreen1",
             "B.1.564" = "darkseagreen1",
             "B.1.565" = "darkseagreen1",
             "B.1.566" = "darkseagreen1",
             "B.1.567" = "darkseagreen1",
             "B.1.570" = "darkseagreen2",
             "B.1.576" = "darkseagreen2",
             "B.1.577" = "darkseagreen2",
             "B.1.578" = "darkseagreen2",
             "B.1.580" = "darkseagreen2",
             "B.1.582" = "darkseagreen2",
             "B.1.588" = "darkseagreen2",
             "B.1.595" = "darkseagreen2",
             "B.1.596" = "darkseagreen2",
             "B.1.599" = "darkseagreen2",
             "B.1.609" = "darkslateblue",
             "B.1.610" = "darkslateblue",
             "B.1.612" = "darkslateblue",
             "B.1.625" = "darkslateblue",
             "B.1.627" = "darkslateblue",
             "B.1.631" = "darkslateblue",
             "B.1.632" = "darkslateblue",
             "B.1.634" = "darkslateblue",
             "B.1.635" = "darkslateblue",
             "B.1.94" = "lavenderblush",
             "B.15" = "lavenderblush",
             "B.55" = "lavenderblush",
             "B.57" = "lavenderblush",
             "C.23" = "lavenderblush",
             "D.2" = "lavenderblush",
             "None" = "black",
             "P.2" = "lavenderblush",
             "Otros" = "magenta"
             ),
    labels=nameLinage,drop = FALSE )
#scale_fill_manual(breaks=c("B.1.1.7","B.1.1.519","B.1.351","B.1.351.1","B.1.427","B.1.429","B.1.526","B.1.526.2","B.1.617.1","B.1.617.2","C.37","P.1","P.1.1","P.2","P.3","Otros"),labels=c("B.1.1.7","B.1.1.519","B.1.351","B.1.351.1","B.1.427","B.1.429","B.1.526","B.1.526.2","B.1.617.1","B.1.617.2","C.37","P.1","P.1.1","P.2","P.3","Otros") )
GraphVariant<- ggpar(GraphVariant,font.label = list(size=18, color = "grey27"),xlab="",ylab="Abundancia relativa")+
  theme_bw() + scale_y_continuous(expand = c(0, 0))+
  scale_x_date(date_breaks = "2 month",date_labels="%b/%y",date_minor_breaks="2 month", limits = as.Date(c(FechaInicio,Fecha)),expand = c(0,0))+
  theme(panel.border = element_blank(),axis.title=element_text(size=18),axis.text= element_text(size=16,margin=margin(c(0,0,0,-5))))+
  theme( plot.title = element_text(hjust = 0.5, size=20,face = "bold", margin=margin(c(0,0,4,0)))) +
  ggtitle("Variantes en la Republica Mexicana")+
  theme(panel.grid.major = element_blank()) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.box.margin=margin(-12,-5,-5,-12),
    legend.key.size = unit(1.5, "cm"))+
  labs(x = paste("n =",nrow(Variant),sep=" "))+
  theme(plot.margin = unit(c(0,0,0,0), "lines"))
#tiff("Test/12Jul2021_VOC_VOI_Total.tiff", width=17, height=11, units="cm", res=300,compression ="lzw")
#jpeg("Test/Regiones/13_Set_VOC_VOI_Total.jpg", width=90, height=40, units="cm", res = 300)
jpeg("Test/Regiones/Test.jpg", width=90, height=40, units="cm", res = 300)
GraphVariant
dev.off()

####################################

####################################
FechaInicio<-'2021-01-01'
Fecha<-'2021-08-31'
test$F<-1
#Obtener todos los registros
Variant<-test
Variant$Fac<-factor(Variant$Lineage,levels=nameLinage)
#Sacar cuantos registros hay por mes a nivel nacional
Ttm<-ddply(Variant, .(M), summarize,  sum(F))
print (Ttm)
GraphVariant<-ggplot(Variant, aes(Collection.date, after_stat(count),fill = Fac)) +
  geom_density(position = "fill",adjust=2.5) +
  scale_fill_manual(breaks=
                      nameLinage,
                    values=c("Alfa" = "darkorange",
                             "Beta" = "indianred3",
                             "Gamma" = "green3",
                             "Delta" = "gold3",
                             "Iota" = "steelblue1",
                             "Kappa" = "cyan2",
                             "Lambda" = "blueviolet",
                             "Mu" = "azure4",
                             "A" = "ivory1",
                             "B.1.1" = "blue",
                             "None" = "black",
                             "B.1.243" = "bisque1",
                             "B.1.609" = "darkkhaki",
                             "B.1.561" = "darkslategray4",
                             "B.1.427" = "darkolivegreen2",
                             "B.1.429" = "darkolivegreen3",
                             "B.1.2" = "mediumpurple1",
                             "B.1.619" = "aquamarine1",
                             "B.1.620" = "aquamarine2",
                             "B.1.628" = "cornflowerblue",
                             "B.1.1.318" = "lightseagreen",
                             "B.1.1.519" = "indianred1",
                             "B.1.1.222" = "firebrick3",
                             "B.1" = "lightpink",
                             
                             
                             "A.1" = "white",
                             "A.2" = "white",
                             "A.2.3" = "white",
                             "A.2.5" = "white",
                             "A.2.5.1" = "white",
                             "A.2.5.2" = "white",
                             "A.3" = "white",
                             "A.5" = "white",
                             "AZ.3" = "lightpink4",
                             "B" = "tan",
                             "B.1.1.10" = "gray21",
                             "B.1.1.101" = "gray21",
                             "B.1.1.117" = "gray21",
                             "B.1.1.128" = "gray21",
                             "B.1.1.129" = "gray21",
                             "B.1.1.133" = "gray21",
                             "B.1.1.161" = "gray21",
                             "B.1.1.182" = "gray21",
                             "B.1.1.198" = "gray21",
                             "B.1.1.205" = "gray61",
                             "B.1.1.207" = "gray61",
                             "B.1.1.220" = "gray61",
                             "B.1.1.231" = "gray61",
                             "B.1.1.236" = "gray61",
                             "B.1.1.244" = "gray61",
                             "B.1.1.274" = "gray61",
                             "B.1.1.28" = "gray61",
                             "B.1.1.285" = "gray61",
                             "B.1.1.307" = "pink2",
                             "B.1.1.310" = "pink2",
                             "B.1.1.311" = "pink2",
                             "B.1.1.316" = "pink2",
                             "B.1.1.317" = "pink2",
                             "B.1.1.322" = "pink2",
                             "B.1.1.329" = "pink2",
                             "B.1.1.33" = "lightsteelblue1",
                             "B.1.1.334" = "lightsteelblue1",
                             "B.1.1.344" = "lightsteelblue1",
                             "B.1.1.348" = "lightsteelblue1",
                             "B.1.1.362" = "lightsteelblue1",
                             "B.1.1.398" = "lightsteelblue1",
                             "B.1.1.343" = "skyblue",
                             "B.1.1.416" = "skyblue",
                             "B.1.1.418" = "skyblue",
                             "B.1.1.422" = "skyblue",
                             "B.1.1.432" = "skyblue",
                             "B.1.1.434" = "skyblue",
                             "B.1.1.47" = "skyblue",
                             "B.1.1.512" = "aliceblue",
                             "B.1.1.517" = "aliceblue",
                             "B.1.1.518" = "aliceblue",
                             "B.1.1.526" = "aliceblue",
                             "B.1.1.70" = "aliceblue",
                             "B.1.1.71" = "aliceblue",
                             "B.1.1.8" = "aliceblue",
                             "B.1.1.93" = "aliceblue",
                             "B.1.105" = "steelblue",
                             "B.1.111" = "steelblue",
                             "B.1.119" = "steelblue",
                             "B.1.126" = "steelblue",
                             "B.1.127" = "steelblue",
                             "B.1.128" = "steelblue",
                             "B.1.160.32" = "steelblue",
                             "B.1.165" = "steelblue",
                             "B.1.177" = "steelblue",
                             "B.1.177.53" = "steelblue",
                             "B.1.179" = "steelblue",
                             "B.1.189" = "steelblue",
                             "B.1.201" = "pink3",
                             "B.1.206" = "pink3",
                             "B.1.208" = "pink3",
                             "B.1.214.3" = "pink3",
                             "B.1.221" = "pink3",
                             "B.1.229" = "pink3",
                             "B.1.232" = "pink3",
                             "B.1.234" = "pink3",
                             "B.1.236" = "pink3",
                             "B.1.239" = "pink3",
                             "B.1.240" = "pink3",
                             "B.1.241" = "pink3",
                             "B.1.243" = "pink3",
                             "B.1.243.1" = "pink3",
                             "B.1.245" = "pink3",
                             "B.1.258.21" = "pink3",
                             "B.1.267" = "pink3",
                             "B.1.319" = "pink4",
                             "B.1.320" = "pink4",
                             "B.1.324" = "pink4",
                             "B.1.333" = "pink4",
                             "B.1.346" = "pink4",
                             "B.1.349" = "pink4",
                             "B.1.36.10" = "pink4",
                             "B.1.36.31" = "pink4",
                             "B.1.361" = "pink4",
                             "B.1.366" = "pink4",
                             "B.1.369" = "pink4",
                             "B.1.375" = "pink4",
                             "B.1.393" = "pink4",
                             "B.1.395" = "pink4",
                             "B.1.396" = "pink4",
                             "B.1.397" = "pink4",
                             "B.1.398" = "pink4",
                             "B.1.399" = "pink4",
                             "B.1.400" = "darkseagreen",
                             "B.1.400.1" = "darkseagreen",
                             "B.1.402" = "darkseagreen",
                             "B.1.404" = "darkseagreen",
                             "B.1.405" = "darkseagreen",
                             "B.1.409" = "darkseagreen",
                             "B.1.411" = "darkseagreen",
                             "B.1.415" = "darkseagreen",
                             "B.1.415.1" = "darkseagreen",
                             "B.1.438.1" = "darkseagreen",
                             "B.1.446" = "darkseagreen",
                             "B.1.451" = "darkseagreen",
                             "B.1.459" = "darkseagreen",
                             "B.1.499" = "darkseagreen",
                             "B.1.503" = "darkseagreen1",
                             "B.1.527" = "darkseagreen1",
                             "B.1.540" = "darkseagreen1",
                             "B.1.551" = "darkseagreen1",
                             "B.1.556" = "darkseagreen1",
                             "B.1.558" = "darkseagreen1",
                             "B.1.561" = "darkseagreen1",
                             "B.1.564" = "darkseagreen1",
                             "B.1.565" = "darkseagreen1",
                             "B.1.566" = "darkseagreen1",
                             "B.1.567" = "darkseagreen1",
                             "B.1.570" = "darkseagreen2",
                             "B.1.576" = "darkseagreen2",
                             "B.1.577" = "darkseagreen2",
                             "B.1.578" = "darkseagreen2",
                             "B.1.580" = "darkseagreen2",
                             "B.1.582" = "darkseagreen2",
                             "B.1.588" = "darkseagreen2",
                             "B.1.595" = "darkseagreen2",
                             "B.1.596" = "darkseagreen2",
                             "B.1.599" = "darkseagreen2",
                             "B.1.609" = "darkslateblue",
                             "B.1.610" = "darkslateblue",
                             "B.1.612" = "darkslateblue",
                             "B.1.625" = "darkslateblue",
                             "B.1.627" = "darkslateblue",
                             "B.1.631" = "darkslateblue",
                             "B.1.632" = "darkslateblue",
                             "B.1.634" = "darkslateblue",
                             "B.1.635" = "darkslateblue",
                             "B.1.94" = "lavenderblush",
                             "B.15" = "lavenderblush",
                             "B.55" = "lavenderblush",
                             "B.57" = "lavenderblush",
                             "C.23" = "lavenderblush",
                             "D.2" = "lavenderblush",
                             "None" = "black",
                             "P.2" = "lavenderblush",
                             "Otros" = "magenta"
                    ),
                    labels=nameLinage,drop = FALSE )
#scale_fill_manual(breaks=c("B.1.1.7","B.1.1.519","B.1.351","B.1.351.1","B.1.427","B.1.429","B.1.526","B.1.526.2","B.1.617.1","B.1.617.2","C.37","P.1","P.1.1","P.2","P.3","Otros"),labels=c("B.1.1.7","B.1.1.519","B.1.351","B.1.351.1","B.1.427","B.1.429","B.1.526","B.1.526.2","B.1.617.1","B.1.617.2","C.37","P.1","P.1.1","P.2","P.3","Otros") )
GraphVariant<- ggpar(GraphVariant,font.label = list(size=18, color = "grey27"),xlab="",ylab="Abundancia relativa")+
  theme_bw() + scale_y_continuous(expand = c(0, 0))+
  scale_x_date(date_breaks = "1 month",date_labels="%b/%y",date_minor_breaks="1 month", limits = as.Date(c(FechaInicio,Fecha)),expand = c(0,0))+
  theme(panel.border = element_blank(),axis.title=element_text(size=18),axis.text= element_text(size=16,margin=margin(c(0,0,0,-5))))+
  theme( plot.title = element_text(hjust = 0.5, size=20,face = "bold", margin=margin(c(0,0,4,0)))) +
  ggtitle("Variantes en la Republica Mexicana")+
  theme(panel.grid.major = element_blank()) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.box.margin=margin(-12,-5,-5,-12),
    legend.key.size = unit(1.5, "cm"))+
  labs(x = paste("n =",nrow(Variant),sep=" "))+
  theme(plot.margin = unit(c(0,0,0,0), "lines"))
#tiff("Test/12Jul2021_VOC_VOI_Total.tiff", width=17, height=11, units="cm", res=300,compression ="lzw")
#jpeg("Test/Regiones/13_Set_VOC_VOI_Total.jpg", width=90, height=40, units="cm", res = 300)
jpeg("Test/Regiones/Test2.jpg", width=90, height=40, units="cm", res = 300)
GraphVariant
dev.off()




# library("ggplot2")
# library("ggpubr")
# library(dplyr)
# library(plyr)
# setwd("C:/Users/beyon/Dropbox/Difusion-CoViGEn/")
# 
# #Tiene dos entradas.
# #Archivo separado por tabulador con cuatro columnas: 1) Accession ID, 2)Collection date, 3) Lineage y 4) Location
# #linCMX= read.table(file= 'gisaid_hcov-19_1Nov20-18Jun21_Metadata.txt', header=T, sep ='\t' )
# #linCMX= read.table(file= 'gisaid_hcov-19_1Dic21-18Jun21_Metadata.txt', header=T, sep ='\t' )
# linCMX= read.table(file= 'Metadata/gisaid_hcov-19_06Septiembre2021_Metadata_TEST.txt', header=T, sep ='\t', stringsAsFactors = FALSE )
# 
# #Archivo de dos columnas separadas por tabulador: 1)Estado y 2) Region
# EstReg= read.table(file= 'Metadata/EstadoRegion_v3.txt', header=T, sep ='\t', stringsAsFactors = FALSE )
# 
# #Poner la fecha en formato y otra variable que sea por mes
# linCMX$Collection.date<- as.Date(linCMX$Collection.date, "%d/%m/%Y")
# linCMX$M<-as.factor(format(as.Date(linCMX$Collection.date), "%Y-%m"))
# 
# #Cambiar Mexico City a español
# Cambiar<-(linCMX$Location %in% c("Mexico City","Distrito Federal"))
# linCMX$Location[Cambiar]<-"Ciudad de Mexico"
# #Cambiar State of Mexico a español
# Cambiar<-(linCMX$Location %in% c("State of Mexico"))
# linCMX$Location[Cambiar]<-"Estado de Mexico"
# #Poner a cada registro a que región pertenece
# linCMX$Region<-linCMX$Location
# for (i in 1:32){
#   linCMX[linCMX$Location==EstReg[i,1],]$Region<-EstReg[i,2] 
# }
# FechaInicio<-'2020-02-27'
# Fecha<-'2021-08-28'
# linCMX<-linCMX[(linCMX$Collection.date>FechaInicio & linCMX$Collection.date<Fecha),]
# 
# test <- linCMX
# test$Lineage[test$Lineage%in% c('B.1.1.7', 'Q.1', 'Q.3')] <- 'Alfa'
# test$Lineage[test$Lineage %in% c('B.1.351', 'B.1.351.1')] <- 'Beta'
# test$Lineage[test$Lineage %in% c("P.1",
#                                  "P.1.1", 
#                                  "P.1.2",
#                                  "P.1.4",
#                                  "P.1.7",
#                                  "P.1.8")] <- 'Gamma'
# test$Lineage[test$Lineage %in% c("B.1.617.2",
#                                  "AY.2", 
#                                  "AY.3", 
#                                  "AY.4", 
#                                  "AY.5",
#                                  "AY.6",
#                                  "AY.9", 
#                                  "AY.10", 
#                                  "AY.11", 
#                                  "AY.12")] <- 'Delta'
# test$Lineage[test$Lineage == 'B.1.526'] <- 'Iota'
# test$Lineage[test$Lineage == 'B.1.617.1'] <- 'Kappa'
# test$Lineage[test$Lineage %in% c("C.37",
#                                  "C.37.1")] <- 'Lambda'
# test$Lineage[test$Lineage %in% c("B.1.621",
#                                  "B.1.621.1")] <- 'Mu'
# 
# test$Lineage[test$Lineage %in% c("A",
#                                  "A.1",
#                                  "A.2",
#                                  "A.2.3",
#                                  "A.2.5",
#                                  "A.2.5.1",
#                                  "A.2.5.2",
#                                  "A.3",
#                                  "A.5")] <- 'A'
# 
# test$Lineage[test$Lineage %in% c("B.1.1", "B.1.1.10", "B.1.1.101", "B.1.1.117", "B.1.1.128", 
#                                  "B.1.1.129", "B.1.1.133", "B.1.1.161", "B.1.1.182", 
#                                  "B.1.1.198", "B.1.1.205", "B.1.1.207", "B.1.1.220", 
#                                  "B.1.1.231", "B.1.1.236", "B.1.1.244", "B.1.1.274", 
#                                  "B.1.1.28", "B.1.1.285", "B.1.1.307", "B.1.1.310", 
#                                  "B.1.1.311", "B.1.1.316", "B.1.1.317", "B.1.1.322", 
#                                  "B.1.1.329", "B.1.1.33", "B.1.1.334", "B.1.1.344", 
#                                  "B.1.1.348", "B.1.1.362", "B.1.1.398", "B.1.1.403", 
#                                  "B.1.1.416", "B.1.1.418", "B.1.1.422", "B.1.1.432", 
#                                  "B.1.1.434", "B.1.1.47", "B.1.1.512", "B.1.1.517", 
#                                  "B.1.1.518", "B.1.1.526", "B.1.1.70", "B.1.1.71",
#                                  "B.1.1.8", "B.1.1.93")] <- 'B.1.1'
# 
# test$Lineage[!test$Lineage %in% c("Alfa", "Beta", "Gamma", "Delta", "Iota",
#                                   "Kappa", "Lambda", "Mu", "A", "B.1.1", "None",
#                                   "B.1.628", "B.1.609", "B.1.561", "B.1.427", "B.1.429",
#                                   "B.1.1.519", "B.1.619", "B.1.620", "B.1.628",
#                                   "B.1.243", "B.1.2", "B.1.1.318", "B.1.1.222", "B.1")] <- "Otros"
# 
# 
# 
# 


#########################
# linajes <- c("B.1.1.7",
#              
#              "B.1.351",
#              
#              "P.1",
#              "P.1.1", 
#              "P.1.2",
#              "P.1.4",
#              "P.1.7",
#              "P.1.8",
#              
#              "B.1.617.2",
#              "AY.2", 
#              "AY.3", 
#              "AY.4", 
#              "AY.5",
#              "AY.6",
#              "AY.9", 
#              "AY.10", 
#              "AY.11", 
#              "AY.12",
#              
#              "B.1.526",
#              
#              "B.1.617.1",
#         
#              "C.37",
#              "C.37.1",
#              
#              "B.1.621",
#              "B.1.621.1",
#              
#              "B.1.1.318",
#              "B.1.1.519",
#              "B.1.427",
#              "B.1.429"
#              )
# 
# #Poner la fecha en formato y otra variable que sea por mes
# linCMX$Collection.date<- as.Date(linCMX$Collection.date, "%d/%m/%Y")
# linCMX$M<-as.factor(format(as.Date(linCMX$Collection.date), "%Y-%m"))
# #Poner como "Otros" a todos los regristros que no estan en esta lista
# Cambiar<-!(linCMX$Lineage %in% linajes)
# linCMX$Lineage[Cambiar]<-"Otros"
# 
# nameLinage<-c("Alfa",
#               "Beta",
#               "Gamma",
#               "Delta",
#               "Iota",
#               "Kappa",
#               "Lambda",
#               "Mu",
#               "B.1.1.318",
#               "B.1.1.519",
#               "B.1.427",
#               "B.1.429",
#               "Otros"
#               )
# #Cambiar Mexico City a español
# Cambiar<-(linCMX$Location %in% c("Mexico City","Distrito Federal"))
# linCMX$Location[Cambiar]<-"Ciudad de Mexico"
# #Cambiar State of Mexico a español
# Cambiar<-(linCMX$Location %in% c("State of Mexico"))
# linCMX$Location[Cambiar]<-"Estado de Mexico"
# #Poner a cada registro a que región pertenece
# linCMX$Region<-linCMX$Location
# for (i in 1:32){
#   linCMX[linCMX$Location==EstReg[i,1],]$Region<-EstReg[i,2] 
# }
# FechaInicio<-'2020-02-27'
# Fecha<-'2021-08-28'
# linCMX<-linCMX[(linCMX$Collection.date>FechaInicio & linCMX$Collection.date<Fecha),]
# 
# 
# test <- linCMX
# test$Lineage[test$Lineage == 'B.1.1.7'] <- 'Alfa'
# test$Lineage[test$Lineage %in% c('B.1.351', 'B.1.351.1')] <- 'Beta'
# test$Lineage[test$Lineage %in% c("P.1",
#                                "P.1.1", 
#                                "P.1.2",
#                                "P.1.4",
#                                "P.1.7",
#                                "P.1.8")] <- 'Gamma'
# test$Lineage[test$Lineage %in% c("B.1.617.2",
#                                  "AY.2", 
#                                  "AY.3", 
#                                  "AY.4", 
#                                  "AY.5",
#                                  "AY.6",
#                                  "AY.9", 
#                                  "AY.10", 
#                                  "AY.11", 
#                                  "AY.12")] <- 'Delta'
# test$Lineage[test$Lineage == 'B.1.526'] <- 'Iota'
# test$Lineage[test$Lineage == 'B.1.617.1'] <- 'Kappa'
# test$Lineage[test$Lineage %in% c("C.37",
#                                  "C.37.1")] <- 'Lambda'
# test$Lineage[test$Lineage %in% c("B.1.621",
#                                  "B.1.621.1")] <- 'Mu'
################

# c("Alfa" = "darkorange",
# "Beta" = "indianred3",
# "Gamma" = "green3",
# "Delta" = "gold3",
# "Iota" = "steelblue1",
# "Kappa" = "cyan2",
# "Lambda" = "blueviolet",
# "Mu" = "azure4",
# "A" = "coral2",
# "B.1.1" = "blue",
# "None" = "black",
# "B.1.609" = "darkslategrey",
# "B.1.561" = "darkslategray4",
# "B.1.427" = "darkolivegreen2",
# "B.1.429" = "darkolivegreen3",
# "B.1.2" = "mediumpurple1",
# "B.1.619" = "aquamarine1",
# "B.1.620" = "aquamarine2",
# "B.1.628" = "azure2",
# "B.1.1.318" = "lightseagreen",
# "B.1.1.519" = "indianred1",
# "B.1.1.222" = "indianred3",
# "B.1" = "mediumorchid1",
# "Otros" = "magenta")



#"Alfa", "Beta", "Gamma", "Delta", "Iota", "Kappa", "Lambda", "Mu", "A", "B.1.1", "None", "B.1.609", "B.1.561", "B.1.427", "B.1.429", "B.1.2", "B.1.619", "B.1.620", "B.1.628", "B.1.1.318", "B.1.1.519", "B.1.1.222", "B.1", "Otros”
##Grafica a nivel nacional

# nameLinage <- c("Alfa", 
#                 "Beta", 
#                 "Gamma", 
#                 "Delta", 
#                 "Iota", 
#                 "Kappa", 
#                 "Lambda", 
#                 "Mu", 
#                 "A", 
#                 "B.1.1", 
#                 "None", 
#                 "B.1.243",
#                 "B.1.609", 
#                 "B.1.561", 
#                 "B.1.427", 
#                 "B.1.429", 
#                 "B.1.2", 
#                 "B.1.619", 
#                 "B.1.620", 
#                 "B.1.628", 
#                 "B.1.1.318", 
#                 "B.1.1.519", 
#                 "B.1.1.222", 
#                 "B.1", 
#                 "Otros")


# FechaInicio<-'2020-02-27'
# Fecha<-'2021-08-16'
# test$F<-1
# #Obtener todos los registros
# Variant<-test
# Variant$Fac<-factor(Variant$Lineage,levels=nameLinage)
# #Sacar cuantos registros hay por mes a nivel nacional
# Ttm<-ddply(Variant, .(M), summarize,  sum(F))
# print (Ttm)
# GraphVariant<-ggplot(Variant, aes(Collection.date, after_stat(count),fill = Fac)) +
#   geom_density(position = "fill",adjust=2.5) +
#   scale_fill_manual(breaks=c(
#                               "Alfa", 
#                               "Beta", 
#                               "Gamma", 
#                               "Delta", 
#                               "Iota", 
#                               "Kappa", 
#                               "Lambda", 
#                               "Mu", 
#                               "A", 
#                               "B.1.1", 
#                               "None", 
#                               "B.1.243",
#                               "B.1.609", 
#                               "B.1.561", 
#                               "B.1.427", 
#                               "B.1.429", 
#                               "B.1.2", 
#                               "B.1.619", 
#                               "B.1.620", 
#                               "B.1.628", 
#                               "B.1.1.318", 
#                               "B.1.1.519", 
#                               "B.1.1.222", 
#                               "B.1", 
#                               "Otros"),
#                     values=c("Alfa" = "darkorange",
#                              "Beta" = "indianred3",
#                              "Gamma" = "green3",
#                              "Delta" = "gold3",
#                              "Iota" = "steelblue1",
#                              "Kappa" = "cyan2",
#                              "Lambda" = "blueviolet",
#                              "Mu" = "azure4",
#                              "A" = "ivory1",
#                              "B.1.1" = "blue",
#                              "None" = "black",
#                              "B.1.243" = "bisque1",
#                              "B.1.609" = "darkkhaki",
#                              "B.1.561" = "darkslategray4",
#                              "B.1.427" = "darkolivegreen2",
#                              "B.1.429" = "darkolivegreen3",
#                              "B.1.2" = "mediumpurple1",
#                              "B.1.619" = "aquamarine1",
#                              "B.1.620" = "aquamarine2",
#                              "B.1.628" = "cornflowerblue",
#                              "B.1.1.318" = "lightseagreen",
#                              "B.1.1.519" = "indianred1",
#                              "B.1.1.222" = "firebrick3",
#                              "B.1" = "lightpink",
#                              "Otros" = "magenta"),
#                     labels=c("Alfa", 
#                              "Beta", 
#                              "Gamma", 
#                              "Delta", 
#                              "Iota", 
#                              "Kappa", 
#                              "Lambda", 
#                              "Mu", 
#                              "A*", 
#                              "B.1.1**", 
#                              "None", 
#                              "B.1.243",
#                              "B.1.609", 
#                              "B.1.561", 
#                              "B.1.427", 
#                              "B.1.429", 
#                              "B.1.2", 
#                              "B.1.619", 
#                              "B.1.620", 
#                              "B.1.628", 
#                              "B.1.1.318", 
#                              "B.1.1.519", 
#                              "B.1.1.222", 
#                              "B.1", 
#                              "Otros"),drop = FALSE )
# #scale_fill_manual(breaks=c("B.1.1.7","B.1.1.519","B.1.351","B.1.351.1","B.1.427","B.1.429","B.1.526","B.1.526.2","B.1.617.1","B.1.617.2","C.37","P.1","P.1.1","P.2","P.3","Otros"),labels=c("B.1.1.7","B.1.1.519","B.1.351","B.1.351.1","B.1.427","B.1.429","B.1.526","B.1.526.2","B.1.617.1","B.1.617.2","C.37","P.1","P.1.1","P.2","P.3","Otros") )
# GraphVariant<- ggpar(GraphVariant,font.label = list(size=9, color = "grey27"),xlab="",ylab="Abundancia relativa")+
#   theme_bw() + scale_y_continuous(expand = c(0, 0))+
#   scale_x_date(date_breaks = "3 month",date_labels="%b/%y",date_minor_breaks="3 month", limits = as.Date(c(FechaInicio,Fecha)),expand = c(0,0))+
#   theme(panel.border = element_blank(),axis.title=element_text(size=16.5),axis.text= element_text(size=9,margin=margin(c(0,0,0,-5))))+
#   theme( plot.title = element_text(hjust = 0.5, size=10,face = "bold", margin=margin(c(0,0,4,0)))) +
#   ggtitle("Variantes en la Republica Mexicana")+
#   theme(panel.grid.major = element_blank()) +
#   theme(
#     legend.title = element_blank(),
#     legend.text = element_text(size = 9),
#     legend.box.margin=margin(-10,-5,-5,-10),
#     legend.key.size = unit(0.4, "cm"))+
#   labs(x = paste("n =",nrow(Variant),sep=" "))+
#   theme(plot.margin = unit(c(0,0,0,0), "lines"))
# #tiff("Test/12Jul2021_VOC_VOI_Total.tiff", width=17, height=11, units="cm", res=300,compression ="lzw")
# jpeg("Test/Regiones/TEST_VOC_VOI_Total.jpg", width=25, height=14, units="cm", res = 300)
# GraphVariant
# dev.off()


# ##Grafica a nivel nacional
# FechaInicio<-'2020-02-27'
# Fecha<-'2021-08-16'
# test$F<-1
# #Obtener todos los registros
# Variant<-test
# Variant$Fac<-factor(Variant$Lineage,levels=nameLinage)
# #Sacar cuantos registros hay por mes a nivel nacional
# Ttm<-ddply(Variant, .(M), summarize,  sum(F))
# print (Ttm)
# GraphVariant<-ggplot(Variant, aes(Collection.date, after_stat(count),fill = Fac)) +
#   geom_density(position = "fill",adjust=3) +
#   scale_fill_manual(breaks=c("Alfa",
#                              "Beta",
#                              "Gamma",
#                              "Delta",
#                              "Iota",
#                              "Kappa",
#                              "Lambda",
#                              "Mu",
#                              "B.1.1.318",
#                              "B.1.1.519",
#                              "B.1.427",
#                              "B.1.429",
#                              "Otros"),
#                     values=c("Alfa"="darkorange",
#                              "Beta"="indianred3",
#                              "Gamma"="green3",
#                              "Delta"="gold3",
#                              "Iota"="steelblue1",
#                              "Kappa"="cyan2",
#                              "Lambda"="blueviolet",
#                              "Mu"="azure4",
#                              "B.1.1.318" = "lightseagreen",
#                              "B.1.1.519"="indianred1",
#                              "B.1.427"="darkolivegreen2",
#                              "B.1.429"="darkolivegreen3",
#                              "Otros"="magenta"),
#                     labels=c("Alfa",
#                              "Beta",
#                              "Gamma",
#                              "Delta",
#                              "Iota",
#                              "Kappa",
#                              "Lambda",
#                              "Mu",
#                              "B.1.1.318",
#                              "B.1.1.519",
#                              "B.1.427",
#                              "B.1.429",
#                              "Otros"),drop = FALSE )
# #scale_fill_manual(breaks=c("B.1.1.7","B.1.1.519","B.1.351","B.1.351.1","B.1.427","B.1.429","B.1.526","B.1.526.2","B.1.617.1","B.1.617.2","C.37","P.1","P.1.1","P.2","P.3","Otros"),labels=c("B.1.1.7","B.1.1.519","B.1.351","B.1.351.1","B.1.427","B.1.429","B.1.526","B.1.526.2","B.1.617.1","B.1.617.2","C.37","P.1","P.1.1","P.2","P.3","Otros") )
# GraphVariant<- ggpar(GraphVariant,font.label = list(size=9, color = "grey27"),xlab="",ylab="Abundancia relativa")+
#   theme_bw() + scale_y_continuous(expand = c(0, 0))+
#   scale_x_date(date_breaks = "2 month",date_labels="%b/%y",date_minor_breaks="2 month", limits = as.Date(c(FechaInicio,Fecha)),expand = c(0,0))+
#   theme(panel.border = element_blank(),axis.title=element_text(size=16.5),axis.text= element_text(size=9,margin=margin(c(0,0,0,-5))))+
#   theme( plot.title = element_text(hjust = 0.5, size=10,face = "bold", margin=margin(c(0,0,4,0)))) +
#   ggtitle("Variantes en la Republica Mexicana")+
#   theme(panel.grid.major = element_blank()) +
#   theme(
#     legend.title = element_blank(),
#     legend.text = element_text(size = 9),
#     legend.box.margin=margin(-10,-5,-5,-10),
#     legend.key.size = unit(0.4, "cm"))+
#   labs(x = paste("n =",nrow(Variant),sep=" "))+
#   theme(plot.margin = unit(c(0,0,0,0), "lines"))
# #tiff("Test/12Jul2021_VOC_VOI_Total.tiff", width=17, height=11, units="cm", res=300,compression ="lzw")
# jpeg("Test/Regiones/TEST_VOC_VOI_Total.jpg", width=25, height=14, units="cm", res = 300)
# GraphVariant
# dev.off()


