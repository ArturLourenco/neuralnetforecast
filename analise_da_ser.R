###############################################################################################################################
################################################# Clear Console and Memory ####################################################
###############################################################################################################################

gc(TRUE) #garbage colector da RAM
rm(list = ls()) #limpar memria do GE
dev.off() # limpar os plots
cat("\014") #limpar console

###############################################################################################################################
########################################### Load necessary libraries and set work directory####################################
###############################################################################################################################

list.of.packages <- c("hydroTSM", "lubridate",'zoo','openxlsx')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
require(hydroTSM)
require(lubridate)
require(openxlsx)
require(reshape2)
Sys.setenv(R_ZIPCMD= "C:/Rtools/bin/zip")
setwd("G:/My Drive/Pesquisa/Coordenador de Projetos Pesquisa 2017/Pr�-tratamento dos dados") ## mudar a o local aonde est�o os arquivos
load("Dados_Zoo.RData") ## Carregar dados gerados pelo hydroweb
path<- getwd();

###############################################################################################################################
########################################### An�lise e pr� tratamento dos dados ################################################
###############################################################################################################################

setAs("character","myDate", function(from) as.Date(from, format="%d/%m/%Y") )

volumejatoba2<- read.csv2("G:/My Drive/Pesquisa/Coordenador de Projetos Pesquisa 2018/Dados Caatinga/volumejatoba.csv",header = T,
                         sep = ",",dec = ",",
                         colClasses=c("character","myDate","character","character","character"))

volumejatoba2$Volume....<- as.numeric(gsub(",",".",volumejatoba2$Volume....))
volumejatoba2$'Volume..m³.'<- as.numeric(gsub(',','.',gsub("\\.","",volumejatoba2$'Volume..m³.')))
volumejatoba2$'Afluência.Defluência..m².'<- as.numeric(gsub(',','.',gsub("\\.","",volumejatoba2$'Afluência.Defluência..m².')))

princesa2017<- read.csv2("princesa2017.csv",header = T,
                         sep = ",",dec = ",",
                         colClasses=c("character","character","myDate","character"))

princesa2017$Precipita�.�.o..mm.<- as.numeric(gsub(",",".",princesa2017$Precipita�.�.o..mm.))

princesa2018<- read.csv2("princesa2018.csv",header = T,
                         sep = ",",dec = ",",
                         colClasses=c("character","character","myDate","character"))

princesa2018$Precipita�.�.o..mm.<- as.numeric(gsub(",",".",princesa2018$Precipita�.�.o..mm.))

princesa2019<- read.csv2("princesa2019.csv",header = T,
                         sep = ",",dec = ",",
                         colClasses=c("character","character","myDate","character"))

princesa2019$Precipita�.�.o..mm.<- as.numeric(gsub(",",".",princesa2019$Precipita�.�.o..mm.))

a2017<-zoo(princesa2017$Precipita�.�.o..mm.,princesa2017$Data.do.registro)

a2018<-zoo(princesa2018$Precipita�.�.o..mm.,princesa2018$Data.do.registro)

a2019<-zoo(princesa2019$Precipita�.�.o..mm.,princesa2019$Data.do.registro)

v_jatobaii_1995_2020<- zoo(volumejatoba2$'Volume..m³.',volumejatoba2$Data.do.registro)
v_jatobaii_1995_2020full<- merge(v_jatobaii_1995_2020,zoo(,seq(start(v_jatobaii_1995_2020),end(v_jatobaii_1995_2020),by="month")), all=TRUE)

v_max<- aggregate(v_jatobaii_1995_2020, format(time(v_jatobaii_1995_2020), "%m/%Y"), max, na.rm = TRUE) # Calcular a m�dias dos dias por m�s

pi_plu_ser_1911_2019<- rbind.zoo(a1911_1991,a1994_2016,a2017,a2018,a2019)
pi_plu_ser_1911_2019_full <- merge(pi_plu_ser_1911_2019,zoo(,seq(start(pi_plu_ser_1911_2019),end(pi_plu_ser_1911_2019),by="day")), all=TRUE)

p_annual<- aggregate(pi_plu_ser_1911_2019_full, format(time(pi_plu_ser_1911_2019_full), "%Y"), sum, na.rm = TRUE) # Calcular a m�dias dos dias por m�s

p_med<- aggregate(pi_plu_ser_1911_2019_full, format(time(pi_plu_ser_1911_2019_full), "%d/%m/%Y"), mean, na.rm = TRUE)

pi_plu_ser_1911_2016<- rbind.zoo(a1911_1991,a1994_2016,a2017,a2018,a2019)
pi_plu_ser_1911_2016_full <- merge(pi_plu_ser_1911_2016,zoo(,seq(start(pi_plu_ser_1911_2016),end(pi_plu_ser_1911_2016),by="day")), all=TRUE)
dip(start(pi_plu_ser_1911_2016),end(pi_plu_ser_1911_2016), out.type = "nmbr")
dip(start(pi_plu_ser_1911_2016_full),end(pi_plu_ser_1911_2016_full), out.type = "nmbr")

# m <- daily2annual(x, FUN=sum, na.rm=TRUE)
# m <- daily2monthly(x, FUN=sum, na.rm=TRUE)
# 
# smry(m)
# 
# hydroplot(x, var.type="Precipitation", main="em Pianc�",pfreq = "ma", from="1985-01-01")

###############################################################################################################################
########################################### An�lise e pr� tratamento dos dados ################################################
###############################################################################################################################


Datazoomonth<- daily2monthly(pi_plu_ser_1911_2019_full, FUN= sum) # totais mensais
Datamatrix <- matrix(Datazoomonth, ncol=12, byrow=TRUE)
monthnames <- c('Jan','Fev','Mar','Abr','Mai','Jun','Jul','Ago','Set','Out','Nov','Dez');
colnames(Datamatrix) <- monthnames
rownames(Datamatrix) <- unique(format(time(Datazoomonth), "%Y"))
print(matrixplot(Datamatrix, ColorRamp="Precipitation", main="Precipita��o Mensal Princesa Isabel, [mm/m�s]"))


###############################################################################################################################
########################################### An�lise e pr� tratamento dos dados ################################################
###############################################################################################################################

dwimatrix<- dwi(pi_plu_ser_1911_2016_full, out.unit="years", dates=1)
matrixplot(as.matrix(dwimatrix), var.type="Days", main="Number of months with info per year")

mm<- aggregate(daily2monthly(pi_plu_ser_1911_2019_full, FUN=sum, na.rm=TRUE), 
          format(time(daily2monthly(pi_plu_ser_1911_2019_full, FUN=sum, na.rm=TRUE)), "%m"), mean) # Calcular as m�dias mensais da s�rie

md<- aggregate(pi_plu_ser_1911_2016_full, format(time(pi_plu_ser_1911_2016_full), "%d"), mean, na.rm = TRUE) # Calcular a m�dias dos dias de cada m�s

ma<- daily2annual(pi_plu_ser_1911_2016_full, FUN=sum, na.rm=TRUE) # Calcular os volumes anuais

mean(daily2annual(pi_plu_ser_1911_2016_full, FUN=sum, na.rm=TRUE)) # Calcula m�dia anual


mdm<- aggregate(pi_plu_ser_1911_2019_full, format(time(pi_plu_ser_1911_2019_full), "%d/%m"), mean, na.rm = TRUE) # Calcular a m�dias dos dias por m�s

a<- as.data.frame(mdm)
a$names <- rownames(a)
a<- cbind(a,substr(a[,2],4,5),substr(a[,2],1,2))
colnames(a)<- c("Valor","Dia_Mes","Mes","Dia")
b<- dcast(a,Dia ~ Mes,value.var = "Valor")
c<-b[,2:ncol(b)]
colnames(c)<- monthnames
matrixplot(as.matrix(c), var.type="Days", main="M�dia de Precipita��o para cada dia e m�s em mm, 1911/2016")


mdm<- aggregate(pi_plu_ser_1911_2016_full, format(time(pi_plu_ser_1911_2016_full), "%d/%m"), mean, na.rm = TRUE) # Calcular a m�dias dos dias por m�s

a<- as.data.frame(mdm)
a$names <- rownames(a)
a<- cbind(a,substr(a[,2],4,5),substr(a[,2],1,2))
colnames(a)<- c("Valor","Dia_Mes","Mes","Dia")
b<- dcast(a,Dia ~ Mes,value.var = "Valor")
c<-b[,2:ncol(b)]
colnames(c)<- monthnames
matrixplot(as.matrix(c), var.type="Days", main="M�dia de Precipita��o para cada dia e m�s em mm, 1911/2016")

mdm<- aggregate(pi_plu_ser_1911_2016_full, format(time(pi_plu_ser_1911_2016_full), "%d/%m"), max, na.rm = TRUE) # Calcular a m�dias dos dias por m�s

a<- as.data.frame(mdm)
a$names <- rownames(a)
a<- cbind(a,substr(a[,2],4,5),substr(a[,2],1,2))
colnames(a)<- c("Valor","Dia_Mes","Mes","Dia")
b<- dcast(a,Dia ~ Mes,value.var = "Valor")
c<-b[,2:ncol(b)]
colnames(c)<- monthnames
matrixplot(as.matrix(c), var.type="Days", main="M�xima de Precipita��o para cada dia e m�s em mm, 1911/2016")

mdm<- aggregate(pi_plu_ser_1911_2016_full, format(time(pi_plu_ser_1911_2016_full), "%d/%m"), min, na.rm = TRUE) # Calcular a m�dias dos dias por m�s

a<- as.data.frame(mdm)
a$names <- rownames(a)
a<- cbind(a,substr(a[,2],4,5),substr(a[,2],1,2))
colnames(a)<- c("Valor","Dia_Mes","Mes","Dia")
b<- dcast(a,Dia ~ Mes,value.var = "Valor")
c<-b[,2:ncol(b)]
colnames(c)<- monthnames
matrixplot(as.matrix(c), var.type="Days", main="M�nima de Precipita��o para cada dia e m�s em mm, 1911/2016")

countday<- function(x){sum((x>0),na.rm = TRUE)}

mdmcount<- aggregate(pi_plu_ser_1911_2016_full, format(time(pi_plu_ser_1911_2016_full), "%d/%m"), countday)

a<- as.data.frame(mdmcount)
a$names <- rownames(a)
a<- cbind(a,substr(a[,2],4,5),substr(a[,2],1,2))
colnames(a)<- c("Valor","Dia_Mes","Mes","Dia")
b<- dcast(a,Dia ~ Mes,value.var = "Valor")
c<-b[,2:ncol(b)]
colnames(c)<- monthnames
matrixplot(as.matrix(c),ColorRamp = "PCPAnomaly", var.type="Days", main="Dias com chuva por dia e m�s em mm, 1911/2016")

mdmmax<- aggregate(pi_plu_ser_1911_2016_full, format(time(pi_plu_ser_1911_2016_full), "%d/%m"), max, na.rm = TRUE)
a<- as.data.frame(mdmmax)
a<- a[order(a$mdmmax,decreasing = TRUE), ,drop = FALSE]
a<- cbind(a,index(a))
a$names <- rownames(a)
# a<- cbind(a,Freq = a$`index(a)`/ (nrow(a)+1))
a<- cbind(a,Freq = a$`index(a)`/ (105+1))
a<- cbind(a,TR = (1/a$Freq) )
a<- cbind(a,Mes = substr(a[,3],4,5),Dia = substr(a[,3],1,2))
colnames(a)<- c("Valor","Ind","Dia_Mes","Freq","TR","Mes","Dia")
b<- dcast(a,Dia ~ Mes,value.var = "TR")
c<-b[,2:ncol(b)]
colnames(c)<- monthnames
matrixplot(as.matrix(c),ColorRamp = "TEMPAnomaly2", var.type="Days", main="Tempo de Retorno em anos para chuvas m�ximas, 1911/2016")


mdm<- aggregate(pi_plu_ser_1911_2019_full, format(time(pi_plu_ser_1911_2019_full), "%m/%Y"), sum, na.rm = TRUE) # Calcular a chuva mensal
a<- as.data.frame(mdm)
a$names <- rownames(a)
a<- cbind(a,substr(a[,2],4,7),substr(a[,2],1,2))
colnames(a)<- c("Valor","Dia_Ano","Ano","Mes")
b<- dcast(a,Ano ~ Mes,value.var = "Valor")
c<-b[,2:ncol(b)]
colnames(c)<- monthnames
rownames(c)<- as.vector(unique(a$Ano))
matrixplot(as.matrix(c[(nrow(c)-10):nrow(c),]), var.type="Days", main="Precipita��o Mensal (mm), 2006-2016 - Princesa Isabel/PB")
e<-rowSums(c[,1:2])

# ggplot(c, aes(rownames(x), colnames(x), fill = c)) + geom_tile() +
#   xlab("X Coordinate (feet)") + ylab("Y Coordinate (feet)") +
#   opts(title = "Surface elevation data") +
#   scale_fill_gradient(limits = c(7000, 10000), low = "black", high = "white") +
#   scale_x_continuous(expand = c(0,0)) +
#   scale_y_continuous(expand = c(0,0))

#################################### Exportar ###################################
# 
# wb2 <- createWorkbook("data")
# addWorksheet(wb2,sheetName = "data");
# writeData(wb2,sheet = "data",x = data,colNames = TRUE, rowNames = TRUE);
# saveWorkbook(wb2, file = paste(path,'data_inputs_and_target.xlsx',sep='/'), overwrite = TRUE) # nome do arquivo de sa�da, ou seja, em uma coluna;
# 
