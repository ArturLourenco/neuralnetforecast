###############################################################################################################################
################################################# Limpar Console e Memória ####################################################
###############################################################################################################################

gc(TRUE) #garbage colector (coletor de lixo) da RAM
rm(list = ls()) #limpar memoria do GE
dev.off() # limpar os plots (gráficos)
cat("\014") #limpar console

mapsd <- function(x) {
  (x - mean(x)) / sd(x)
}

revmapsd <- function(normaldata,origdata) {
  (normaldata * sd(origdata)) + mean(origdata)
}

###############################################################################################################################
########################################### Instalar e carregar pacotes (bibliotecas) necessários ####################################
###############################################################################################################################

list.of.packages <- c("hydroTSM","corrplot","lubridate",'openxlsx','reshape2','plyr',
                      'neuralnet',"stringr","hydroGOF","pracma","RSNNS","doParallel","ggcorrplot",
                      "corrplot") #list pacotes
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])] #instala se não houver
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

###############################################################################################################################
########################################### Definindo diretório de trabalho(pasta) e carregando funções #######################
###############################################################################################################################

Sys.setenv(R_ZIPCMD= "C:/Rtools/bin/zip")
setwd("G:/My Drive/Pesquisa/Coordenador de Projetos Pesquisa 2017/Pré-tratamento dos dados") ## mudar a o local aonde estão os arquivos
path<- getwd();
source("lagpadf.r")
registerDoParallel(cores = 4)

###############################################################################################################################
########################################### Carregando dados externos e definindo dados de entrada para a RNA #################
###############################################################################################################################

load("Dados_Zoo.RData") ## Carregar dados gerados pelo hydroweb

pi_plu_ser_1911_2016<- rbind.zoo(a1911_1991,a1994_2016) #junta as duas séries
pi_plu_ser_1911_2016_full <- merge(pi_plu_ser_1911_2016,zoo(,seq(start(pi_plu_ser_1911_2016),
                                                                 end(pi_plu_ser_1911_2016),by="day")), all=TRUE)
wb <- loadWorkbook("Vol_Jatoba_II.xlsx") # colocar o nome do arquivo do excel
pi_mensal<- daily2monthly(pi_plu_ser_1911_2016_full, FUN=sum, na.rm=TRUE)
vol_jatoba = read.xlsx(wb, sheet = 1,detectDates = TRUE)  
jatoba_ser<- zoo(as.numeric(gsub(",",".",vol_jatoba$`Volume.(%)`)),format(as.Date(vol_jatoba$Data.do.registro,"%d/%m/%Y"),"%Y/%m"))
posto_plu<- window(pi_mensal,index. = index(pi_mensal), as.Date("2008-01-01"),as.Date("2016-12-01")) #delimitar intervalo de chuva
posto_plu_2<- zoo(coredata(posto_plu),format(as.Date(time(posto_plu),"%d/%m/%Y"),"%Y/%m")) #transformar para zoo
EV<- c(207,183,191,169,152,132,143,167,185,208,210,215)
EV<- rep(EV,yip(start(posto_plu),end(posto_plu),out.type = "nmbr"))
EV<- zoo(EV,index(posto_plu_2))
dataset<- merge(posto_plu_2,jatoba_ser,EV,all = FALSE)
dataset<- as.data.frame(dataset)
data<- dataset
rown<- rownames(data)

picos<- data[,2]
vales<- data[,2]

seqds<- cumsum(c(1, sign(diff(data[,2])))) #seríe sitética 

# seqpv<- c(1,rbind(match(findpeaks(picos)[,1],picos),match(abs(findpeaks(-vales)[,1]),vales)),nrow(data))
# seqds<- c()
# 
# for(i in 1:length(seqpv/2)){
#   if(i==length(seqpv/2)){break}
#   ifelse(i %% 2 == 0,seqds[seqpv[i]:seqpv[i+1]]<--1,seqds[seqpv[i]:seqpv[i+1]]<-1)
# }
# 
# seqds[seqpv[length(seqpv)-1]:seqpv[length(seqpv)]]<- 1
# 
# # seqds[match(findpeaks(picos)[,1],picos)]<-2
# # seqds[match(abs(findpeaks(-vales)[,1]),vales)]<-4

picos[-match(findpeaks(picos)[,1],picos)]<-0
vales[-match(abs(findpeaks(-vales)[,1]),vales)]<-0

picos[picos>0]<- max(seqds)
vales[vales>0]<- min(seqds)
# data<- cbind(data,data$jatoba_ser,Picos=picos[dropout=T])
data<- cbind(data,SS=seqds[dropout=T],Picos=picos[dropout=T],Vales=vales[dropout=T])


# data<- cbind(data[,"jatoba_ser"],data[,"jatoba_ser"],data[,"jatoba_ser"],
#              data[,"posto_plu_2"],data[,"posto_plu_2"],data[,"posto_plu_2"],data[,"Com"],data[,"Picos"],data[,"Vales"])

# data<- cbind(data[,"jatoba_ser"],data[,"jatoba_ser"],data[,"jatoba_ser"],
#              data[,"posto_plu_2"],data[,"posto_plu_2"],data[,"Com"],data[,"Picos"],data[,"Vales"])

data<- cbind(data[,"jatoba_ser"],data[,"jatoba_ser"],data[,"jatoba_ser"],
             data[,"posto_plu_2"],data[,"posto_plu_2"],data[,"EV"],data[,"SS"])

# data<- cbind(data[,"jatoba_ser"],data[,"jatoba_ser"],data[,"jatoba_ser"],
#              data[,"posto_plu_2"],data[,"posto_plu_2"],data[,"Com"],data[,"Com"])

# colnames(data)<- c("V_t","V_t_1","V_t_2","p_t","p_t_1","p_t_2","Com","Picos","Vales")
# colnames(data)<- c("V_t","V_t_1","V_t_2","p_t","p_t_2","Com","Picos","Vales")
colnames(data)<- c("V_t","V_t_1","V_t_2","p_t","p_t_2","EV","SS")
# colnames(data)<- c("V_t","V_t_1","V_t_2","p_t","p_t_2","Com_t_1","Com_t_2")
coln<- colnames(data)


###############################################################################################################################
########################################### Definindo os atrasos (lag) para a RNA #############################################
###############################################################################################################################

# data<- cbind(data[,"V_t"],lagpad(data[,"V_t_1"],1),lagpad(data[,"V_t_2"],2),data[,"p_t"],lagpad(data[,"p_t_1"],1)
#              ,lagpad(data[,"p_t_2"],2),data[,"Com"],data[,"Picos"],data[,"Vales"])

# data<- cbind(data[,"V_t"],lagpad(data[,"V_t_1"],1),lagpad(data[,"V_t_2"],2),data[,"p_t"]
#              ,lagpad(data[,"p_t_2"],2),data[,"Com"],data[,"Picos"],data[,"Vales"])

data<- cbind(data[,"V_t"],lagpad(data[,"V_t_1"],1),lagpad(data[,"V_t_2"],2),data[,"p_t"]
             ,lagpad(data[,"p_t_2"],2),lagpad(data[,"EV"],1),lagpad(data[,"SS"],1))

# data<- cbind(data[,"V_t"],lagpad(data[,"V_t_1"],1),lagpad(data[,"V_t_2"],2),data[,"p_t"]
#              ,lagpad(data[,"p_t_2"],2),lagpad(data[,"Com_t_1"],1),lagpad(data[,"Com_t_2"],2))

colnames(data)<-coln
rownames(data)<- rown
data<- as.data.frame(data)

###############################################################################################################################
########################################### Correlação entre as entradas ##################################################
###############################################################################################################################

M <- cor(data[,-1]) # get correlations
row.names(M)<- c("V(t-1)","V(t-2)", "P(t)","P(t-2)", "EV(t-1)","SS")
col.nam(M)<- c("V(t-1)","V(t-2)", "P(t)","P(t-2)", "EV(t-1)","SS")
corrplot(M, method = "circle") #plot matrix

ggcorrplot(M, method = "circle",legend.title = "r",ggtheme = theme_gray) #plot matrix


###############################################################################################################################
########################################### Checar se falta dados nas séries ##################################################
###############################################################################################################################

apply(data,2,function(x) sum(is.na(x)))

###############################################################################################################################
########################################### Definindo o conjunto de dados para treino e teste ##################################
###############################################################################################################################

# index <- sample(1:nrow(data),round(0.75*nrow(data)))
index <- 1:round(0.75*nrow(data))
train <- data[index,]
test <- data[-index,]
# train <- data[1:73,]
# test <- data[75:98,]


###############################################################################################################################
########################################### Modelo linear de comparação ao modelo RNA (não linear) #############################################
###############################################################################################################################

lm.fit <- glm(V_t~., data=train)

###############################################################################################################################
########################################### Previsão do modelo linear #############################################
###############################################################################################################################

pr.lm <- predict(lm.fit,test)
MSE.lm <- sum((pr.lm - test$V_t)^2)/nrow(test) #desempenho do modelo linear (Mean Square Error  = Erro do Quadradro Médio)

###############################################################################################################################
########################################### Normalizado (scalonando) os dados de entrada da RNA #############################################
###############################################################################################################################

# maxs <- apply(data, 2, max) 
# mins <- apply(data, 2, min)
# scaled <- as.data.frame(scale(data, center = mins, scale = maxs - mins))
# scaled <- as.data.frame(scale(data))
scaled<- data
# colnames(scaled)<- coln

###############################################################################################################################
########################################### Definindo os dados de treino e teste escalonados #############################################
###############################################################################################################################

# maxst <- apply(scaled[index,], 2, max)
# minst <- apply(scaled[index,], 2, min)
# maxste <- apply(scaled[-index,], 2, max)
# minste <- apply(scaled[-index,], 2, min)
# train_n <- as.data.frame(scale(scaled[index,], center = minst, scale = maxst - minst))
# test_n <- as.data.frame(scale(scaled[-index,], center = minste, scale = maxste - minste))
# 
# train_ <- normalizeData(train_n)
# test_ <- normalizeData(test_n)

train_ <- normalizeData(scaled[index,])
test_ <- normalizeData(scaled[-index,])

# maxst <- apply(train_n, 2, max)
# minst <- apply(train_n, 2, min)
# maxste <- apply(test_n, 2, max)
# minste <- apply(test_n, 2, min)
# train_ <- as.data.frame(scale(train_n, center = minst, scale = maxst - minst))
# test_ <- as.data.frame(scale(test_n, center = minste, scale = maxste - minste))

# train_ <- scaled[1:73,]
# test_ <- scaled[75:98,]
colnames(train_)<- coln
colnames(test_)<- coln


###############################################################################################################################
########################################### Treinando a Rede Neural Artificial #############################################
###############################################################################################################################


n <- colnames(train_)
f <- as.formula(paste("V_t ~", paste(n[!n %in% "V_t"], collapse = " + ")))

#parametros que podem ser modificados threshold, hidden, learnigrate

nn<- neuralnet(f, data=train_, hidden = c(15,15), threshold = 0.01,        
               stepmax = 1e+05, rep = 1, startweights = NULL, 
               learningrate.limit = NULL, 
               learningrate.factor = list(minus = 0.5, plus = 1.2), 
               learningrate=0.0115, lifesign = "full", 
               lifesign.step = 1000, algorithm = "backprop", 
               err.fct = "sse", act.fct = "logistic", 
               linear.output = F, exclude = NULL, 
               constant.weights = NULL, likelihood = FALSE) #Aqui pode-se mudar os parâmetros da rede neural

###############################################################################################################################
########################################### Vizualizar o Modelo #################################################################
###############################################################################################################################

#plot(nn)

###############################################################################################################################
########################################### Previsão do Modelo RNA ##########################################################
###############################################################################################################################

pr.nn <- compute(nn,test_[,-1])

###############################################################################################################################
########################################### "Desnormalizando" os dados de saída da RNA #############################################
###############################################################################################################################

#normalizando depois
# 
# pr.nn_ <- denormalizeData(pr.nn$net.result,getNormParameters(test_))
# test.r <- denormalizeData(test_,getNormParameters(test_))[,1]
# 
# pr.nn_ <- (pr.nn_)*(max(scaled[-index,]$V_t)-min(scaled[-index,]$V_t))+min(scaled[-index,]$V_t)
# test.r <- (test.r)*(max(scaled[-index,]$V_t)-min(scaled[-index,]$V_t))+min(scaled[-index,]$V_t)

# pr.nn_ <- pr.nn$net.result*(max(data$V_t)-min(data$V_t))+min(data$V_t)
# test.r <- (test_$V_t)*(max(data$V_t)-min(data$V_t))+min(data$V_t)

# pr.nn_ <- pr.nn$net.result*(max(data$V_t)-min(data$V_t))+min(data$V_t)
# test.r <- (test_$V_t)*(max(data$V_t)-min(data$V_t))+min(data$V_t)

# normalizando primeiro
# pr.nn_ <- pr.nn$net.result*(max(test_n[,1])-min(test_n[,1]))+min(test_n[,1])
# test.r <- (test_$V_t)*(max(test_n[,1])-min(test_n[,1]))+min(test_n[,1])
# 
# pr.nn_ <- denormalizeData(pr.nn$net.result,getNormParameters(test_n))
# test.r <- denormalizeData(test.r,getNormParameters(test_n))[,1]

pr.nn_ <- denormalizeData(pr.nn$net.result,getNormParameters(test_)) #só normalizando
test.r <- denormalizeData(test_,getNormParameters(test_))[,1] #só normalizando

# pr.nn_ <- revmapsd(pr.nn$net.result,data[,1])
# test.r <- revmapsd(test_[,1],data[,1])

###############################################################################################################################
########################################### Calculando o desempenho do modelo de RNA e do Linear #############################################
###############################################################################################################################

gofrna<- gof(sim= as.vector(t(pr.nn_)), obs = as.vector(t(test.r))) #GOF = Goodness of fit = Eficiência do ajuste
goflm<- gof(sim= as.vector(t(pr.lm)), obs = as.vector(t(test.r)))

###############################################################################################################################
########################################### Comparação dos Modelos #############################################
###############################################################################################################################

cbind(RNA=gofrna,LM=goflm)#mostrar resultados na tela

Qob= as.vector(t(test.r))

ggof(sim=as.vector(t(pr.nn_)), obs=Qob, dates = str_replace_all(paste(as.character(row.names(pr.nn_)),"/01"), fixed(" "), "")
     ,date.fmt = "%Y/%m/%d",
     main = "Volume Observado vs. Simulado",xlab = "Tempo", ylab=c("V, %"),
     gofs=c("r","PBIAS","NSE"))

ggof(sim=as.vector(t(pr.lm)), obs=as.vector(t(test.r)), dates = str_replace_all(paste(as.character(row.names(pr.nn_)),"/01"), fixed(" "), ""),
     ,date.fmt = "%Y/%m/%d",main = "Volume Observado vs. Simulado",xlab = "Tempo", ylab=c("V, %"),
     gofs=c("r","PBIAS","NSE"))

data<-cbind(date=as.Date(paste(row.names(data),"/01",sep="")),data)
datax<- data[,c(-3,-4,-6)]
datax<-melt(datax,"date")

theme_set(theme_bw())
ggplot(data = data, aes(x=as.Date(paste(row.names(data),"/01",sep="")))) + 
  geom_line(aes(y=normalizeData(V_t), col="V_t")) + 
  geom_line(aes(y=normalizeData(p_t), col="p_t")) + 
  geom_line(aes(y=normalizeData(EV), col="EV")) + 
  geom_line(aes(y=normalizeData(SS), col="SS")) +   
  labs(title="Comportamento das Séries 2008-2016", 
       subtitle="Valores Normalizados", y="Valor",x = "Tempo") +  # title and caption
  #scale_x_date(labels = lbls, breaks = brks) +  # change to monthly ticks and labels
  scale_color_manual(name="", 
                     values = c("V_t"="#a6c732", "p_t"="#ff7f00", "EV"="#019583", "SS"="#191970"), labels = c("V(t)", "P(t)","EV(t)","SS")) +  # line color
  theme(panel.grid.minor = element_blank())  # turn off minor grid

varnames<- c("V_t"="Volume (m³)", "p_t"="Precipitação (mm)", "EV"="Evaporação (mm)", "SS"="Série Sintética Monotônica")
ggplot(datax, aes(date, value, colour = variable)) + geom_line() +
  facet_wrap(~ variable, ncol = 1, scales = "free_y",labeller = as_labeller(varnames)) + 
  labs(title="Comportamento das Séries 2008-2016", subtitle="T vs. V, P, EV e SS", y="Valor",x = "Tempo (mês)") +
  scale_color_manual(name="", 
                     values = c("V_t"="#a6c732", "p_t"="#ff7f00", "EV"="#019583", "SS"="#191970"), labels = c("V(t)", "P(t)","EV(t)","SS")) +  # line color
  theme(panel.grid.minor = element_blank(),text=element_text(size=14,  family="serif"))  # turn off minor grid
  ggsave("myplot.pdf")

###############################################################################################################################
########################################### Validação Cruzada dos Modelos #############################################
###############################################################################################################################


library(boot)
set.seed(200)

# Linear model cross validation
lm.fit <- glm(V_t~.,data=data)
cv.glm(data,lm.fit,K=10)$delta[1]


# Neural net cross validation
set.seed(450)
cv.error <- NULL
k <- 10

# Initialize progress bar
library(plyr) 
pbar <- create_progress_bar('text')
pbar$init(k)

for(i in 1:k){
  index <- sample(1:nrow(data),round(0.9*nrow(data)))
  train.cv <- scaled[index,]
  test.cv <- scaled[-index,]
  
  nn <- neuralnet(f,data=train.cv,hidden=c(5,5),linear.output=F)
  
  pr.nn[[i]] <- compute(nn,test.cv[,2:6])
  
  pr.nn[[i]] <- pr.nn[[i]]$net.result*(max(data$V_t)-min(data$V_t))+min(data$V_t)
  
  test.cv.r <- (test.cv$V_t)*(max(data$V_t)-min(data$V_t))+min(data$V_t)
  
  cv.error[i] <- sum((test.cv.r - pr.nn[[i]])^2)/nrow(test.cv)
  
  pbar$step()
}

# Average MSE
mean(cv.error)

# MSE vector from CV
cv.error

# Graphical results
ggof(sim=as.vector(t(pr.nn[[4]])), obs=as.vector(t(test.cv.r)), dates = as.character(row.names(test.cv)),
     main = "Vazão Observada vs. Simulada",xlab = "Tempo", ylab=c("Q, [m3/s]"),
     gofs=c("r","PBIAS","NSE"))

# Visual plot of CV results
boxplot(cv.error,xlab='MSE CV',col='cyan',
        border='blue',names='CV error (MSE)',
        main='CV error (MSE) for NN',horizontal=TRUE)
