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

list.of.packages <- c("hydroTSM", "lubridate",'openxlsx','reshape2','plyr','neuralnet',"stringr","hydroGOF","pracma","RSNNS") #list pacotes
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])] #instala se não houver
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

###############################################################################################################################
########################################### Definindo diretório de trabalho(pasta) e carregando funções #######################
###############################################################################################################################

Sys.setenv(R_ZIPCMD= "C:/Rtools/bin/zip")
setwd("~/Universidade/IFPB/Pesquisa/Coordenador de Projetos Pesquisa/Pré-tratamento dos dados") ## mudar a o local aonde estão os arquivos
path<- getwd();
source("lagpadf.r")

###############################################################################################################################
########################################### Carregando dados externos e definindo dados de entrada para a RNA #################
###############################################################################################################################

load("~/Universidade/IFPB/Pesquisa/Coordenador de Projetos Pesquisa/Pré-tratamento dos dados/Dados_Zoo.RData") ## Carregar dados gerados pelo hydroweb

pi_plu_ser_1911_2016<- rbind.zoo(a1911_1991,a1994_2016) #junta as duas séries
pi_plu_ser_1911_2016_full <- merge(pi_plu_ser_1911_2016,zoo(,seq(start(pi_plu_ser_1911_2016),
                                                                 end(pi_plu_ser_1911_2016),by="day")), all=TRUE)
wb <- loadWorkbook("Vol_Jatoba_II.xlsx") # colocar o nome do arquivo do excel
pi_mensal<- daily2monthly(pi_plu_ser_1911_2016_full, FUN=sum, na.rm=TRUE)
vol_jatoba = read.xlsx(wb, sheet = 1,detectDates = TRUE)  
jatoba_ser<- zoo(as.numeric(gsub(",",".",vol_jatoba$`Volume.(%)`)),format(as.Date(vol_jatoba$Data.do.registro,"%d/%m/%Y"),"%Y/%m"))
posto_plu<- window(pi_mensal,index. = index(pi_mensal), as.Date("2008-01-01"),as.Date("2016-12-01"))
posto_plu_2<- zoo(coredata(posto_plu),format(as.Date(time(posto_plu),"%d/%m/%Y"),"%Y/%m"))
dataset<- merge(posto_plu_2,jatoba_ser,all = FALSE)
dataset<- as.data.frame(dataset)
data<- dataset
rown<- rownames(data)

picos<- data[,2]
vales<- data[,2]

seqds<- cumsum(c(1, sign(diff(data[,2]))))

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

picos[picos>0]<- 3
vales[vales>0]<- -54
# data<- cbind(data,data$jatoba_ser,Picos=picos[dropout=T])
data<- cbind(data,Com=seqds[dropout=T],Picos=picos[dropout=T],Vales=vales[dropout=T])


# data<- cbind(data[,"jatoba_ser"],data[,"jatoba_ser"],data[,"jatoba_ser"],
#              data[,"posto_plu_2"],data[,"posto_plu_2"],data[,"posto_plu_2"],data[,"Com"],data[,"Picos"],data[,"Vales"])

# data<- cbind(data[,"jatoba_ser"],data[,"jatoba_ser"],data[,"jatoba_ser"],
#              data[,"posto_plu_2"],data[,"posto_plu_2"],data[,"Com"],data[,"Picos"],data[,"Vales"])

data<- cbind(data[,"jatoba_ser"],data[,"jatoba_ser"],data[,"jatoba_ser"],
             data[,"posto_plu_2"],data[,"posto_plu_2"],data[,"Com"])

# colnames(data)<- c("V_t","V_t_1","V_t_2","p_t","p_t_1","p_t_2","Com","Picos","Vales")
# colnames(data)<- c("V_t","V_t_1","V_t_2","p_t","p_t_2","Com","Picos","Vales")
colnames(data)<- c("V_t","V_t_1","V_t_2","p_t","p_t_2","Com")
coln<- colnames(data)


###############################################################################################################################
########################################### Definindo os atrasos (lag) para a RNA #############################################
###############################################################################################################################

# data<- cbind(data[,"V_t"],lagpad(data[,"V_t_1"],1),lagpad(data[,"V_t_2"],2),data[,"p_t"]
#              ,lagpad(data[,"p_t_2"],2),data[,"Com"],data[,"Picos"],data[,"Vales"])

data<- cbind(data[,"V_t"],lagpad(data[,"V_t_1"],1),lagpad(data[,"V_t_2"],2),data[,"p_t"]
             ,lagpad(data[,"p_t_2"],2),data[,"Com"])

colnames(data)<-coln
rownames(data)<- rown
data<- as.data.frame(data)

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

maxs <- apply(data, 2, max)
mins <- apply(data, 2, min)
scaled <- as.data.frame(scale(data, center = mins, scale = maxs - mins))
# scaled <- as.data.frame(scale(data))
# scaled<- data
# colnames(scaled)<- coln

###############################################################################################################################
########################################### Definindo os dados de treino e teste escalonados #############################################
###############################################################################################################################

# train_ <- normalizeData(scaled[index,])
# test_ <- normalizeData(scaled[-index,])
train_ <- scaled[1:73,]
test_ <- scaled[75:98,]
# colnames(train_)<- coln
# colnames(test_)<- coln

###############################################################################################################################
########################################### Treinando a Rede Neural Artificial #############################################
###############################################################################################################################


n <- names(train_)
f <- as.formula(paste("V_t ~", paste(n[!n %in% "V_t"], collapse = " + ")))

nn<- neuralnet(f, data=train_, hidden = c(10,10), threshold = 0.01,        
               stepmax = 1e+05, rep = 1, startweights = NULL, 
               learningrate.limit = NULL, 
               learningrate.factor = list(minus = 0.5, plus = 1.2), 
               learningrate=0.0015, lifesign = "full", 
               lifesign.step = 1000, algorithm = "backprop", 
               err.fct = "sse", act.fct = "logistic", 
               linear.output = F, exclude = NULL, 
               constant.weights = NULL, likelihood = FALSE) #Aqui pode-se mudar os parâmetros da rede neural

###############################################################################################################################
########################################### Vizualizar o Modelo #################################################################
###############################################################################################################################

# plot(nn)

###############################################################################################################################
########################################### Previsão do Modelo RNA ##########################################################
###############################################################################################################################

pr.nn <- compute(nn,test_[,-1])

###############################################################################################################################
########################################### "Desnormalizando" os dados de saída da RNA #############################################
###############################################################################################################################

pr.nn_ <- pr.nn$net.result*(max(data$V_t)-min(data$V_t))+min(data$V_t)
test.r <- (test_$V_t)*(max(data$V_t)-min(data$V_t))+min(data$V_t)

# pr.nn_ <- denormalizeData(pr.nn$net.result,getNormParameters(test_))
# test.r <- denormalizeData(test_,getNormParameters(test_))[,1]

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

cbind(RNA=gofrna,LM=goflm)

ggof(sim=as.vector(t(pr.nn_)), obs=as.vector(t(test.r)), dates = str_replace_all(paste(as.character(row.names(pr.nn_)),"/01"), fixed(" "), "")
     ,date.fmt = "%Y/%m/%d",
     main = "Vazão Observada vs. Simulada",xlab = "Tempo", ylab=c("Q, [m3/s]"),
     gofs=c("r","PBIAS","NSE","RMSE"))

ggof(sim=as.vector(t(pr.lm)), obs=as.vector(t(test.r)), dates = str_replace_all(paste(as.character(row.names(pr.nn_)),"/01"), fixed(" "), ""),
     ,date.fmt = "%Y/%m/%d",main = "Vazão Observada vs. Simulada",xlab = "Tempo", ylab=c("Q, [m3/s]"),
     gofs=c("r","PBIAS","NSE","RMSE"))

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
