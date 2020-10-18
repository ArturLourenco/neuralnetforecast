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

list.of.packages <- c("hydroTSM", "lubridate",'openxlsx')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
require(hydroTSM)
require(lubridate)
require(openxlsx)
require(reshape2)
library(lattice)
#library(chron)
require(quantmod)
require(ggplot2)
require(reshape2)
library(plyr)
require(scales)
library(colorRamps)

Sys.setenv(R_ZIPCMD= "C:/Rtools/bin/zip")
setwd("G:/Meu Drive/Pesquisa/Coordenador de Projetos Pesquisa/Pré-tratamento dos dados") ## mudar a o local aonde estão os arquivos
load("G:/Meu Drive/Pesquisa/Coordenador de Projetos Pesquisa/Pré-tratamento dos dados/Dados_Zoo.RData") ## Carregar dados gerados pelo hydroweb
source("calendarplot.R")
path<- getwd();

###############################################################################################################################
########################################### Análise e pré tratamento dos dados ################################################
###############################################################################################################################

pi_plu_ser_1911_2016<- rbind.zoo(a1911_1991,a1994_2016)
pi_plu_ser_1911_2016_full <- merge(pi_plu_ser_1911_2016,zoo(,seq(start(pi_plu_ser_1911_2016),end(pi_plu_ser_1911_2016),by="day")), all=TRUE)

# Make a dataframe
df<-data.frame(date=date(pi_plu_ser_1911_2016_full),coredata(pi_plu_ser_1911_2016_full))

df<- df[35065:35794,]
# df<- df[!(weekdays(as.Date(df$date)) %in% c('sábado','domingo')),]

# We will facet by year ~ month, and each subgraph will
# show week-of-month versus weekday
# the year is simple
df$year<-as.numeric(as.POSIXlt(df$date)$year+1900)
# the month too 
df$month<-as.numeric(as.POSIXlt(df$date)$mon+1)
# but turn months into ordered facors to control the appearance/ordering in the presentation
df$monthf<-factor(df$month,levels=as.character(1:12),labels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),ordered=TRUE)
# the day of week is again easily found
df$weekday = as.POSIXlt(df$date)$wday
# again turn into factors to control appearance/abbreviation and ordering
# I use the reverse function rev here to order the week top down in the graph
# you can cut it out to reverse week order
df$weekdayf<-factor(df$weekday,levels=rev(1:7),labels=rev(c("Mon","Tue","Wed","Thu","Fri","Sat","Sun")),ordered=TRUE)
# the monthweek part is a bit trickier 
# first a factor which cuts the dfa into month chunks
df$yearmonth<-as.yearmon(df$date)
df$yearmonthf<-factor(df$yearmonth)
# then find the "week of year" for each day
df$week <- as.numeric(format(df$date,"%W"))
# and now for each monthblock we normalize the week to start at 1 
df<-ddply(df,.(yearmonthf),transform,monthweek=1+week-min(week))

# Now for the plot
P<- ggplot(df, aes(monthweek, weekdayf, fill = coredata.pi_plu_ser_1911_2016_full.)) + 
  geom_tile(colour = "white") + facet_grid(year~monthf) + scale_fill_gradient(low="red", high="blue") +
  labs(title = "Time-Series Calendar Heatmap") +  xlab("Week of Month") + ylab("")
P

P2<- ggplot(df, aes(monthweek, weekdayf, fill = coredata.pi_plu_ser_1911_2016_full.)) + 
  geom_tile(colour = "white") + facet_grid(year~monthf) + scale_color_gradientn(colours = matlab.like2(5),name = 'Cº') +
  labs(title = "Time-Series Calendar Heatmap") +  xlab("Week of Month") + ylab("")
P2

# Plot as calendar heatmap
calendarHeat(df$date, df$coredata.pi_plu_ser_1911_2016_full., 
             varname="Pricipitação Princesa Isabel")
