## IT data

library(PerformanceAnalytics)
library(astsa)
library(itsmr)
library(lubridate)
library(zoo)
library(randtests)
library(forecast)
library(urca)
library(aTSA)
library(ggplot2)
library(tsoutliers)
library(gridExtra)
library(rugarch)
library(tseries)
library(quantmod)
library(dplyr)
library(reshape2)

#----------loading data  IT
#####
nifty_it <- read.csv('/home/mahendra/Downloads/sem_3/TSA/project/data/it_data.csv')
nifty_it <- nifty_it[852:1233,c(3,7)]
#nifty_it <- nifty_it[1222:1233,c(3,7)]
#View(nifty_it)
dim(nifty_it) # 382   2
nifty_it[,1] <- dmy(nifty_it[,1])
#plot(nifty_it$Close, ylab="Stock Prices",main="Figure : Closing prices of the stocks",type = 'l')
tso_it <- zoo(nifty_it$Close, nifty_it$Date)
#plot(tso_it,type = 'l')
######################## EDA ############################
it <- data.frame(xts(nifty_it$Close, order.by=as.POSIXct(nifty_it$Date)))
names(it) <- "it closed"
chartSeries(it, type = "line", show.grid = TRUE,name = "CLOSING Price of NIFTY-IT")
########acf and pacf of the original data just to see
#m<- ggAcf(nifty_it$Close, col='red',main='Acf of NIFTY-IT original stock price')
#n<- ggPacf(nifty_it$Close,col='steelblue',main='PAcf of NIFTY-IT original stock price')
#grid.arrange(m,n, ncol = 2, nrow = 1)
########
####################
### log Return #####
### sqr return #####
####################
Return_it=CalculateReturns(tso_it, method = 'log')
return_it <- data.frame(xts(Return_it, order.by=as.POSIXct(nifty_it$Date)))
chartSeries(return_it, type = "line", show.grid = TRUE,name = "Log-returns of NIFTY-IT")

rtrn_it=Return_it[-c(1),] # remove the first row as it does not contain a value
chart_Series(rtrn_it)
          # histogram of the returns
chart.Histogram(return_it,methods = c("add.density","add.normal"),
                colorset = c("blue","red","black"),
                main = "histogram of the log-returns of Nifty-IT data")
legend("topright",legend = c("return","kernel","normal dist"),fill = c("blue","red","black"))
#plot.ts(Return_it,type="o",xlab="Date",ylab="log return ", main="Log return of IT")
#sum(na.omit(Return_it))/length(na.omit(Return_it)) # mean= 0.003004869

sqr_Return_it = Return_it^2
sqr_return_it <- data.frame(xts(sqr_Return_it, order.by=as.POSIXct(nifty_it$Date)))
chartSeries(sqr_return_it, type = "line", show.grid = TRUE,name = "square of Log-returns of NIFTY-IT")
#plot.ts(sqr_Return_it,type="o",xlab="Date",ylab="square log return ", main="Squared Log return of IT")
#####################################
#### Augmented Dickey Fuller Test ###
#######  ADF of returns  ############
#####################################
#summary(ur.df(logret,type='drift'))
summary(ur.df(na.omit(Return_it)))

#####################################
##### ACF of return #####
#### PACF of return #####
#####################################
a<- ggAcf(na.omit(as.vector(Return_it)), col='red',main='Acf of  Log-Return of NIFTY-IT data')
p<- ggPacf(na.omit(as.vector(Return_it)),col='steelblue',main='PAcf of  Log-Return of NIFTY-IT data')
grid.arrange(a,p, ncol = 2, nrow = 1)
############# Identifying the mean model by ARIMA #########################
arima_it <- auto.arima(na.omit(as.vector(Return_it)))
arima_it
checkresiduals(arima_it)

# adf test of the residual
summary(ur.df(resid(arima_it),type="none",lag=1))
#autoplot(arima_it)


#  ********************************************************* NOW GARCH 

#################################################################
#Absolute Return or Squared of Return are auto correlated.
#Absolute Return or Squared of Return acf an pacf
##################################################################
#a<- ggAcf(abs(na.omit(as.vector(Return_it))), col='red',main='Acf of Absolute Return_it of NIFTY')
#p<- ggPacf(abs(na.omit(as.vector(Return_it))),col='steelblue',main='PAcf of Absolute Return_it of NIFTY')
#grid.arrange(a,p, ncol = 2, nrow = 1)


c <- ggAcf(na.omit(as.vector(Return_it))^2, lag.max = 40, col='red', main='ACF of squared of log-Return Values of the IT data')
d<- ggPacf(na.omit(as.vector(Return_it))^2,lag.max = 40, col='steelblue',main= 'PACF of squared of log-Return Values of the IT data')
grid.arrange(c,d, ncol = 2, nrow = 1)


############################################
# Testing ARCH   ##########################
############################################
library(FinTS)
ArchTest(Return_it,lags=1,demean = TRUE)

###############################
#### Volatility Clustering ####
###############################
#arima_res_it <- arima_it$residuals
#ggtsdisplay(arima_res_it,main="Residuals after fitting best-ARIMA model")

sq_residual_it <- arima_res_it^2
ggtsdisplay(sq_residual_it,main="Squared Residuals after fitting best-ARIMA model")

chart.RollingPerformance(na.omit(Return_it),width = 22,FUN = 'sd.annualized',scale=252, main = 'Rolling 1 month Volatility of the log-return of It data')

######Skewness Kurtois ##############
ggplot(aes(as.vector(na.omit(Return_it))), data=na.omit(Return_it)) + 
  geom_histogram(bins = 100,col='black',fill='red') + 
  ggtitle('Return_it of MSFt')
skewness=skewness((as.vector(na.omit(Return_it))))
kurtosis=kurtosis((as.vector(na.omit(Return_it))))
sprintf("skewness= %f kurtosis= %f",skewness,kurtosis)
############## QQ Plot ##############
ggplot(data=nifty_it, aes(sample = as.vector(Return_it))) +
  stat_qq() +
  stat_qq_line(col='red') + ggtitle('QQ plot of Nifty-IT Returns')

######################################
Box.test(na.omit(as.vector(Return_it)),  lag = 1, type = "Ljung-Box", fitdf = 0)
######################################


##################################################################################
################################# GARCH Model ####################################
##################################################################################


NIFTY_IT_MODELS_p<-list()
NIFTY_IT_MODELS_q<-list()
NIFTY_IT_MODELS_P<-list()
NIFTY_IT_MODELS_Q<-list()
NIFTY_IT_MODELS_AIC<-list()
NIFTY_IT_MODELS_BIC<-list()
NIFTY_IT_MODELS_AICC<-list()

ind=0
for (p in seq(0,5)){
  for (q in seq(0,5)){
    for (P in seq(0,5)){
      for (Q in seq(0,5)){
        try({
          spec <- ugarchspec(mean.model = list(armaOrder=c(p,q)),
                             variance.model = list(model = 'eGARCH',
                                                   garchOrder = c(P,Q)),distribution = 'std')
          fit <- ugarchfit(spec = spec, data= na.omit(Return_it)) 
          k=p+q+P+Q
          n=382
          
          AICind<-infocriteria(fit)[1]
          BICind<-infocriteria(fit)[2] })
        AICcind <- AICind + (2*k*(k+1)/(n-k-1))
        
        ind=ind+1
        NIFTY_IT_MODELS_p[[ind]]<-p
        NIFTY_IT_MODELS_q[[ind]]<-q
        NIFTY_IT_MODELS_P[[ind]]<-P
        NIFTY_IT_MODELS_Q[[ind]]<-Q
        try({
          NIFTY_IT_MODELS_AIC[[ind]]<-AICind
          NIFTY_IT_MODELS_BIC[[ind]]<-BICind
          NIFTY_IT_MODELS_AICC[[ind]]<-AICcind
        })
        
        print(ind)
        #it_aic[i,j] <- infocriteria(r)[1]
        #it_bic <- infocriteria(r)[2]
        #NIFTY_IT_MODELS[nrow(NIFTY_IT_MODELS)+1,"pp"]<-p
        #NIFTY_IT_MODELS[nrow(NIFTY_IT_MODELS)+1,"qq"]<-q
        #NIFTY_IT_MODELS[nrow(NIFTY_IT_MODELS)+1,"PP"]<-P
        #NIFTY_IT_MODELS[nrow(NIFTY_IT_MODELS)+1,"QQ"]<-Q
        #NIFTY_IT_MODELS[nrow(NIFTY_IT_MODELS)+1,"AIC"]<-infocriteria(r)[1]
        #NIFTY_IT_MODELS[nrow(NIFTY_IT_MODELS)+1,"BIC"]<-infocriteria(r)[2]"
        
        
        #print(summary(k))
        #print("--------------------------------")
        #NIF_AIC[i,j]<- stats::AIC(k)
        #plot.ts(k$fitted.values[-1,1]**2,main=paste("Estimated GARCH(",paste(i,j,sep=","),") variance for the 'Nifty' dataset"))
      }
    }
  }
}
NIFTY_IT_MODELS<-data.frame(matrix(nrow=1296,ncol=7))#1296
columns<-c("pp","qq","PP","QQ","AIC","BIC","AICC")
colnames(NIFTY_IT_MODELS)<-columns

NIFTY_IT_MODELS$pp<-as.character(NIFTY_IT_MODELS_p)
NIFTY_IT_MODELS$qq<-as.character(NIFTY_IT_MODELS_q)
NIFTY_IT_MODELS$PP<-as.character(NIFTY_IT_MODELS_P)
NIFTY_IT_MODELS$QQ<-as.character(NIFTY_IT_MODELS_Q)
NIFTY_IT_MODELS$AIC<-as.character(NIFTY_IT_MODELS_AIC)
NIFTY_IT_MODELS$BIC<-as.character(NIFTY_IT_MODELS_BIC)
NIFTY_IT_MODELS$AICC<-as.character(NIFTY_IT_MODELS_AICC)
#View(NIFTY_IT_MODELS)
write.csv(NIFTY_IT_MODELS,file = "IT_score.csv",sep=",")

#************************************

dat<-read.csv('/home/mahendra/Downloads/sem_3/TSA/project/data/IT_score.csv')
df<-dat%>%select(X,AIC,AICC,BIC)%>%filter(AIC<0)%>%filter(BIC<0)%>%filter(AICC<0)
d <- melt(df, id.vars="X")

ggplot(data=d,
       aes(x=X, y=value, colour=variable)) +
  geom_line()+ labs(x="sl no of different combination of ARIMA and GARCH model", y="score",title = "AIC,BIC and AICc score of different model")



###*************** Best model specification and fitting 
garch_it <- ugarchspec(mean.model = list(armaOrder=c(4,2)),
                             variance.model = list(model = 'eGARCH', 
                                                   garchOrder = c(4,3)),distribution = 'std')
fit_garch_it <- ugarchfit(spec = garch_it, data= na.omit(as.vector(Return_it)))
#fit_garch_it <- ugarchfit(garch_it, data= na.omit(Return_it))
#fit_garch_it
#plot(fit_garch_it,which='all')
## forecasting
forecast_it<- ugarchforecast(fit_garch_it,n.ahead = 30)
#forecast_it
#forecast_it@forecast$seriesFor

par(mfrow=c(1,2))
plot(forecast_it,which=1)
plot(forecast_it,which=3)

########### going back to original data 


nifty_It <- read.csv('/home/mahendra/Downloads/sem_3/TSA/project/data/It_data.csv')
nifty_It <- nifty_It[852:1251,c(3,7)]
nifty_It[,1] <- dmy(nifty_It[,1])
original_It <- nifty_It$Close
Update <- c()
end=original_It[382]
for (i in seq(1,18)){
  end= end*exp(forecast_it@forecast$seriesFor[i])
  print(end)
  Update <- c(Update,end)
}
par(mfrow=c(1,1))
plot(c(1:382),original_It[1:382],type="l",col="black",xlim=c(1,420),
     ylim=c(10000,45000),main="Forcasting the original stock value",
     xlab="time point",ylab="close price of nifty-it ",xaxt='n')
lines(c(383:400),original_It[383:400],type="l",col="green")
lines(c(383:400),Update,type="p",col="red")
legend("bottomright",legend = c("forecasted stock values","original privious values",
                                "original future ground truths"),
       fill = c("red","black","green"))

## RMSE
nifty_It<- read.csv('/home/mahendra/Downloads/sem_3/TSA/project/data/It_data.csv')
nifty_It <- nifty_It[852:1251,c(3,7)]
nifty_It[,1] <- dmy(nifty_It[,1])
tso_It <- zoo(nifty_It$Close, nifty_It$Date)
Return_It=CalculateReturns(tso_It, method = 'log')
true_returns_it <- na.omit(as.vector(Return_It))

predicted_returns_it <- c()
predicted_stocks_it <- c()
total_sqr_loss_in_returns_it <- 0
total_sqr_loss_in_stock_it <- 0
for (i in seq(1,18)){
  fit_garch_it <- ugarchfit(spec = garch_it, data = true_returns_it[1:(381-1+i)] )
  forecast_it<- ugarchforecast(fit_garch_it,n.ahead = 1 )
  pred_return=forecast_it@forecast$seriesFor[1]
  predicted_returns_it= c(predicted_returns_it,pred_return)
  sqr_loss_return= (pred_return - true_returns_it[381+i])^2
  total_sqr_loss_in_returns_it = total_sqr_loss_in_returns_it + sqr_loss_return
  
  
  previous_stock = nifty_It$Close[382-1+i]
  pred_stock = previous_stock*exp(pred_return)
  predicted_stocks_it <- c( predicted_stocks_it, pred_stock)
  print(pred_stock)
  sqr_loss_stock= (pred_stock - nifty_It$Close[382+i])^2
  total_sqr_loss_in_stock_it = total_sqr_loss_in_stock_it + sqr_loss_stock
}

predicted_returns_it
predicted_stocks_it
(total_sqr_loss_in_returns_it^0.5)/length(predicted_returns_it)
(total_sqr_loss_in_stock_it^0.5)/length(predicted_stocks_it)

plot(c(1:382),original_It[1:382],type="l",col="black",xlim=c(1,420),ylim=c(10000,45000),
     main="Forcasting the original stock value",xlab="time point",
     ylab="close price of nifty-it ",xaxt='n')
lines(c(383:400),original_It[383:400],type="l",col="green")
lines(c(383:400),predicted_stocks_it,type="l",col="red")
legend("bottomright",legend = c("forecasted stock values","original privious values",
                                "original future ground truths"),
       fill = c("red","black","green"))

## zoom in
plot(c(1:382),original[1:382],type="l",col="black",xlim=c(370,420),ylim=c(17000,45000),
     main="Forcasting the original stock value of it",xlab="time point",
     ylab="close price of nifty-it ",xaxt='n')
lines(c(383:400),original[383:400],type="b",col="green")
lines(c(383:400),predicted_stocks_it,type="b",col="red")
legend("bottomright",legend = c("forecasted stock values","original privious values",
                                "original future ground truths"),
       fill = c("red","black","green"))


#################

#data<-data.frame(NIFTY-IT=as.numeric(nifty_it$Close),NIFTY-BANK=as.numeric(nifty_bank$Close),NIFTY-OIL=as.numeric(nifty_oil$Close),NIFTY-METAL=as.numeric(nifty_metal$Close))

data<- cbind(nifty_it$Close,nifty_bank$Close,nifty_oil$Close,nifty_metal$Close)

colnames(data)<-c("NIFTY-IT","NIFTY-BANK","NIFTY-OIL","NIFTY-METAL")
correl<-cor(data)
library(corrplot)
corrplot(correl, type = "upper", method="square", order = "hclust",
         tl.col = "black", tl.srt = 30,addCoef.col = "white")


