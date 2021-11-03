## bank data

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

#----------loading data  bank
#####
nifty_bank <- read.csv('/home/mahendra/Downloads/sem_3/TSA/project/data/bank_data.csv')
nifty_bank <- nifty_bank[852:1233,c(3,7)]
dim(nifty_bank) # 382   2
nifty_bank[,1] <- dmy(nifty_bank[,1])
tso_bank <- zoo(nifty_bank$Close, nifty_bank$Date)
######################## EDA ############################
bank <- data.frame(xts(nifty_bank$Close, order.by=as.POSIXct(nifty_bank$Date)))
names(bank) <- "bank closed"
chartSeries(bank, type = "line", show.grid = TRUE,name = "CLOSING Price of NIFTY-bank")
####################
### log Return #####
### sqr return #####
####################
Return_bank=CalculateReturns(tso_bank, method = 'log')
return_bank <- data.frame(xts(Return_bank, order.by=as.POSIXct(nifty_bank$Date)))
chartSeries(return_bank, type = "line", show.grid = TRUE,name = "Log-returns of NIFTY-bank")

rtrn_bank=Return_bank[-c(1),] # remove the first row as bank does not contain a value
chart_Series(rtrn_bank)
# histogram of the returns
chart.Histogram(return_bank,methods = c("add.density","add.normal"),colorset = c("blue","red","black")
                ,main = "histogram of the log-returns of bank data")
legend("topright",legend = c("return","kernel","normal dist"),fill = c("blue","red","black"))

sqr_Return_bank = Return_bank^2
sqr_return_bank <- data.frame(xts(sqr_Return_bank, order.by=as.POSIXct(nifty_bank$Date)))
chartSeries(sqr_return_bank, type = "line", show.grid = TRUE,name = "square of Log-returns of NIFTY-bank")
#####################################
#### Augmented Dickey Fuller Test ###
#######  ADF of returns  ############
#####################################
summary(ur.df(na.omit(Return_bank)))

#####################################
##### ACF of return #####
#### PACF of return #####
#####################################
a<- ggAcf(na.omit(as.vector(Return_bank)), col='red',main='Acf of  Log-Return of NIFTY-bank data')
p<- ggPacf(na.omit(as.vector(Return_bank)),col='steelblue',main='PAcf of  Log-Return of NIFTY-bank data')
grid.arrange(a,p, ncol = 2, nrow = 1)
############# Identifying the mean model by ARIMA #########################
arima_bank <- auto.arima(na.omit(as.vector(Return_bank)))
arima_bank
checkresiduals(arima_bank)

# adf test of the residual
summary(ur.df(resid(arima_bank),type="none",lag=1))


#  ********************************************************* NOW GARCH 

#################################################################
# Squared of Return are auto correlated.
# Squared of Return acf an pacf
##################################################################
c <- ggAcf(na.omit(as.vector(Return_bank))^2, lag.max = 40, col='green', main='ACF of squared Return Values of the bank data')
d<- ggPacf(na.omit(as.vector(Return_bank))^2,lag.max = 40, col='steelblue',main= 'PACF of squared Return Values of the bank data')
grid.arrange(c,d, ncol = 2, nrow = 1)


############################################
# Testing ARCH   ##########################
############################################
library(FinTS)
ArchTest(Return_bank,lags=1,demean = TRUE)

###############################
#### Volatility Clustering ####
###############################
sq_residual_bank <- arima_res_bank^2
ggtsdisplay(sq_residual_bank,main="Squared Residuals after fitting best ARIMA model")

chart.RollingPerformance(na.omit(Return_bank),width = 22,FUN = 'sd.annualized',scale=252, main = 'Rolling 1 month Volatility of the log-return of bank data')


############## QQ Plot ##############
ggplot(data=nifty_bank, aes(sample = as.vector(Return_bank))) +
  stat_qq() +
  stat_qq_line(col='red') + ggtitle('QQ plot of Nifty-bank Returns')

######################################
Box.test(na.omit(as.vector(Return_bank)),  lag = 1, type = "Ljung-Box", fitdf = 0)
######################################



#*************************
NIFTY_BANK_MODELS_p<-list()
NIFTY_BANK_MODELS_q<-list()
NIFTY_BANK_MODELS_P<-list()
NIFTY_BANK_MODELS_Q<-list()
NIFTY_BANK_MODELS_AIC<-list()
NIFTY_BANK_MODELS_BIC<-list()
NIFTY_BANK_MODELS_AICC<-list()

ind=0
for (p in seq(0,5)){
  for (q in seq(0,5)){
    for (P in seq(0,5)){
      for (Q in seq(0,5)){
        try({
          spec <- ugarchspec(mean.model = list(armaOrder=c(p,q)),
                             variance.model = list(model = 'eGARCH',
                                                   garchOrder = c(P,Q)),distribution = 'std')
          fit <- ugarchfit(spec = spec, data= na.omit(Return_bank)) 
          k=p+q+P+Q
          n=382
          
          AICind<-infocriteria(fit)[1]
          BICind<-infocriteria(fit)[2]
          AICcind <- AICind + (2*k*(k+1)/(n-k-1))
        })
        
        
        ind=ind+1
        NIFTY_BANK_MODELS_p[[ind]]<-p
        NIFTY_BANK_MODELS_q[[ind]]<-q
        NIFTY_BANK_MODELS_P[[ind]]<-P
        NIFTY_BANK_MODELS_Q[[ind]]<-Q
        try({
          NIFTY_BANK_MODELS_AIC[[ind]]<-AICind
          NIFTY_BANK_MODELS_BIC[[ind]]<-BICind
          NIFTY_BANK_MODELS_AICC[[ind]]<-AICcind
        })
        
        print(ind)
      }
    }
  }
}

NIFTY_BANK_MODELS<-data.frame(matrix(nrow=1296,ncol=7))#1296
columns<-c("pp","qq","PP","QQ","AIC","BIC","AICC")
colnames(NIFTY_BANK_MODELS)<-columns

NIFTY_BANK_MODELS$pp<-as.character(NIFTY_BANK_MODELS_p)
NIFTY_BANK_MODELS$qq<-as.character(NIFTY_BANK_MODELS_q)
NIFTY_BANK_MODELS$PP<-as.character(NIFTY_BANK_MODELS_P)
NIFTY_BANK_MODELS$QQ<-as.character(NIFTY_BANK_MODELS_Q)
NIFTY_BANK_MODELS$AIC<-as.character(NIFTY_BANK_MODELS_AIC)
NIFTY_BANK_MODELS$BIC<-as.character(NIFTY_BANK_MODELS_BIC)
NIFTY_BANK_MODELS$AICC<-as.character(NIFTY_BANK_MODELS_AICC)
#View(NIFTY_BANK_MODELS)
write.csv(NIFTY_BANK_MODELS,file = "BANK_score.csv",sep=",")
#*************************
#*

#************************************

dat<-read.csv('/home/mahendra/Downloads/sem_3/TSA/project/data/BANK_score.csv')
df<-dat%>%select(X,AIC,AICC,BIC)%>%filter(AIC<0)%>%filter(BIC<0)%>%filter(AICC<0)
d <- melt(df, id.vars="X")
d

ggplot(data=d,
       aes(x=X, y=value, colour=variable)) +
  geom_line()+ labs(x="sl no of different combination of ARIMA and GARCH model",
                    y="score",title = "AIC,BIC and AICc score of different model")



###*************** Best model specification and fitting 
garch_bank <- ugarchspec(mean.model = list(armaOrder=c(0,0)),
                         variance.model = list(model = 'eGARCH', 
                                               garchOrder = c(1, 1)),distribution = 'std')
fit_garch_bank <- ugarchfit(spec = garch_bank, data= na.omit(as.vector(Return_bank)))
## forecasting
forecast_bank<- ugarchforecast(fit_garch_bank,n.ahead = 30)

par(mfrow=c(1,2))
plot(forecast_bank,which=1)
plot(forecast_bank,which=3)

########### going back to original data 


nifty_Bank <- read.csv('/home/mahendra/Downloads/sem_3/TSA/project/data/Bank_data.csv')
nifty_Bank <- nifty_Bank[852:1251,c(3,7)]
nifty_Bank[,1] <- dmy(nifty_Bank[,1])
original_Bank <- nifty_Bank$Close
Update <- c()
end=original_Bank[382]
for (i in seq(1,18)){
  end= end*exp(forecast_bank@forecast$seriesFor[i])
  print(end)
  Update <- c(Update,end)
}
par(mfrow=c(1,1))
plot(c(1:382),original_Bank[1:382],type="l",col="black",xlim=c(1,420),ylim=c(10000,45000),main="Forcasting the original stock value",xlab="time point",ylab="close price of nifty-it ",xaxt='n')
lines(c(383:400),original_Bank[383:400],type="l",col="green")
lines(c(383:400),Update,type="p",col="red")
legend("bottomright",legend = c("forecasted stock values","original privious values","original future ground truths"),
       fill = c("red","black","green"))

## RMSE
nifty_Bank<- read.csv('/home/mahendra/Downloads/sem_3/TSA/project/data/Bank_data.csv')
nifty_Bank <- nifty_Bank[852:1251,c(3,7)]
nifty_Bank[,1] <- dmy(nifty_Bank[,1])
tso_Bank <- zoo(nifty_Bank$Close, nifty_Bank$Date)
Return_Bank=CalculateReturns(tso_Bank, method = 'log')
true_returns_bank <- na.omit(as.vector(Return_Bank))

predicted_returns_bank <- c()
predicted_stocks_bank <- c()
total_sqr_loss_in_returns_bank <- 0
total_sqr_loss_in_stock_bank <- 0
for (i in seq(1,18)){
  fit_garch_bank <- ugarchfit(spec = garch_bank, data = true_returns_bank[1:(381-1+i)] )
  forecast_bank<- ugarchforecast(fit_garch_bank,n.ahead = 1 )
  pred_return=forecast_bank@forecast$seriesFor[1]
  predicted_returns_bank= c(predicted_returns_bank,pred_return)
  sqr_loss_return= (pred_return - true_returns_bank[381+i])^2
  total_sqr_loss_in_returns_bank = total_sqr_loss_in_returns_bank + sqr_loss_return
  
  
  previous_stock = nifty_Bank$Close[382-1+i]
  pred_stock = previous_stock*exp(pred_return)
  predicted_stocks_bank <- c( predicted_stocks_bank, pred_stock)
  print(pred_stock)
  sqr_loss_stock= (pred_stock - nifty_Bank$Close[382+i])^2
  total_sqr_loss_in_stock_bank = total_sqr_loss_in_stock_bank + sqr_loss_stock
}

predicted_returns_bank
predicted_stocks_bank
(total_sqr_loss_in_returns_bank^0.5)/length(predicted_returns_bank)
(total_sqr_loss_in_stock_bank^0.5)/length(predicted_stocks_bank)

plot(c(1:382),original_Bank[1:382],type="l",col="black",xlim=c(1,420),ylim=c(10000,45000),main="Forcasting the original stock value",xlab="time point",ylab="close price of nifty-it ",xaxt='n')
lines(c(383:400),original_Bank[383:400],type="l",col="green")
lines(c(383:400),predicted_stocks_bank,type="l",col="red")
legend("bottomright",legend = c("forecasted stock values","original privious values","original future ground truths"),
       fill = c("red","black","green"))

## zoom in
plot(c(1:382),original[1:382],type="l",col="black",xlim=c(370,420),ylim=c(17000,45000),main="Forcasting the original stock value of bank",xlab="time point",ylab="close price of nifty-it ",xaxt='n')
lines(c(383:400),original[383:400],type="b",col="green")
lines(c(383:400),predicted_stocks_bank,type="b",col="red")
legend("bottomright",legend = c("forecasted stock values","original privious values","original future ground truths"),
       fill = c("red","black","green"))



