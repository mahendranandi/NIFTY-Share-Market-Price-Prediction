## metal data


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

#----------loading data  metal
#####
nifty_metal <- read.csv('/home/mahendra/Downloads/sem_3/TSA/project/data/metal_data.csv')
nifty_metal <- nifty_metal[852:1233,c(3,7)]
dim(nifty_metal) # 382   2
nifty_metal[,1] <- dmy(nifty_metal[,1])
tso_metal <- zoo(nifty_metal$Close, nifty_metal$Date)
######################## EDA ############################
metal <- data.frame(xts(nifty_metal$Close, order.by=as.POSIXct(nifty_metal$Date)))
names(metal) <- "metal closed"
chartSeries(metal, type = "line", show.grid = TRUE,name = "CLOSING Price of NIFTY-metal")
####################
### log Return #####
### sqr return #####
####################
Return_metal=CalculateReturns(tso_metal, method = 'log')
return_metal <- data.frame(xts(Return_metal, order.by=as.POSIXct(nifty_metal$Date)))
chartSeries(return_metal, type = "line", show.grid = TRUE,name = "Log-returns of NIFTY-metal")

rtrn_metal=Return_metal[-c(1),] # remove the first row as metal does not contain a value
chart_Series(rtrn_metal)
# histogram of the returns
chart.Histogram(return_metal,methods = c("add.density","add.normal"),colorset = c("blue","red","black"),
                main = "histogram of log-returns of the NIFTY_OIL data")
legend("topright",legend = c("return","kernel","normal dist"),fill = c("blue","red","black"))

sqr_Return_metal = Return_metal^2
sqr_return_metal <- data.frame(xts(sqr_Return_metal, order.by=as.POSIXct(nifty_metal$Date)))
chartSeries(sqr_return_metal, type = "line", show.grid = TRUE,name = "square of Log-returns of NIFTY-metal")
#####################################
#### Augmented Dickey Fuller Test ###
#######  ADF of returns  ############
#####################################
summary(ur.df(na.omit(Return_metal)))

#####################################
##### ACF of return #####
#### PACF of return #####
#####################################
a<- ggAcf(na.omit(as.vector(Return_metal)), col='blue',main='Acf of  Log-Return of NIFTY-metal data')
p<- ggPacf(na.omit(as.vector(Return_metal)),col='green',main='PAcf of  Log-Return of NIFTY-metal data')
grid.arrange(a,p, ncol = 2, nrow = 1)
############# Identifying the mean model by ARIMA #########################
arima_metal <- auto.arima(na.omit(as.vector(Return_metal)))
arima_metal
checkresiduals(arima_metal)

# adf test of the residual
summary(ur.df(resid(arima_metal),type="none",lag=1))


#  ********************************************************* NOW GARCH 

#################################################################
# Squared of Return are auto correlated.
# Squared of Return acf an pacf
##################################################################
c <- ggAcf(na.omit(as.vector(Return_metal))^2, lag.max = 40, col='red', main='ACF of squared Return Values of the metal data')
d<- ggPacf(na.omit(as.vector(Return_metal))^2,lag.max = 40, col='steelblue',main= 'PACF of squared Return Values of the metal data')
grid.arrange(c,d, ncol = 2, nrow = 1)


############################################
# Testing ARCH   ##########################
############################################
library(FinTS)
ArchTest(Return_metal,lags=1,demean = TRUE)

###############################
#### Volatility Clustering ####
###############################
sq_residual_metal <- arima_res_metal^2
ggtsdisplay(sq_residual_metal,main="Squared Residuals after fitting best ARIMA model")

chart.RollingPerformance(na.omit(Return_metal),width = 22,FUN = 'sd.annualized',scale=252, main = 'Rolling 1 month Volatility of the log-return of metal data')


############## QQ Plot ##############
ggplot(data=nifty_metal, aes(sample = as.vector(Return_metal))) +
  stat_qq() +
  stat_qq_line(col='red') + ggtitle('QQ plot of Nifty Returns')

######################################
Box.test(na.omit(as.vector(Return_metal)),  lag = 1, type = "Ljung-Box", fitdf = 0)
######################################



##################################################################################
################################# GARCH Model ####################################
##################################################################################
#************************************
NIFTY_METAL_MODELS_p<-list()
NIFTY_METAL_MODELS_q<-list()
NIFTY_METAL_MODELS_P<-list()
NIFTY_METAL_MODELS_Q<-list()
NIFTY_METAL_MODELS_AIC<-list()
NIFTY_METAL_MODELS_BIC<-list()
NIFTY_METAL_MODELS_AICC<-list()

ind=0
for (p in seq(0,5)){
  for (q in seq(0,5)){
    for (P in seq(0,5)){
      for (Q in seq(0,5)){
        try({
          spec <- ugarchspec(mean.model = list(armaOrder=c(p,q)),
                             variance.model = list(model = 'eGARCH',
                                                   garchOrder = c(P,Q)),distribution = 'std')
          fit <- ugarchfit(spec = spec, data= na.omit(Return_metal)) 
          k=p+q+P+Q
          n=382
          
          AICind<-infocriteria(fit)[1]
          BICind<-infocriteria(fit)[2]
          AICcind <- AICind + (2*k*(k+1)/(n-k-1))
        })
        
        
        ind=ind+1
        NIFTY_METAL_MODELS_p[[ind]]<-p
        NIFTY_METAL_MODELS_q[[ind]]<-q
        NIFTY_METAL_MODELS_P[[ind]]<-P
        NIFTY_METAL_MODELS_Q[[ind]]<-Q
        try({
          NIFTY_METAL_MODELS_AIC[[ind]]<-AICind
          NIFTY_METAL_MODELS_BIC[[ind]]<-BICind
          NIFTY_METAL_MODELS_AICC[[ind]]<-AICcind
        })
        
        print(ind)
      }
    }
  }
}

NIFTY_METAL_MODELS<-data.frame(matrix(nrow=1296,ncol=7))#1296
columns<-c("pp","qq","PP","QQ","AIC","BIC","AICC")
colnames(NIFTY_METAL_MODELS)<-columns

NIFTY_METAL_MODELS$pp<-as.character(NIFTY_METAL_MODELS_p)
NIFTY_METAL_MODELS$qq<-as.character(NIFTY_METAL_MODELS_q)
NIFTY_METAL_MODELS$PP<-as.character(NIFTY_METAL_MODELS_P)
NIFTY_METAL_MODELS$QQ<-as.character(NIFTY_METAL_MODELS_Q)
NIFTY_METAL_MODELS$AIC<-as.character(NIFTY_METAL_MODELS_AIC)
NIFTY_METAL_MODELS$BIC<-as.character(NIFTY_METAL_MODELS_BIC)
NIFTY_METAL_MODELS$AICC<-as.character(NIFTY_METAL_MODELS_AICC)
#View(NIFTY_METAL_MODELS)
write.csv(NIFTY_METAL_MODELS,file = "METAL_score.csv",sep=",")
#************************************
#*

#************************************

dat<-read.csv('/home/mahendra/Downloads/sem_3/TSA/project/data/METAL_score.csv')
df<-dat%>%select(X,AIC,AICC,BIC)%>%filter(AIC<0)%>%filter(BIC<0)%>%filter(AICC<0)
d <- melt(df, id.vars="X")

ggplot(data=d,
       aes(x=X, y=value, colour=variable)) +
  geom_line()+ labs(x="sl no of different combination of ARIMA and GARCH model",
                    y="score",title = "AIC,BIC and AICc score of different model")



###*************** Best model specification and fitting 
garch_metal <- ugarchspec(mean.model = list(armaOrder=c(0,0)),
                          variance.model = list(model = 'eGARCH', 
                                                garchOrder = c(0,1)),distribution = 'sstd')
fit_garch_metal <- ugarchfit(spec = garch_metal, data= na.omit(as.vector(Return_metal)))
## forecasting
forecast_metal<- ugarchforecast(fit_garch_metal,n.ahead = 30)

par(mfrow=c(1,2))
plot(forecast_metal,which=1)
plot(forecast_metal,which=3)


########### going back to original data 


nifty_Metal <- read.csv('/home/mahendra/Downloads/sem_3/TSA/project/data/Metal_data.csv')
nifty_Metal <- nifty_Metal[852:1251,c(3,7)]
nifty_Metal[,1] <- dmy(nifty_Metal[,1])
original_Metal <- nifty_Metal$Close
Update <- c()
end=original_Metal[382]
for (i in seq(1,18)){
  end= end*exp(forecast_metal@forecast$seriesFor[i])
  print(end)
  Update <- c(Update,end)
}
par(mfrow=c(1,1))
plot(c(1:382),original_Metal[1:382],type="l",col="black",xlim=c(1,420),
     ylim=c(1000,6500),main="Forcasting the original stock value",
     xlab="time point",ylab="close price of nifty-it ",xaxt='n')
lines(c(383:400),original_Metal[383:400],type="l",col="green")
lines(c(383:400),Update,type="p",col="red")
legend("bottomright",legend = c("forecasted stock values","original privious values",
                                "original future ground truths"),
       fill = c("red","black","green"))

## RMSE
nifty_Metal<- read.csv('/home/mahendra/Downloads/sem_3/TSA/project/data/Metal_data.csv')
nifty_Metal <- nifty_Metal[852:1251,c(3,7)]
nifty_Metal[,1] <- dmy(nifty_Metal[,1])
tso_Metal <- zoo(nifty_Metal$Close, nifty_Metal$Date)
Return_Metal=CalculateReturns(tso_Metal, method = 'log')
true_returns_metal <- na.omit(as.vector(Return_Metal))

predicted_returns_metal <- c()
predicted_stocks_metal <- c()
total_sqr_loss_in_returns_metal <- 0
total_sqr_loss_in_stock_metal <- 0
for (i in seq(1,18)){
  fit_garch_metal <- ugarchfit(spec = garch_metal, data = true_returns_metal[1:(381-1+i)] )
  forecast_metal<- ugarchforecast(fit_garch_metal,n.ahead = 1 )
  pred_return=forecast_metal@forecast$seriesFor[1]
  predicted_returns_metal= c(predicted_returns_metal,pred_return)
  sqr_loss_return= (pred_return - true_returns_metal[381+i])^2
  total_sqr_loss_in_returns_metal = total_sqr_loss_in_returns_metal + sqr_loss_return
  
  
  previous_stock = nifty_Metal$Close[382-1+i]
  pred_stock = previous_stock*exp(pred_return)
  predicted_stocks_metal <- c( predicted_stocks_metal, pred_stock)
  print(pred_stock)
  sqr_loss_stock= (pred_stock - nifty_Metal$Close[382+i])^2
  total_sqr_loss_in_stock_metal = total_sqr_loss_in_stock_metal + sqr_loss_stock
}

predicted_returns_metal
predicted_stocks_metal
(total_sqr_loss_in_returns_metal^0.5)/length(predicted_returns_metal)
(total_sqr_loss_in_stock_metal^0.5)/length(predicted_stocks_metal)

plot(c(1:382),original_Metal[1:382],type="l",col="black",xlim=c(1,420),
     ylim=c(1000,6500),main="Forcasting the original stock value",
     xlab="time point",ylab="close price of nifty-it ",xaxt='n')
lines(c(383:400),original_Metal[383:400],type="l",col="green")
lines(c(383:400),predicted_stocks_metal,type="l",col="red")
legend("bottomright",legend = c("forecasted stock values",
                                "original privious values","original future ground truths"),
       fill = c("red","black","green"))

## zoom in
plot(c(1:382),original_Metal[1:382],type="l",col="black",xlim=c(370,420),
     ylim=c(1000,6500),main="Forcasting the original stock value of metal",
     xlab="time point",ylab="close price of nifty-it ",xaxt='n')
lines(c(383:400),original_Metal[383:400],type="b",col="green")
lines(c(383:400),predicted_stocks_metal,type="b",col="red")
legend("bottomright",legend = c("forecasted stock values",
                                "original privious values","original future ground truths"),
       fill = c("red","black","green"))




