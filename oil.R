## oil data

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

#----------loading data  oil

#####
nifty_oil <- read.csv('/home/mahendra/Downloads/sem_3/TSA/project/data/oil_data.csv')
nifty_oil <- nifty_oil[c(849:1200,seq(1202,1260,2)),c(3,7)]         # 849 : 1200 1202,1204,..,1260
dim(nifty_oil) # 382   2
nifty_oil[,1] <- dmy(nifty_oil[,1])
tso_oil <- zoo(nifty_oil$Close, nifty_oil$Date)
######################## EDA ############################
oil <- data.frame(xts(nifty_oil$Close, order.by=as.POSIXct(nifty_oil$Date)))
names(oil) <- "oil closed"
chartSeries(oil, type = "line", show.grid = TRUE,name = "CLOSING Price of NIFTY-OIL data")
####################
### log Return #####
### sqr return #####
####################
Return_oil=CalculateReturns(tso_oil, method = 'log')
return_oil <- data.frame(xts(Return_oil, order.by=as.POSIXct(nifty_oil$Date)))
chartSeries(return_oil, type = "line", show.grid = TRUE,name = "Log-returns of NIFTY-oil")

rtrn_oil=Return_oil[-c(1),] # remove the first row as oil does not contain a value
chart_Series(rtrn_oil)
# histogram of the returns
chart.Histogram(return_oil,methods = c("add.density","add.normal"),colorset = c("blue","red","black")
                ,main = "histogram of the log-returns of oil data")
legend("topright",legend = c("return","kernel","normal dist"),fill = c("blue","red","black"))

sqr_Return_oil = Return_oil^2
sqr_return_oil <- data.frame(xts(sqr_Return_oil, order.by=as.POSIXct(nifty_oil$Date)))
chartSeries(sqr_return_oil, type = "line", show.grid = TRUE,name = "square of Log-returns of NIFTY-oil")
#####################################
#### Augmented Dickey Fuller Test ###
#######  ADF of returns  ############
#####################################
summary(ur.df(na.omit(Return_oil)))

#####################################
##### ACF of return #####
#### PACF of return #####
#####################################
a<- ggAcf(na.omit(as.vector(Return_oil)), col='blue',main='Acf of  Log-Return of NIFTY-oil data')
p<- ggPacf(na.omit(as.vector(Return_oil)),col='magenta',main='PAcf of  Log-Return of NIFTY-oil data')
grid.arrange(a,p, ncol = 2, nrow = 1)
############# Identifying the mean model by ARIMA #########################
arima_oil <- auto.arima(na.omit(as.vector(Return_oil)))
arima_oil
checkresiduals(arima_oil)

# adf test of the residual
summary(ur.df(resid(arima_oil),type="none",lag=1))

#  ********************************************************* NOW GARCH 

#################################################################
# Squared of Return are auto correlated.
# Squared of Return acf an pacf
##################################################################
c <- ggAcf(na.omit(as.vector(Return_oil))^2, lag.max = 40, col='green', main='ACF of squared Return Values of the oil data')
d<- ggPacf(na.omit(as.vector(Return_oil))^2,lag.max = 40, col='steelblue',main= 'PACF of squared Return Values of the oil data')
grid.arrange(c,d, ncol = 2, nrow = 1)


############################################
# Testing ARCH   ##########################
############################################
library(FinTS)
ArchTest(Return_oil,lags=1,demean = TRUE)

###############################
#### Volatility Clustering ####
###############################
sq_residual_oil <- arima_res_oil^2
ggtsdisplay(sq_residual_oil,main="Squared Residuals after fitting best ARIMA")

chart.RollingPerformance(na.omit(Return_oil),width = 22,FUN = 'sd.annualized',scale=252, main = 'Rolling 1 month Volatility of the log-return of oil data')


############## QQ Plot ##############
ggplot(data=nifty_oil, aes(sample = as.vector(Return_oil))) +
  stat_qq() +
  stat_qq_line(col='red') + ggtitle('QQ plot of Nifty-oil Returns')

######################################
Box.test(na.omit(as.vector(Return_oil)),  lag = 1, type = "Ljung-Box", fitdf = 0)
######################################


##################################################################################
################################# GARCH Model ####################################
##################################################################################


#**********************************
NIFTY_OIL_MODELS_p<-list()
NIFTY_OIL_MODELS_q<-list()
NIFTY_OIL_MODELS_P<-list()
NIFTY_OIL_MODELS_Q<-list()
NIFTY_OIL_MODELS_AIC<-list()
NIFTY_OIL_MODELS_BIC<-list()
NIFTY_OIL_MODELS_AICC<-list()

ind=0
for (p in seq(0,5)){
  for (q in seq(0,5)){
    for (P in seq(0,5)){
      for (Q in seq(0,5)){
        try({
          spec <- ugarchspec(mean.model = list(armaOrder=c(p,q)),
                             variance.model = list(model = 'eGARCH',
                                                   garchOrder = c(P,Q)),distribution = 'std')
          fit <- ugarchfit(spec = spec, data= na.omit(Return_oil)) 
          k=p+q+P+Q
          n=382
          
          AICind<-infocriteria(fit)[1]
          BICind<-infocriteria(fit)[2]
          AICcind <- AICind + (2*k*(k+1)/(n-k-1))
        })
        
        
        ind=ind+1
        NIFTY_OIL_MODELS_p[[ind]]<-p
        NIFTY_OIL_MODELS_q[[ind]]<-q
        NIFTY_OIL_MODELS_P[[ind]]<-P
        NIFTY_OIL_MODELS_Q[[ind]]<-Q
        try({
          NIFTY_OIL_MODELS_AIC[[ind]]<-AICind
          NIFTY_OIL_MODELS_BIC[[ind]]<-BICind
          NIFTY_OIL_MODELS_AICC[[ind]]<-AICcind
        })
        
        print(ind)
      }
    }
  }
}

NIFTY_OIL_MODELS<-data.frame(matrix(nrow=1296,ncol=7))#1296
columns<-c("pp","qq","PP","QQ","AIC","BIC","AICC")
colnames(NIFTY_OIL_MODELS)<-columns

NIFTY_OIL_MODELS$pp<-as.character(NIFTY_OIL_MODELS_p)
NIFTY_OIL_MODELS$qq<-as.character(NIFTY_OIL_MODELS_q)
NIFTY_OIL_MODELS$PP<-as.character(NIFTY_OIL_MODELS_P)
NIFTY_OIL_MODELS$QQ<-as.character(NIFTY_OIL_MODELS_Q)
NIFTY_OIL_MODELS$AIC<-as.character(NIFTY_OIL_MODELS_AIC)
NIFTY_OIL_MODELS$BIC<-as.character(NIFTY_OIL_MODELS_BIC)
NIFTY_OIL_MODELS$AICC<-as.character(NIFTY_OIL_MODELS_AICC)
#View(NIFTY_OIL_MODELS)
write.csv(NIFTY_OIL_MODELS,file = "OIL_score.csv",sep=",")
#********************************
#*

#************************************

dat<-read.csv('/home/mahendra/Downloads/sem_3/TSA/project/data/OIL_score.csv')
dat<- dat[2:864,]
length(dat$X)
df<-dat%>%select(X,AIC,AICC,BIC)#%>%filter(AIC<0)#%>%filter(BIC<0)%>%filter(AICC<0)
d <- melt(df, id.vars="X")
d

ggplot(data=d,
       aes(x=X, y=value, colour=variable)) +
  geom_line()+ labs(x="sl no of different combination of ARIMA and GARCH model",
                    y="score",title = "AIC,BIC and AICc score of different model")



###*************** Best model specification and fitting 
garch_oil <- ugarchspec(mean.model = list(armaOrder=c(2,3)),
                        variance.model = list(model = 'sGARCH', 
                                              garchOrder = c(4,3)),distribution = 'std')
fit_garch_oil <- ugarchfit(spec = garch_oil, data= na.omit(as.vector(Return_oil)))

## forecasting
forecast_oil<- ugarchforecast(fit_garch_oil,n.ahead = 30)

par(mfrow=c(1,2))
plot(forecast_oil,which=1)
plot(forecast_oil,which=3)


########### going back to original data 


nifty_Oil <- read.csv('/home/mahendra/Downloads/sem_3/TSA/project/data/Oil_data.csv')
nifty_Oil <- nifty_Oil[849:1248,c(3,7)]
nifty_Oil[,1] <- dmy(nifty_Oil[,1])
original_Oil <- nifty_Oil$Close
Update <- c()
end=original_Oil[382]
for (i in seq(1,18)){
  end= end*exp(forecast_oil@forecast$seriesFor[i])
  print(end)
  Update <- c(Update,end)
}
par(mfrow=c(1,1))
plot(c(1:382),original_Oil[1:382],type="l",col="black",xlim=c(1,420),
     ylim=c(1000,10000),main="Forcasting the original stock value",
     xlab="time point",ylab="close price of nifty-it ",xaxt='n')
lines(c(383:400),original_Oil[383:400],type="l",col="green")
lines(c(383:400),Update,type="p",col="red")
legend("bottomright",legend = c("forecasted stock values","original privious values",
                                "original future ground truths"),
       fill = c("red","black","green"))

## RMSE
nifty_Oil<- read.csv('/home/mahendra/Downloads/sem_3/TSA/project/data/Oil_data.csv')
nifty_Oil <- nifty_Oil[852:1251,c(3,7)]
nifty_Oil[,1] <- dmy(nifty_Oil[,1])
tso_Oil <- zoo(nifty_Oil$Close, nifty_Oil$Date)
Return_Oil=CalculateReturns(tso_Oil, method = 'log')
true_returns_oil <- na.omit(as.vector(Return_Oil))

predicted_returns_oil <- c()
predicted_stocks_oil <- c()
total_sqr_loss_in_returns_oil <- 0
total_sqr_loss_in_stock_oil <- 0
for (i in seq(1,18)){
  fit_garch_oil <- ugarchfit(spec = garch_oil, data = true_returns_oil[1:(381-1+i)] )
  forecast_oil<- ugarchforecast(fit_garch_oil,n.ahead = 1 )
  pred_return=forecast_oil@forecast$seriesFor[1]
  predicted_returns_oil= c(predicted_returns_oil,pred_return)
  sqr_loss_return= (pred_return - true_returns_oil[381+i])^2
  total_sqr_loss_in_returns_oil = total_sqr_loss_in_returns_oil + sqr_loss_return
  
  
  previous_stock = nifty_Oil$Close[382-1+i]
  print(previous_stock)
  pred_stock = previous_stock*exp(pred_return)
  predicted_stocks_oil <- c( predicted_stocks_oil, pred_stock)
  print(pred_stock)
  sqr_loss_stock= (pred_stock - nifty_Oil$Close[382+i])^2
  total_sqr_loss_in_stock_oil = total_sqr_loss_in_stock_oil + sqr_loss_stock
}
nifty_Oil$Close[382]

predicted_returns_oil
predicted_stocks_oil
(total_sqr_loss_in_returns_oil^0.5)/length(predicted_returns_oil)
(total_sqr_loss_in_stock_oil^0.5)/length(predicted_stocks_oil)

plot(c(1:382),original_Oil[1:382],type="l",col="black",xlim=c(1,420),
     ylim=c(10000,45000),main="Forcasting the original stock value",
     xlab="time point",ylab="close price of nifty-it ",xaxt='n')
lines(c(383:400),original_Oil[383:400],type="l",col="green")
lines(c(383:400),predicted_stocks_oil,type="l",col="red")
legend("bottomright",legend = c("forecasted stock values",
                                "original privious values","original future ground truths"),
       fill = c("red","black","green"))

## zoom in
plot(c(1:382),original_Oil[1:382],type="l",col="black",xlim=c(370,420),
     ylim=c(10000,45000),main="Forcasting the original stock value of oil",
     xlab="time point",ylab="close price of nifty-it ",xaxt='n')
lines(c(383:400),original_Oil[383:400],type="b",col="green")
lines(c(383:400),predicted_stocks_oil,type="b",col="red")
legend("bottomright",legend = c("forecasted stock values",
                                "original privious values","original future ground truths"),
       fill = c("red","black","green"))

