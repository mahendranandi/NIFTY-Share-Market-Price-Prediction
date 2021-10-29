
# `[h] Data visualization [EDA]:`
Finance is a field where time series arises naturally from the evolution of indexes and prices.
So the financial data which I have worked on the Nifty daily stock index (i.e. closing
prices) of IT,METAL,BANK,OIL. So here we have 4 time series. Now we take one
stock at a time and try to analyze and then fit an appropriate model to forecast.

* Libraries used:

```
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
```

* Data Loading and ploting the time series data:
```
nifty_it <- read.csv('/home/mahendra/Downloads/sem_3/TSA/project/data/it_data.csv')
nifty_it <- nifty_it[852:1233,c(3,7)]
dim(nifty_it)          #### (382 x 2)
nifty_it[,1] <- dmy(nifty_it[,1])
#plot(nifty_it$Close, ylab="Stock Prices",main="Figure : Closing prices of the stocks",type = 'l')
tso_it <- zoo(nifty_it$Close, nifty_it$Date)
#plot(tso_it)
it <- data.frame(xts(nifty_it$Close, order.by=as.POSIXct(nifty_it$Date)))
names(it) <- "it closed"
chartSeries(it, type = "line", show.grid = TRUE,name = "CLOSING Price of NIFTY-IT")
```

<img src="./Images/plot_it" align="middle" >

Data taken are on daily basis. Closed prices are only taken for analysis and the currency in
which the stock prices are recorded are in rupees. The time stamp on the data is from 23rd
March, 2020 to 4th october 2021. 




# `[i] Log returns:`
The next step is the calculation of the daily return of the price and display it. For the return calculation we use the function CalculateReturns(). Here is the related code.We use again the function charSeries() in order to display the time series of the returns. 
<img src="./Images/return_it" align="middle" >


Here is the graph of the returns.As we can see there the time series of  returns is almost zero mean(0.003) and the returns displays for some random day very high volatility, meaning that the standard stationarity won't work here.


Now we can display the histogram of returns and try to see if the normal distribution could be used for the conditional error term.







# `[j] Analysis on Log-Returns:`


### `[j][a] Augmented Dicky Fuller test [Unit root test]:`
As we aren’t sure about the stationarity of a model, a hypothesis test can help us,


### `[j][b] ACF, PACF of Log-returns :`


### `[j][c] Mean model [ARIMA] selection: `


### `[j][D] Observation of the residuals after fitting ARIMA model:`







Here Figure 1.1 shows the original time series of the SBI stock index for the period from
January 3, 2005 to December 28, 2020. Note that this index seems to increase with time,
but there are some downward periods commonly denoted as bear markets. As from the
above plot, it is also observed that data is much volatile and there is a sudden fall (crash) in
stock prices starting from March due to COVID.

In order to study these indices, it is customary in finance to consider the logarithm return,
which is defined as ****************&****
where denotes the price or the index value at time t. These returns (i.e. the log returns)
are displayed in Figure 1.2. where we can observe some great drops/bursts during say 2014
and some abrupt changes or great volatility during say 2008 and 2018.

Another look at the volatility is shown in Figure 1.3 where the squared log returns,
are
plotted. From this graph the high volatility of the SBI is evident during the periods
mentioned above.
################################



One of the assumptions of the ARMA model is that the error term are either strongly or weakly stationary. 


The problem is that in real life. This assumption is not always satisfied. Indeed, when looking financial data such as stock market data (AAPL, TESLA, GOOGL) or currency data (EUR/USD, GBP/USD), even indices data ( S&P 500, DAX 30, US30, NASDAQ 100 etc.) and cryptocurrency. Usually these data display an error term which presents a sort of stochastic variation of their volatility over time, meaning that considering the stationarity assumption will lead to a misspecification of the model estimation and therefore will lead to a bad forecast. 


GARCH models are usually the one considering such heteroskedasticity of the error terms and the stochastic change of their related volatility. In this post I will describe a simplified version of the GARCH model, also I will show how to estimate such model setting, how to interpret or read the results and how to find the optimal setting.


### GARCH Model Setting

GARCH stands for Generalized Autoregressive Conditional Heteroskedasticity Models. GARCH models are commonly used to estimate the volatility of returns for stocks, currencies, indices cryptocurrencies. Professional traders use this tool to price assets and detect which asset will potentially provide the best return in their portfolio. Also they can use this tool to adjust their portfolio allocation and risk management. 


There exist a large variety of GARCH model : Standard GARCH (SGARCH), Nonlinear GARCH (NGARCH), Nonlinear Asymetric GARCH (NAGARCH), Integrated GARCH (IGARCH), Exponential GARCH (EGARCH), GARCH in Mean (GARCH-M), Quadratic GARCH (QGARCH), Glosten-Jagannathan-Runke GARCH (GJR-GARCH), Treshold GARCH (TGARCH), Family GARCH (FGARCH), Continuous-time GARCH (COGARCH), Zero drift GARCH (ZDGARCH) etc. I will present only two of these variants : the standard GARCH and the GJR-GARCH models.


### The standard GARCH Model

To model the GARCH model, we need to know first how the ARCH model is set. So let us consider the error term e[t] or the residual from the demeaned return. Then the error term is decomposed into two main terms that are the stochastic term z[t] and the time-dependent standard deviation s[t] such that :


R[t] = mu + e[t]

e[t] = s[t]*z[t]. 

R[t] is the variable representing the time series of the return of the stock considered, mu is the mean and e[t] is the error term. The variable z[t] is assumed to be a strong white noise process. If we consider that q is the number of lags for the volatility modelling (ARCH(q)), then, we have




Therefore, an ARCH(q) model means that the time-dependent volatility depends on the first q lag squared values of the error term.


Then, based on the ARCH(q) model, we can define the model setting of the GARCH. Indeed, the GARCH model is considered when the error variance s[t] is assumed to follow an ARMA process. In that situation, the GARCH(p,q) model with p the number of lags of the s[t] terms and q the number of lags for the ARCH terms e[t]^2. 


Therefore, the main difference between the GARCH model and the ARCH model is that the GARCH model consider also the volatility of the previous period, while the ARCH model do not. This is truly important as in the financial market we can usually observe mean reverting patterns of the instruments/variables and this mean-reverting pattern can in some case could happen by respecting a certain average range, meaning that the volatility of the previous periods should be considered.




Then, a GARCH(1,1) is such that 



and the ARCH(1) model is nothing else than the GARCH(0,1) model.


The particularity of the standard GARCH model is that we consider that the conditional error term follows a normal distribution. This is not always the case for all types of data. We usually observe in the financial data more skewed data.  Therefore, we should also consider checking if the residuals follow that pattern. The GARCH model with skewed student t-distribution (STTD) is usually considered as an alternative to the normal distribution in order to check if we have a better model fitting. 


### Model Estimation

The estimation of the GARCH model is very simple. Indeed considering a GARCH(p,q) model, we have 4 steps :

Estimate the AR(q) model for the returns. and get the residuals e[t]

Construct the time series of the squared residuals, e[t]^2.

Compute and plot the autocorrelation of the squared rediduals e[t]^2.

Estimate  the ARMA (p,q) model for the volatility  s[t] of the residuals based on one of the specified model.




# References
- https://www.idrisstsafack.com/post/garch-models-with-r-programming-a-practical-example-with-tesla-stock
- https://www.tableau.com/learn/articles/time-series-analysis
- https://quantivity.wordpress.com/2011/02/21/why-log-returns/
- https://medium.datadriveninvestor.com/when-is-log-transformation-necessary-for-financial-returns-4b3f5bb58e62
- http://eclr.humanities.manchester.ac.uk/index.php/R_GARCH



******************************************************************************************************************************************************** 
<!-- 

# Content:
- `[a] Introduction `
- `[b] What is Time Series Analysis`
- `[c] Difference from Regression analysis`
- `[d] Stationarity, White Noise, IID`
- `[e] Steps to follow serially`
- `[f] About Finance Data and Datasets`
- `[g] Which model and why?`
- `[h] Data visualization [EDA]`
- `[i] Log returns`
- `[j] Analysis on Log-Returns`
    - `[a] Augmented Dicky Fuller test [Unit root test]`
    - `[b] ACF, PACF of Log-returns `
    - `[c] Mean model [ARIMA] selection `
    - `[D] Observation of the residuals after fitting ARIMA model`
- `[k] Square Log-Returns to observe Volatility`
- `[l] ACF and PACF of the sqr-Log-Returns`
- `[m] Check if Volatility [ARCH effect] is present`
    - `[a] ARCH test`
    - `[b] Monthly rolling average volatility`
- `[n] GARCH model selection:`
    - `[a] Following the distribution of the log-returns`
    - `[b] Guess about the order of the model ( IF POSSIBLE )`
    - `[c] AIC,AICc,BIC value`
    - `[d] choosing the best model `
- `[o] Forcasting with the best model`
- `[p] References`

            “I will tell you how to become rich. Close the doors. Be fearful when others are greedy.
            Be greedy when others are fearful.”                                – By Warren Buffett

###

# `[a] Introduction: `
 In time series analysis, time is a significant variable of the data. Times series analysis helps us study our world and learn how we progress within it. Time series analysis can indeed be used to predict stock trends. Stock markets are where individual and institutional investors come together to buy and sell shares in a public venue. Nowadays these exchanges exist as electronic marketplaces. The supply and demand helps to determine the price for each security or the levels at which stock market participants - investors and traders - are willing to buy and sell. A stock or share (also known as a company’s “equity”) is a financial instrument that represents ownership in a company. There are many indexes out of which NIFTY is a diversified  stock index. It is used for a variety of purposes such as benchmarking fund portfolios, index based derivatives and index funds. There are two main stock exchanges in India that make up the stock markets. One of them is Bombay Stock Exchange (BSE) and the other one is the National Stock Exchange (NSE). NIFTY is owned and managed by NSE Indices Limited (formerly known as India Index Services & Product Limited) (NSE Indices). NSE Indices is India’s specialized company focused upon the index as a core product. In this project we are going to analyze and implement different models step by step in order to get a model which would be best suited for prediction or forecasting purpose. We use the log return of the stock prices of some stock of the NIFTY, i.e. BANK, OIL,IT and METAL Banks and  we are going to take only the daily closing prices of then and then try to fit a traditional model i.e. ARMA model and found the best model according to the AIC and BIC values and again check whether there is any ARCH effect or not i.e. to check for the presence of Heteroskedasticity in the data. If present then we model the variance part through ARCH and GARCH model and found the best mean and variance model which would capture all the cluster volatility and the bursts in the data and would forecast appropriately. The caveat out here is 100% accuracy in prediction is not possible but still using time series analysis we can develop some model which will give us an idea or a prediction of how the next few days stock price would be.

# `[b] What is Time Series Analysis?`

Time series analysis is a statistical technique that deals with time series data, or trend analysis.  Time series data means that data is in a series of  particular time periods or intervals. This is a specific way in which analysts record data points at consistent intervals over a set period of time rather than just recording the data points intermittently or randomly. What sets time series data apart from other data is that the analysis can show how variables change over time. In other words, time is a crucial variable because it shows how the data adjusts over the course of the data points as well as the final results. It provides an additional source of information and a set order of dependencies between the data.
    Time series analysis typically requires a large number of data points to ensure consistency and reliability. An extensive data set ensures you have a representative sample size and that analysis can cut through noisy data. It also ensures that any trends or patterns discovered are not outliers and can account for seasonal variance. Additionally, time series data can be used for forecasting—predicting future data based on historical data.

    Examples of time series analysis in action include:
    Weather data
    Rainfall measurements
    Temperature readings
    Heart rate monitoring (EKG)
    Brain monitoring (EEG)
    Quarterly sales
    Stock prices
    Automated stock trading
    Industry forecasts
    Interest rates

 
# `[c] Difference from Regression analysis :`
We have to follow first that Regression is a mathematical model to make relation between variables, and it is being used in Time series analysis also to remove trend. So, regression helps us to get the trend, and after removing trend( with the help of regression) and seasonality by some other method( like differencing) we have to check wheathere the series is Stationary(discussed later) or not and then we can approach for the time series model accordingly. we can even think of time series as an extension of linear regression. Time series uses terms such as autocorrelation and moving average to summarize historical information of the y variable with the hope that these features better predict future y. So, there is a difference one should be clear about.
Again, Regression can also be applied to non-ordered series where a target variable is dependent on values taken by other variables. These other variables are called as Features. When making a prediction, new values of Features are provided and Regression provides an answer for the Target variable.So unlike time series analysis, we need some feature to predic. Essentially, Regression is a kind of intrapolation technique. Regression can be applied to Time-series problems as well. e.g. Auto-regression.
So, in a simple way,
* Time-series forecast is Extrapolation [ out side the given data ] and realation between target variable and time.
* Regression is Intrapolation [ inside the given data also ] and relation between target variable and features.

# `[d] Stationarity, White Noise, IID`

- [ ] `Stationarity:`
A time series has stationarity if a shift in time doesn't cause a change in the shape of the distribution. Basic properties of the distribution like the mean , variance and covariance are constant over time. In the most intuitive sense, stationarity means that the statistical properties of a process generating a time series do not change over time . It does not mean that the series does not change over time, just that the way it changes does not itself change over time. A common assumption in many time series techniques is that the data are stationary.Stationarity can be defined in precise mathematical terms, but for our purpose we mean a flat looking series, without trend, constant variance over time, a constant autocorrelation structure over time and no periodic fluctuations (seasonality).\
- [ ] Data points are often non-stationary or have means, variances, and covariances that change over time. Non-stationary behaviors can be trends, cycles, random walks, or combinations of the three.Non-stationary data, as a rule, are unpredictable and cannot be modeled or forecasted. In order to receive consistent, reliable results, the non-stationary data needs to be transformed into stationary data with one of the following techniques:
```
** non stationarity in mean:
        * deterministic trend
               * detranding
        * stochastic trend
               * differencing
```
* We can difference the data. That is, given the series Zt, we create the new series
Yi=Zi−Zi−1.
The differenced data will contain one less point than the original data. Although you can difference the data more than once, one difference is usually sufficient.
* If the data contain a trend, we can fit some type of curve to the data and then model the residuals from that fit. Since the purpose of the fit is to simply remove long term trend, a simple fit, such as a straight line, is typically used.
* For non-constant variance, taking the logarithm or square root of the series may stabilize the variance. For negative data, you can add a suitable constant to make all the data positive before applying the transformation. This constant can then be subtracted from the model to obtain predicted (i.e., the fitted) values and forecasts for future points.
The above techniques are intended to generate series with constant location and scale. Although seasonality also violates stationarity, this is usually explicitly incorporated into the time series model.

- [ ] Types of Stationary
Models can show different types of stationarity:

**Strict stationarity** means that the joint distribution of any moments of any degree (e.g. expected values, variances, third order and higher moments) within the process is never dependent on time. This definition is in practice too strict to be used for any real-life model.
**First-order** stationarity series have means that never changes with time. Any other statistics (like variance) can change.
**Second-order** stationarity (also called weak stationarity) time series have a constant mean, variance and an autocovariance that doesn’t change with time. Other statistics in the system are free to change over time. This constrained version of strict stationarity is very common.
**Trend-stationary** models fluctuate around a deterministic trend (the series mean). These deterministic trends can be linear or quadratic, but the amplitude (height of one oscillation) of the fluctuations neither increases nor decreases across the series.
**Difference-stationary** models are models that need one or more differencings to become stationary (see Transforming Models below).

- [ ] It can be difficult to tell if a model is stationary or not. Unlike some visible seasonality , you usually can’t tell by looking at a graph. If we aren’t sure about the stationarity of a model, a hypothesis test can help. we have several options for testing, including:
* Unit root tests (e.g. Augmented Dickey-Fuller (ADF) test or Zivot-Andrews test),
* A KPSS test (run as a complement to the unit root tests).
* A run sequence plot,
* The Priestley-Subba Rao (PSR) Test or Wavelet-Based Test, which are less common tests based on spectrum analysis.\
Though we will use only the Unit root test here. (To know more about it)[https://www.investopedia.com/articles/trading/07/stationary.asp]

- [ ] `White Noise and IID:`
A white noise process is only defined by the first 2 moments. A noise sequence (et) is a white noise sequence if\
the expectation of each element is zero, E(et) = 0\
the variance of each element is finite, Var(et) < infinity\
the elements are uncorrelated, Cor(et, es) = 0\
But it does not specify higher moments of the distribution, like skewness and kurtosis. IID white noise provides that the sample has the same distribution, so also higher moments have to be the same. The noise sequence would be an iid noise if in addition the elements are not just uncorrelated but also independet. So therefore every iid noise is also white noise, but the reverse is just true for Gaussian white noise sequence. A Gaussian white noise implies a normal distribution of et and a normal distribution is completely defined by the first 2 moments. So in this case: White noise process = Iid white noise.  IID is a special case of white noise. So, the difference is that for iid noise we assume each sample has the same probability distribution while, white noise samples could follow different probability distribution.
The concept of iid is used when we make assumptions about the error, e.g. in regression analysis when we say that the error terms are iid following normal with mean 0 and a common variance sigma^2. However the concept of white noise is used in time series analysis, when we make more complicated models like random walk or ARMA or ARIMA models.

# `[e] Steps to follow serially:`
The following steps are to be followed:
1. Visualization of data
2. Removing trend and seasonality
3. chechiking stationarity of the residuals
4. fitting best ARIMA model
5. if residuals are stationary but there is still volatility is present check if ARCH effect is present or not
6. if present then to model the variance use GARCH model
7. choosing best ARMIA + GARCH model to model mean and variance at the same time.


# NIFTY-Share-Market-Price-Prediction
Time series analysis on NIFTY data ( bank,oil,metal,it ) using GARCH model in R.




# `[f] About Finance Data and Datasets`
A financial information is a formal record of the financial activities of a business, person, or other entity. Relevant financial information is presented in a structured manner and in a form easy to understand. Financial data consists of pieces or sets of information related to the financial health of a business. The pieces of data are used by internal management to analyze business performance and determine whether tactics and strategies must be altered. People and organizations outside a business will also use financial data reported by the business to judge its credit worthiness, decide whether to invest in the business, and determine whether the business is complying with government regulations.  
NFI(non financial information) is associated with information that is not expressed in financial terms. NFI is a system of information that does not necessarily derive from the accounting system. NFI is not related to financial and economic data. 

### When is Log Transformation Necessary for Financial Returns?
There are a lot of benefits in using the log return or compounded return over the simple one. Some of these are lognormality, raw log equality, and low algorithmic complexity. Finally, use log return when temporal behavior of return is the focus of interest.
Returns for stock prices are normally distributed but prices are not. They are lognormally distributed (assumed at least and require verification per case).
let us discuss when should you use simple returns over compounded returns

Using non-stationary time series data in financial models produces unreliable and spurious results and leads to poor understanding and forecasting. The solution to the problem is to transform the time series data so that it becomes stationary. If the non-stationary process is a random walk with or without a drift, it is transformed to stationary process by differencing. On the other hand, if the time series data analyzed exhibits a deterministic trend, the spurious results can be avoided by detrending.
Sometimes the non-stationary series may combine a stochastic and deterministic trend at the same time and to avoid obtaining misleading results both differencing and detrending should be applied, as differencing will remove the trend in the variance and detrending will remove the deterministic trend.
### Why Log Returns: 
- [ ] `LINEAR RETURN VS. COMPOUND RETURNS:`
The linear return is defined as-

<!-- ![\Large f(x,y)=$\frac{sin(10(x^2+y^2))}{10}](https://latex.codecogs.com/svg.latex?\Large&space;f(x,y)=\frac{sin(10(x^2+y^2))}{10})  -->

![ L_t=$\frac{P_{t+1}}{P_t}-1](https://latex.codecogs.com/svg.latex?&space;L_t=\frac{P_{t+1}}{P_t}-1) 

Linear Return has the property of being asset-additive that is, you can aggregate the returns more easily. In equation form, if we denote w_1,w_2,w_3,...,w_n as the corresponding weights of n securities, the portfolio returns are simply:
![ L_t=$\frac{P_{t+1}}{P_t}-1](https://latex.codecogs.com/svg.latex?&space;L_{t,p}=w_1P_{t,1}+w_2P_{t,2}+...+w_nP_{t,n})
Where L(t,p) = Portfolio Return using Linear Returns\
Linear returns are therefore being used by risk and portfolio managers for risk analysis, performance attribution, and portfolio optimization.
Compound Returns, however, are calculated using the following formula:

<!-- ![ ](https://latex.codecogs.com/svg.latex?&space;C_t=ln(\frac{P_{t+1}}{P_t})= C_t+C_{t+1}+...+C_{t+k-1} ) -->

![image](https://user-images.githubusercontent.com/76218559/139386164-4b16e198-ada0-4ab1-817c-2e5705c01ab7.png)


Compounded returns are time-additive that is, you can add them across time to get the total return over a specified period. Others benifits are-\
**LOGNORMALITY**
A common assumption for most assets or stocks is that their prices are log-normally distributed. One reason for this is that prices cannot assume a negative value and as the stock price goes closer to 0, the movement becomes smaller.

![image](https://user-images.githubusercontent.com/76218559/139386410-d17f537d-ec11-4fa1-85a5-86edafd0d78e.png)

Formulation of lognormal distribution of prices. The reason for this is that prices change over time and so is the variance of the price. So instead of using the marginal distribution, we are using the conditional distribution.
What does this mean? If the log of price is normally distributed, then:
![image](https://user-images.githubusercontent.com/76218559/139386823-e316792b-9297-4f57-9fde-004c8e40cd75.png)

A happy benefit of this is that there are a lot of tried-and-tested tools, theories, and methods that can be applied when a variable is normally distributed depending on your objective.\
**RAW-LOG RETURN EQUALITY**
When (ri) is small, then:
![image](https://user-images.githubusercontent.com/76218559/139387389-c99b35d9-8f6b-4a58-9459-8851713bcfe8.png)

The approximation is considered good and is relevant for small trades.\
**ALGORITHMIC COMPLEXITY**
Assume for example that a stock is traded, n-times, on a daily basis. To get the compound return:
![image](https://user-images.githubusercontent.com/76218559/139387534-08959af5-d7a9-460f-8a11-f1da44110fae.png)

*We know that the product of two normal variables is not normally distributed.
The log of return, however, resolves these difficulties naturally.*
Using the formula for the log of returns:
![image](https://user-images.githubusercontent.com/76218559/139387579-deea9261-ff80-4c48-8648-17edefef4604.png)

So, to calculate an n-week (or daily) return, we can apply the formula:
![image](https://user-images.githubusercontent.com/76218559/139387631-f5771e9a-e9b4-4fdb-99d3-9bc52faa94be.png)

We know that we can decompose the equation to:
![image](https://user-images.githubusercontent.com/76218559/139387648-b049db8d-1aac-4b26-906e-b212ab00bc91.png)

Thus, the algorithmic complexity is reduced from `O(n)` multiplications to just an `O(1)` addition.
This is extremely helpful for large (n)s.
Furthermore, probability theory tells us that the sum of normally distributed variables is a normally distributed variable itself.


### `About The Datasets:`
We shall take 4 stock data named NIFTY-bank, NIFTY-oil, NIFTY-metal and NIFTY-it. The data is collected from [here](https://www1.nseindia.com/products/content/equities/indices/historical_index_data.htm). Each of them contains 7 different column as shown:
[image here]<img src="./Images/data_it.png" align="middle">

we will take data from 23rd march 2020 to 4th october 2021. Because there due to codiv-19 there was a sudden fall in all stock. So modeling before that is less significant. We will take only closiing price for modeling. 


# `[g] Which model and why?`
One of the assumptions of the ARMA model is that the error term are either strongly or weakly stationary. But the problem is that in real life this assumption is not always satisfied. Indeed, when looking financial data such as stock market data (AAPL, TESLA, GOOGL) or currency data (EUR/USD, GBP/USD), even indices data ( S&P 500, DAX 30, US30, NASDAQ 100 etc.) and cryptocurrency. Usually these data display an error term which presents a sort of stochastic variation of their volatility over time, meaning that considering the stationarity assumption will lead to a misspecification of the model estimation and therefore will lead to a bad forecast. GARCH models are usually the one considering such heteroskedasticity of the error terms and the stochastic change of their related volatility.
GARCH stands for Generalized Autoregressive Conditional Heteroskedasticity Models. GARCH models are commonly used to estimate the volatility of returns for stocks, currencies, indices cryptocurrencies.

There exist a large variety of GARCH model : Standard GARCH (SGARCH), Nonlinear GARCH (NGARCH), Nonlinear Asymetric GARCH (NAGARCH), Integrated GARCH (IGARCH), Exponential GARCH (EGARCH), GARCH in Mean (GARCH-M), Quadratic GARCH (QGARCH), Glosten-Jagannathan-Runke GARCH (GJR-GARCH), Treshold GARCH (TGARCH), Family GARCH (FGARCH), Continuous-time GARCH (COGARCH), Zero drift GARCH (ZDGARCH) etc. I will present only two of these variants : the standard GARCH and the GJR-GARCH models.

**The standard GARCH Model**

To model the GARCH model, we need to know first how the ARCH model is set. So let us consider the error term e[t] or the residual from the demeaned return. Then the error term is decomposed into two main terms that are the stochastic term z[t] and the time-dependent standard deviation s[t] such that :
R[t] = mu + e[t]
e[t] = s[t]*z[t]. 
Where R[t] is the variable representing the time series of the return of the stock considered, mu is the mean and e[t] is the error term. The variable z[t] is assumed to be a strong white noise process. If we consider that q is the number of lags for the volatility modelling (ARCH(q)), then, we have

![](https://latex.codecogs.com/svg.latex?\Large&space;s_{t}^{2}=\alpha_{0}+\alpha_{1}e_{t-1}^{2}+...+\alpha_{q}e_{t-q}^{2})
where alpha_0 >0 and alpha_i >= 0 for i>0

Therefore, an ARCH(q) model means that the time-dependent volatility depends on the first q lag squared values of the error term.


Then, based on the ARCH(q) model, we can define the model setting of the GARCH. Indeed, the GARCH model is considered when the error variance s[t] is assumed to follow an ARMA process. In that situation, the GARCH(p,q) model with p the number of lags of the s[t] terms and q the number of lags for the ARCH terms e[t]^2. 


Therefore, the main difference between the GARCH model and the ARCH model is that the GARCH model consider also the volatility of the previous period, while the ARCH model do not. This is truly important as in the financial market we can usually observe mean reverting patterns of the instruments/variables and this mean-reverting pattern can in some case could happen by respecting a certain average range, meaning that the volatility of the previous periods should be considered.

![](https://latex.codecogs.com/svg.latex?&space;R_t=x_t^`b+e_t)

![](https://latex.codecogs.com/svg.latex?&space;e_t/b=N(0,s_t))

![](https://latex.codecogs.com/svg.latex?&space;s_{t}^{2}=\omega+\alpha_{1}e_{t-1}^{2}+...+\alpha_{q}e_{t-q}^{2}+\beta_{1}s_{t-1}^{2}+...+\beta_{p}s_{t-p}^{2})

Then, a GARCH(1,1) is such that 
![](https://latex.codecogs.com/svg.latex?&space;s_{t}^{2}=\omega+\alpha_{1}e_{t-1}^{2}+\beta_{1}s_{t-1}^{2})
and the ARCH(1) model is nothing else than the GARCH(0,1) model.

The particularity of the standard GARCH model is that we consider that the conditional error term follows a normal distribution. This is not always the case for all types of data. We usually observe in the financial data more skewed data.  Therefore, we should also consider checking if the residuals follow that pattern. The GARCH model with skewed student t-distribution (STTD) is usually considered as an alternative to the normal distribution in order to check if we have a better model fitting. 

 -->

