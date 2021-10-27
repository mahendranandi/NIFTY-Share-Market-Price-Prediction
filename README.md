# NIFTY-Share-Market-Price-Prediction
Time series analysis on NIFTY data ( bank,oil,metal,it ) using GARCH model in R.



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

- [ ] It can be difficult to tell if a model is stationary or not. Unlike some visible seasonality , you usually can’t tell by looking at a graph. If you aren’t sure about the stationarity of a model, a hypothesis test can help. You have several options for testing, including:
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

# `[f] About Finance Data and Datasets`

NFI(non financial information) is associated with information that is not expressed in financial terms. NFI is a system of information that does not necessarily derive from the accounting system. NFI is not related to financial and economic data. 

Using non-stationary time series data in financial models produces unreliable and spurious results and leads to poor understanding and forecasting. The solution to the problem is to transform the time series data so that it becomes stationary. If the non-stationary process is a random walk with or without a drift, it is transformed to stationary process by differencing. On the other hand, if the time series data analyzed exhibits a deterministic trend, the spurious results can be avoided by detrending.

Sometimes the non-stationary series may combine a stochastic and deterministic trend at the same time and to avoid obtaining misleading results both differencing and detrending should be applied, as differencing will remove the trend in the variance and detrending will remove the deterministic trend.
### Why Log Returns: 
https://quantivity.wordpress.com/2011/02/21/why-log-returns/
https://medium.datadriveninvestor.com/when-is-log-transformation-necessary-for-financial-returns-4b3f5bb58e62

# `[g] Which model and why?`
GARCH model but why?


2. Data Description
Data taken are on daily basis. Closed prices are only taken for analysis and the currency in
which the stock prices are recorded are in rupees. The time stamp on the data is from 3 rd
January, 2005 to 28 th December, 2020. The data is collected from the url:
http://in.yahoo.finance.com
Finance is a field where time series arises naturally from the evolution of indexes and prices.
So the financial data which I have worked on the Nifty50 daily stock index (i.e. closing
prices) of SBI, HDFC Bank and AXIS Bank. So here we have 3 time series. Now we take one
stock at a time and try to analyze and then fit an appropriate model to forecast.

3. Objective
This report is mainly focused on building up a good model for each stocks i.e. for SBI, HDFC
and AXIS Bank so that we can use these fitted model to forecast.



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





