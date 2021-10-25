# NIFTY-Share-Market-Price-Prediction
Time series analysis on NIFTY data ( bank,oil,metal,it ) using GARCH model in R.



# Content:
- `[a] Introduction `
- `[b] What is Time Series Analysis`
- `[c] Difference from Regression analysis`
- `[d] About Finance Data and Datasets`
- `[e] Stationarity, White Noise, IID`
- `[f] Which model and why?`
- `[g] Steps to follow serially`
- `[h] Data visualization [EDA]`
- `[i] Log returns`
- `[j] Analysis on Log-Returns`
    - `[a] Augmented Dicky Fuller test [Unit root test]`
    - `[b] ACF, PACF of Log-returns `
    - `[c] Mean model [ARIMA] selection `
    - `[D] Observation of the residuals after fitting ARIMA model`
- `[k] Square Log-Returns to observe Volatility`
- `[l] ACF and PACF of the sqr-Log-Returns`
- `[m] Check if Volatility is present`
    - `[a] ARCH test`
    - `[b] Monthly rolling average volatility`
- `[n] GARCH model selection:`
    - `[a] Following the distribution of the log-returns`
    - `[b] Guess about the order of the model ( IF POSSIBLE )`
    - `[c] AIC,AICc,BIC value`
    - `[d] choosing the best model `
- `[o] Forcasting with the best model`
- `[p] References`

# Introduction:
Time series analysis can indeed be used to predict stock trends. The caveat out here is 100%
accuracy in prediction is not possible but still using time series analysis we can develop some
model which will give us an idea or a prediction of how the next few days stock price would
be. In this project we are going to analyze and implement different models step by step in
order to get a model which would be best suited for prediction or forecasting purpose. We
use the log return of the stock prices of each stock of the NIFTY50 (Banking Sector), i.e. SBI,
HDFC and AXIS Banks and then try to fit a traditional model i.e. ARMA model and found the
best model according to the AIC and BIC values and again check whether there is any ARCH
effect or not i.e. to check for the presence of Heteroskedasticity in the data. If present then
we model the variance part through ARCH and GARCH model and found the best mean and
variance model which would capture all the cluster volatility and the bursts in the data and
would forecast appropriately. These are done for each of the stocks and found that SBI and
AXIS bank uses ARCH model to model the variance but in not dependent on the previous
variances in the final model. While in case of HDFC, we have to use GARCH model to model
the variance and it depends on the previous variance to get a good fit. Hence at the end of
the project we get the best models for each stock which are ready to be used to get a good
forecast.
Stock markets are where individual and institutional investors come together to buy and sell
shares in a public venue. Nowadays these exchanges exist as electronic marketplaces. The
supply and demand helps to determine the price for each security or the levels at which
stock market participants - investors and traders - are willing to buy and sell. A stock or
share (also known as a company’s “equity”) is a financial instrument that represents
ownership in a company. There are many indexes out of which NIFTY 50 is a diversified 50
stock index accounting or 13 sectors of the economy. It is used for a variety of purposes
such as benchmarking fund portfolios, index based derivatives and index funds. There are
two main stock exchanges in India that make up the stock markets. One of them is Bombay
Stock Exchange (BSE) and the other one is the National Stock Exchange (NSE). NIFTY 50 is
owned and managed by NSE Indices Limited (formerly known as India Index Services &
Product Limited) (NSE Indices). NSE Indices is India’s specialized company focused upon the
index as a core product. So in this project we are going to do a time series analysis on one of
the sector of the NIFTY 50 i.e. the Banking sector and we are going to take only the daily
closing prices of 3 banks, namely SBI, HDFC Bank and AXIS Bank.


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

# References
- https://www.idrisstsafack.com/post/garch-models-with-r-programming-a-practical-example-with-tesla-stock
