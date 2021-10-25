# NIFTY-Share-Market-Price-Prediction
Time series analysis on NIFTY data ( bank,oil,metal,it ) using GARCH model in R.



Content:
- [a] Introduction 
- [b] What is Time Series Analysis
- [c] Difference from Regression analysis
- [d] About Finance Data and Datasets
- [e] Stationarity, White Noise, IID
- [f] Which model and why?
- [g] Steps to follow serially
- [h] Data visualization [EDA]
- [i] Log returns
- [j] Analysis on Log-Returns
    - [a] Augmented Dicky Fuller test [Unit root test]
    - [b] ACF, PACF of Log-returns 
    - [c] Mean model [ARIMA] selection 
    - [D] Observation of the residuals after fitting ARIMA model
- [k] Square Log-Returns to observe Volatility
- [l] ACF and PACF of the sqr-Log-Returns
- [m] Check if Volatility is present
    - [a] ARCH test
    - [b] Monthly rolling average volatility
- [n] GARCH model selection:
    - [a] Following the distribution of the log-returns
    - [b] Guess about the order of the model ( IF POSSIBLE )
    - [c] AIC,AICc,BIC value
    - [d] choosing the best model 
    
- [o] Forcasting with the best model

