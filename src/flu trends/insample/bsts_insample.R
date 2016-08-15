#loading required libraries

library(plyr)
#library(glmnet)
library(bsts)
library(BoomSpikeSlab)
library(ggplot2)
library(reshape2)

#loading fludata
flu_data_untouched = read.csv("/Users/f/projects/masters-thesis/src/flu trends/InSample/correlate-ILI-FullSample_lags.csv")
flu_data_raw = read.csv("/Users/f/projects/masters-thesis/src/flu trends/InSample/correlate-ILI-FullSample_lags.csv")
#clean up data/delete unnecessary variables
delete1 = names(flu_data_raw) %in% c("index", "Date")


flu_data_raw = flu_data_raw[!delete1]

#get rid of first column which has another date index
numvar = ncol(flu_data_raw)
numobs = nrow(flu_data_raw)
flu_data = flu_data_raw[54:numobs,2:numvar]

#Here we are going to get predictions for the entire sample
predict_data = flu_data_raw[54:numobs,]

#set up model to be used for BSTS
y = flu_data$ILI
flu_data$ILI <- NULL
flu_data_x = flu_data
x = data.matrix(flu_data_x)



#Number of nowcasts to make is N
# N = 10
ss <- AddLocalLinearTrend(state.specification=NULL, y)
bsts.model <- bsts(y~x, state.specification = ss, niter = 100, 		ping=0, seed=1)	

burn = SuggestBurn(0.1, bsts.model)
time = index(bsts.model$original.series)
errors <- bsts.prediction.errors(bsts.model, burn = burn)
forecast <- t(as.numeric(bsts.model$original.series) - t(errors))
                     
ave_errors = colMeans(errors)
ave_forecasts = colMeans(forecast)
errors_squared = ave_errors^2
MSE = mean(errors_squared)

numfitted = numobs - 54
dates = flu_data_untouched$Date[55:numobs]

df.graph = data.frame(flu_data_untouched$Date[55:numobs], bsts.model$original.series[1:numfitted], ave_forecasts[1:numfitted])
df.graph = rename(df.graph, c("flu_data_untouched.Date.55.numobs."="Date", "bsts.model.original.series.1.numfitted."="ILI", "ave_forecasts.1.numfitted."="Fitted"))

write.csv(df.graph, file = "/Users/f/projects/masters-thesis/src/flu trends/InSample/BSTS_results.csv")
difference = abs((bsts.model$original.series[400:500] - ave_forecasts[400:500]))
prediction = predict(bsts.model, newdata = predict_data)
difference <- prediction$original.series - prediction$mean



