library(glmnet)
library(plyr)
set.seed(-1)


# NA'S IN CSV FILE CAUSE BIG BIG BIG PROBLEMS - CHANGE IN EXCEL TO 0'S!!!!!!!!!!
flu_data_raw = read.csv("/Users/f/projects/masters-thesis/src/flu trends/InSample/correlate-ILI-FullSample_lags.csv")
numobs = nrow(flu_data_raw)
fitteddates = flu_data_raw$Date[54:numobs]

#First prepare data set for the lasso - i.e. so that lasso can be implemented

#Delete unnecessary variables
delete = names(flu_data_raw) %in% c("X","Date")
flu_data_1 = flu_data_raw[!delete]

#Delete observations with N/As  // only select observations which have values for all varaibles
#since there are 53 lags, observations should start at 54
#there are numobs in total

flu_data_2 = flu_data_1[54:numobs,]

#Select observations to be used in the lasso - i.e. x and y
y = flu_data_2$ILI
flu_data_3 = flu_data_2
flu_data_2$ILI <- NULL
x = data.matrix(flu_data_2)


#Find the lasso coefficients
#Returns a vector of coefficients including a constant
lasso.mod = glmnet(x, y, alpha=1)
#find best lambda
train = sample(1:nrow(x),nrow(x)/2)
test = (-train)
y.test=y[test]
lasso.mod=glmnet(x[train,],y[train],alpha=1)
set.seed(1)
cv.out=cv.glmnet(x[train,], y[train],alpha=1)
bestlam = cv.out$lambda.min
lasso.pred = predict(lasso.mod, s=bestlam, newx=x[test,])
mean((lasso.pred-y.test)^2)
out = glmnet(x, y, alpha =1)
lasso.coef=predict(out,type="coefficients", s=bestlam)
lasso.coef = matrix(lasso.coef)



#Find fitted values

#First subset lasso.coef to exclude constant
lasso.coef2 = lasso.coef[2:nrow(lasso.coef)]
lasso.coef3 = data.matrix(lasso.coef2)
#so lasso.coef3 is what needs to be multiplied by the observation!


#Create matrix for fitted and true values

results = matrix(0,nrow(flu_data_2),2)

#Now get fill matrix by selecting each observation!
for(j in 1:nrow(results))
{
obs_i = flu_data_2[j,]
obs_i_m = data.matrix(obs_i)

#Now take inner product
dotproduct = obs_i_m%*%lasso.coef3

#Add to lasso.coef constant to get fitted
fitted = lasso.coef[1]+dotproduct

#true value
truth = flu_data_3$ILI[j]

results[j,1] = fitted
results[j,2] = truth

}


df.results = data.frame(results)
df.results = rename(df.results, c("X1"="fitted", "X2"="truth"))
df.results$unformateddates = fitteddates


df.results$Date = strptime(as.character(df.results$unformateddates),"%d/%m/%Y")
df.results$Date = format(df.results$Date, "%Y-%m-%d")
df.results = df.results[c("Date", "fitted", "truth")]

write.csv(df.results, file = "/Users/f/projects/masters-thesis/src/flu trends/InSample/Lasso_results_insample.csv")








