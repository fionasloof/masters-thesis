## This R file produces nowcasts using the lasso algorithm ##
# NOTE: NA'S IN CSV FILE CAUSE BIG BIG BIG PROBLEMS - CHANGE IN EXCEL TO 0'S!!!!!!!!!!

# load required libraries
	library(glmnet)
	library(plyr)
	set.seed(-1)


	flu_data_raw = read.csv("/Users/f/projects/masters-thesis/src/flu trends/nowcast/correlate-ILI-ExcludingHO_lags.csv")
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


#Create results matrix for fitted and true values
	results = matrix(0,52,2)

#Select observations to be used in the lasso - i.e. x and y
#When doing nowcasting can only use up until t-2 observations to select the model to be used
#Selecting model ONCE - then using those same variables re-estimate with lambda equal=0 
#Use recursive estimation - so that each new model estimated contains an additional observation 
#Not a rolling window
#We want nowcasts for a year's worth of observations
#The holdout sample is: week52015-week42016 (Feb 01 2015 - Jan 24 2016)
#So the first prediction is for Feb 01 2015.
#Because data comes out at a two week lag, use the observations up unti 2 weeks earlier to estimate model
#Total number of forecasts to make is 52 
#So numobs-52 is the first prediction and uses the data up until numobs-52-2

#First use Lasso to select relevant variables
#New parameters will be estimated for these variables for the remaining nowcasts
	lastobs_formodel = numobs-52+1-2
#52 because of 52 predicitons, +1 because each time a new model is estimated which includes a new observation, -2 because
#only observations up until 2 periods prior can be used
#Adding 1 is correct - double and triple checked!


	flu_data_model = flu_data_1[54:lastobs_formodel,]
	y = flu_data_model$ILI
	flu_data_3 = flu_data_2
	flu_data_model$ILI <- NULL
	x = data.matrix(flu_data_model)

#Do similar steps but to include the data for the observation we are trying to fit
	nowcastpoint = numobs-52+1
	flu_data_4 = flu_data_1
	flu_data_4$ILI <- NULL
	flu_data_data = flu_data_4[nowcastpoint,]


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
	lasso.coef_raw=predict(out,type="coefficients", s=bestlam)
	lasso.coef = matrix(lasso.coef_raw)



#Find fitted values

#First subset lasso.coef to exclude constant
	lasso.coef2 = lasso.coef[2:nrow(lasso.coef)]
	lasso.coef3 = data.matrix(lasso.coef2)
#so lasso.coef3 is what needs to be multiplied by the observation!

#obs_i = flu_data_2[1,]
	obs_i_m = data.matrix(flu_data_data)

#Now take inner product
	dotproduct = obs_i_m%*%lasso.coef3

#Add to lasso.coef constant to get fitted
	fitted = lasso.coef[1]+dotproduct

#true value
	truth = flu_data_3$ILI[numobs-54+1-52+1]

	results[1,1] = fitted
	results[1,2] = truth


#Now take those variables which were selected and basically just use OLS on them
	paranames =dimnames(lasso.coef_raw)

	chosenvarlist = character()

	for(h in 2:nrow(lasso.coef))
	{
		if(lasso.coef[h] > 0 || lasso.coef[h])
		{
		chosenvarlist = c (chosenvarlist, paranames[[1]][h])
		}
		else
		{
		}
	}


	for(k in 1:52)
	{
	newdf = flu_data_2[chosenvarlist]
	newdf$ILI = flu_data_2$ILI
#select sample from which to use OLS to model
#260 observations in total
	lastvalidobs = 260-52+k-2
	df.modelwith = newdf[1:lastvalidobs,]
	model = lm(ILI~., data=df.modelwith)
	newD =newdf[260-52+k,]
	fitted = predict(model, newdata=newD)
	truth = newdf$ILI[260-52+k]
	results[k,1] = fitted
	results[k,2] = truth
	}


	mse = mean((results[,1]-results[,2])^2)
	df.results = data.frame(results)
	df.results = rename(df.results, c("X1"="fitted", "X2"="truth"))
	
# Date formating - for graphs when necessary
	#df.results$unformateddates = fitteddates
	#df.results$Date = strptime(as.character(df.results$unformateddates),"%d/%m/%Y")
	#df.results$Date = format(df.results$Date, "%Y-%m-%d")
	#df.results = df.results[c("Date", "fitted", "truth")]

write.csv(df.results, file = "/Users/f/projects/masters-thesis/src/flu trends/nowcast/Lasso_results_samevars_reestimate.csv")








