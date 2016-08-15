getwd()

experiment = 6

#sink(paste("C:/Users/Fiona.Sloof/Desktop/Thesis/Feb09/Results/Lasso/Lasso_", experiment,".txt", sep = ''))

cRep=10

for (u in 1)
{
case = u

setwd(paste("/Users/f/projects/masters-thesis/data/simulated data/Case", case, sep = ''))
library(glmnet)
set.seed(-1)


#simulation design - open first spreadsheet to create results matrix
data = read.csv("repnum1.csv")
#delete x variables we don't want to include
delete = names(data) %in% c("X","Constant")
#put data into matrix form
newdata = data[!delete]
nregressors = ncol(newdata)

results = matrix(, nrow=nregressors, ncol =cRep)
print(paste("This is LASSO case number", case))
print(paste("There are", cRep, "repetitions in this simulation"))

truepar = read.csv(paste("/Users/f/projects/masters-thesis/data/simulated data/Case", case, "/TrueParValues.csv", sep = ''), stringsAsFactors = FALSE)
ptm <- proc.time()

#create matrix for calculating MSEs

MSEresults <- matrix(nrow = nrow(truepar), ncol = cRep)
CMSEresults <- matrix(nrow = nrow(truepar), ncol = cRep)

#beginning loop for cRep i
for (i in 1:cRep)
{
	csvi=paste("repnum",i,".csv", sep = '')
	data = read.csv(paste("repnum",i,".csv", sep = ''))
	#delete x variables we don't want to include
	delete = names(data) %in% c("X","Constant")
	#put data into matrix form
	newdata = data[!delete]	
	x=model.matrix(y~.,newdata)[,-1]
	y=data$y
	#fit model
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
	#print(lasso.coef)
	
	#create indicator matrix
	Indicators=mat.or.vec(nrow(lasso.coef),1)
	for (j in 1:nrow(lasso.coef))
	{
	if(lasso.coef[j]>0 | lasso.coef[j]<0)
		{
		Indicators[j]=1
		}
	else
		{
		Indicators[j]=0
		}
	}
	arrayIndicators = matrix(Indicators)
	#print(arrayIndicators)
	results[,i]<-arrayIndicators
	
	
	#create df for MSE calculations
	MSE <- data.frame(truepar, lasso.coef)
	MSE$SquaredErrors <- round((MSE$lasso.coef - MSE$TrueParameterValues)^2, digits = 4)
	MSE$ConSquaredErrors <- 0
	for(j in 1:nrow(lasso.coef))
	{
	if(MSE$lasso.coef[j] > 0 | MSE$lasso.coef[j] <0)
	{
	MSE$ConSquaredErrors[j] <- MSE$SquaredErrors[j]	
	}
	}
	
	MSEresults[,i]<-MSE$SquaredErrors
	CMSEresults[,i] <- MSE$ConSquaredErrors
	
	
	

}

#statistics
RetRates = rowMeans(results)
UMSE = rowMeans(MSEresults)
CMSEinc <- rowSums(CMSEresults != 0)
CMSEsum <- rowSums(CMSEresults)
CMSE = CMSEsum/CMSEinc

Root_UMSE = round(sqrt(UMSE), digits = 4)
Root_CMSE = round(sqrt(CMSE), digits = 4)


 #Calculate total MSE of retained irrelevant variables --- AGAIN MISSING SQUARE ROOT
 RetIrrVarMSE <- matrix(0, nrow(MSEresults), 1)
 for (j in 1:nrow(results))                        
 {
   
  if (MSE$TrueParameterValues[j]== 0)
  {
   
   if (sum(results[j,]) > 0)
   {
   RetIrrVarMSE[j] = sum(MSEresults[j,])/sum(results[j,])
   }
   else 
   {
    RetIrrVarMSE[j] = 0
   }
  }
    if   (MSE$TrueParameterValues[j] != 0)
  {
      RetIrrVarMSE[j] = 0
  }
 
 }
 
RIV_MSE <- sum(RetIrrVarMSE)





#print(RetRates)
regressors = names(newdata)
df = data.frame(regressors,RetRates,truepar[2], Root_UMSE, Root_CMSE)
#drop reporting of constant
df <- df[-c(1),]
#drop y since not including in gauge/potency calculation
df2 <- df[!(df$regressors=="Ly"),]




gaugematrix <- subset(df2, TrueParameterValues == 0)
potencymatrix <- subset(df2, TrueParameterValues > 0 | TrueParameterValues < 0)

#print(gaugematrix)
#print(potencymatrix)

gauge = mean(gaugematrix$RetRates)
potency = mean(potencymatrix$RetRates)

summary <- df[c(-3)]
print("The retention rates are:")
print(df)
print("The gauge is:")
print(gauge)
print("The potency is:")
print(potency)
print("The total RMSE of retained irrelevant variables is")
print(RIV_MSE)
}

time = proc.time()-ptm
print(time)





