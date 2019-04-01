# ################################################
# File: simulation.R
# Scope: functions used for simulation
# ################################################
source("helper.R")

getSDs = function(designMatrix, means, data.origin) {
  residuals = as.vector(data.origin$N - means)
  sd = sqrt(as.numeric(t(residuals) %*% residuals / (dim(data.origin)[1] - dim(designMatrix)[2])))
  
  return(sd)
}

simulate = function(model, data, num, times) {
  simdata = data.frame(matrix(nrow=num, ncol=times));
  
  coeffs = as.vector(coef(model))
  designmatrix = model.matrix(model, data=data)

  means = designmatrix %*% coeffs
  sds = getSDs(designmatrix, means, data)
  for (i in 1:times) {
    simdata[i] = rnorm(num, mean=means, sd = sds)  
  }
  
  return(rowMeans(simdata))
}

getDataSubset = function(data.origin, numRows) {
  if (numRows < nrow(data.origin)) {
    return(nirs.data[sample(nrow(data.origin), numRows), ])
  }
  
  return(nirs.data)
}

simulateOnDataSubset = function(model, data, fmForm, fmFeat, numRows, times) {
  simdata = getDataSubset(data, numRows)
  simdata = subset(simdata, select=c(c("SOC", "N", "pH"), as.character(fmFeat[,"nm"])))
  simdata$N = simulate(model, simdata, numRows, times)
  
  simsets = regsubsets( fmForm, data=simdata, really.big=T, nvmax=nrow(fmFeat)+1, method="backward")
  simOptModelId = which.min(summary(simsets)$cp)
  simOptModelCpValue = min(summary(simsets)$cp)

  simResult = c(simOptModelId, simOptModelCpValue)
  names(simResult) = c("id", "cp")
  
  return(simResult)
}

calculateR2 = function(model, data, num, times){
  y_exp = simulate(model, data, num, times)
  y_mean = sum(data[,2])/dim(data)[1]
  y_obs = c(data[,2])
  var_exp = y_exp - y_mean
  var_true = y_obs - y_mean
  R2 = (t(var_exp) %*% var_exp) / (t(var_true) %*% var_true)
  
  return(R2)
}