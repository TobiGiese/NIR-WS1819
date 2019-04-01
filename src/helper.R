# ###############################
# FUNCTIONS: helper
# ###############################
buildFormula = function(features, colName, respVar) {
  featString = paste(features[,colName], collapse="+")  
  featString = paste(1, featString, sep="+")
  formulaString = paste(respVar, featString, sep = " ~ ")
  
  form = as.formula(formulaString)
  
  return (form)
}

addSelectedFeatureToPlot = function(val, color) {
  abline(v=val, col=color)
}


getModelFromSubsets = function(subsets, modelId, responsevar, data) {
  X <- summary(subsets)$which
  xvars <- dimnames(X)[[2]][-1]
  id <- X[modelId,]
  form <- reformulate(xvars[which(id[-1])], responsevar, id[1])
  lm <- lm(form, data)
  
  return(lm)
}