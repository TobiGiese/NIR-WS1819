# ###############################
# FUNCTIONS: spse & rss
# ###############################

calculateRss = function(model) {
  rss = sum(residuals(model)^2)
}

calcauletSigma2TildeFull = function(fullModel, fullModelSize, sampleSize) {
  sigma2.tilde.full = calculateRss(fullModel) / (sampleSize - fullModelSize)
  
  return(sigma2.tilde.full)
}

calculateTrueSpse = function(model, modelSize, fullModel, sampleSize, fullModelSize) {
  sigma2.tilde.full = calcauletSigma2TildeFull(fullModel, fullModelSize, sampleSize)
  spse = calculateRss(model) + 2*sigma2.tilde.full*modelSize
  
  return(spse)
}

calculateEstimatedSpse = function(cp, fullModel, fullModelSize, sampleSize) {
  sigma_2_tilde_full = calcauletSigma2TildeFull(fullModel, fullModelSize, sampleSize)
  spse = cp * sigma_2_tilde_full + sampleSize * sigma_2_tilde_full
  
  return(spse)
}
