# ##################################################
# File: modelselect.R
# Scope: functions used for initial model selection
# ##################################################

getSlope = function(data, col_from, col_to) {
  slope = matrix(nrow = nrow(data), ncol=col_to-col_from)
  for (i in 1:(col_to-col_from)) {
    slope[,i] = data[,col_from+i] - data[,col_from+i-1]
  }
  
  return(slope)
}

# find wavelengths with highest variablilty as the distance of max and min slope
criterionSlopeDist = function(data, col_from, col_to) {
  slope = getSlope(data, col_from, col_to)
  min = apply(slope, 2, min)
  max = apply(slope, 2, max)
  dist = max - min
  
  return (dist)
}

selectFeatures = function(data, feat.all, th, crit) {
  df = data.frame(
    WAVELENGTHS[1:length(WAVELENGTHS)-1], 
    colnames(data[, COLUMN_ID_NM_FROM:(COLUMN_ID_NM_TO-1)]), 
    crit
  )
  colnames(df) = c("x", "nm", "y")
  feat.select = subset(df, y >= th, select=c("x", "nm", "y"));
  
  return (feat.select)
}