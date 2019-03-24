getDiff = function(row, n) {
  x=c()
  for (i in 1:317) {
    x[i] = abs(row[i]-row[i+1])
  }
  
  return(x)
}

getDiffMatrix = function(M, rows, cols) {
  x1 = getDiff(t(M)[,1], cols)
  x2 = getDiff(t(M)[,2], cols)
  x=rbind(x1,x2)
  for (i in 3:rows) {
    x = rbind(x, getDiff(t(M)[,i], cols))
  }
  
  return(x)
}

nirs.data = read.csv("../NIR.csv", sep=";")
head(nirs.data)

# nir data only
nir.data = nirs.data[, 4:322]

# wavelengths
x=seq(1400, 2664, 4)

plot(seq(1400, 2672, 4), nir.data[1,], type='l', col=1, ylim=c(0.3,0.7))
lines(seq(1400, 2672, 4), nir.data[2,], type='l', col=2, ylim=c(0.3,0.7))
lines(seq(1400, 2672, 4), nir.data[3,], type='l', col=3, ylim=c(0.3,0.7))
lines(seq(1400, 2672, 4), nir.data[4,], type='l', col=4, ylim=c(0.3,0.7))
lines(seq(1400, 2672, 4), nir.data[5,], type='l', col=5, ylim=c(0.3,0.7))

# plot first deviation for first 6 records
plot(x, getDiff(t(nir.data)[,1], 318), type='l', col=1)
lines(x, getDiff(t(nir.data)[,2], 318), type='l', col=2)
lines(x, getDiff(t(nir.data)[,3], 318), type='l', col=3)
lines(x, getDiff(t(nir.data)[,4], 318), type='l', col=4)
lines(x, getDiff(t(nir.data)[,5], 318), type='l', col=5)
lines(x, getDiff(t(nir.data)[,6], 318), type='l', col=6)

# calculate first deviation for the whole dataset
M = getDiffMatrix(nir.data, 533, 318)
# first deviation average
Mavg = colMeans(M)
plot(x, Mavg, type='l', col=1)
# threshold 
THRESHOLD = 0.0015
abline(h=THRESHOLD, col=2)


# filter dataset for colunms having an average value (1. deviation) above threshold
df = data.frame(x, colnames(nirs.data[, 4:320]), Mavg)
colnames(df) = c("x", "nm", "y")
filteredAvg1stDev = subset(df, y >= THRESHOLD, select=c("x", "nm", "y"));
nrow(filteredAvg1stDev)
plot(filteredAvg1stDev$x, filteredAvg1stDev$y, type='p', col=1)

nirs_filtered.data = subset(nirs.data, select=c(c("SOC", "N", "pH"), as.character(filteredAvg1stDev[,"nm"])))

### Mallows CP
require(leaps)

# Select model
subsets = regsubsets( N ~ nm2144+nm2148+nm2152+nm2156+nm2160+nm2164+nm2168+nm2172+nm2176+nm2180+nm2216+nm2220+nm2408+nm2412+nm2416+nm2424+nm2428+nm2476+nm2480+nm2504+nm2508+nm2512+nm2552+nm2556+nm2560+nm2564+nm2568+nm2572+nm2576+nm2580+nm2584+nm2588+nm2592+nm2596+nm2600+nm2656+nm2660+nm2664, nirs_filtered.data, really.big=F, nvmax=38)
optModelId = which.min(summary(subsets)$cp)
coeffs = coef(subsets, optModelId)

# ##########################################################
# ## SIMULATION
# ##########################################################
getColumnMean = function(colname, data) {
  return(mean(data[,colname]))
}

getColumnSD = function(colname, data) {
  return(sd(data[,colname]))
}


simulateColumn = function(num, colname, data.origin, TIMES) {
  data = data.frame(matrix(nrow=num, ncol=TIMES));
  for (i in 1:TIMES) {
    data[i] = rnorm(num, mean=getColumnMean(colname, data.origin), sd = getColumnSD(colname, data.origin))  
  }
  
  dataAVG = rowMeans(data)
  
  return(dataAVG)
}

simulateDataset = function(numRows, coeffs, data.origin, TIMES) {
  sim.data = data.frame(matrix(nrow=numRows));
  
  # simulate nir data
  for (feature in colnames(data.origin)[-1:-3]) {
    sim.data[,feature] = simulateColumn(numRows, feature, data.origin, TIMES)
  }
  
  # use estimated coefficients to calculate value for N
  sim.data$N = as.vector(coeffs)[1];  
  for (feature in names(coeffs[-1])) {
    sim.data$N = sim.data$N + coeffs[feature] + sim.data[,feature]
  }
  
  return(sim.data[,-1])
}


set.seed(13)
simdata = simulateDataset(30, coeffs, nirs.data, 100);
simdata
#plot(seq(1400, 2672, 4), simdata[1,], type='l', col=1, ylim=c(0.3,0.7))
#lines(seq(1400, 2672, 4), simdata[2,], type='l', col=2, ylim=c(0.3,0.7))
#lines(seq(1400, 2672, 4), simdata[3,], type='l', col=3, ylim=c(0.3,0.7))
#lines(seq(1400, 2672, 4), simdata[4,], type='l', col=4, ylim=c(0.3,0.7))
#lines(seq(1400, 2672, 4), simdata[5,], type='l', col=5, ylim=c(0.3,0.7))

