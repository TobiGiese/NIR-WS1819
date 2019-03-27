getDiff = function(row, n) {
  x=c()
  for (i in 1:(length(row)-1)) {
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
#nirs.dataNew = nirs.data[,seq(2, ncol(nirs.data), 2)]
nirs.dataNew = nirs.data[,1:3]
nirs.dataNew = cbind(nirs.dataNew, nirs.data[,sample(150:ncol(nirs.data)-3, 39)])
nirs.dataNew$SOC = nirs.data$SOC
nirs.dataNew$pH = nirs.data$pH
nirs.data = nirs.dataNew

# nir data only
nir.data = nirs.data[, 2:(ncol(nirs.data)-2)]
ncol(nir.data)

# wavelengths
x=seq(1400, 2664, 4)
length(x)
length(nir.data[1,])

plot(seq(1400, 2672, 8), nir.data[1,], type='l', col=1, ylim=c(0.3,0.7))
lines(seq(1400, 2672, 8), nir.data[2,], type='l', col=2, ylim=c(0.3,0.7))
lines(seq(1400, 2672, 8), nir.data[3,], type='l', col=3, ylim=c(0.3,0.7))
lines(seq(1400, 2672, 8), nir.data[4,], type='l', col=4, ylim=c(0.3,0.7))
lines(seq(1400, 2672, 8), nir.data[5,], type='l', col=5, ylim=c(0.3,0.7))

# plot first deviation for first 6 records
length(x)
foo = getDiff(t(nir.data)[,1], 159)
length(foo)
plot(x, getDiff(t(nir.data)[,1], 159), type='l', col=1, ylim = c(0, 0.010))
lines(x, getDiff(t(nir.data)[,2], 159), type='l', col=2)
lines(x, getDiff(t(nir.data)[,3], 159), type='l', col=3)
lines(x, getDiff(t(nir.data)[,4], 159), type='l', col=4)
lines(x, getDiff(t(nir.data)[,5], 159), type='l', col=5)
lines(x, getDiff(t(nir.data)[,6], 159), type='l', col=6)

# calculate first deviation for the whole dataset
M = getDiffMatrix(nir.data, 533, 159)
# first deviation average
Mavg = colMeans(M)
plot(x, Mavg, type='l', col=1)
# threshold 
THRESHOLD = 0.0015
abline(h=THRESHOLD, col=2)


# filter dataset for colunms having an average value (1. deviation) above threshold
df = data.frame(x, colnames(nirs.data[, 2:(ncol(nirs.data)-3)]), Mavg)
colnames(df) = c("x", "nm", "y")
filteredAvg1stDev = subset(df, y >= THRESHOLD, select=c("x", "nm", "y"));
nrow(filteredAvg1stDev)
plot(filteredAvg1stDev$x, filteredAvg1stDev$y, type='p', col=1)

nirs_filtered.data = subset(nirs.data, select=c(c("SOC", "N", "pH"), as.character(filteredAvg1stDev[,"nm"])))

nirs_filtered.data = nirs.dataNew;
### Mallows CP
require(leaps)

# Select model
#subsets = regsubsets( N ~ 1+nm1512+nm1520+nm1528+nm1536+nm1544+nm1552+nm2144+nm2152+nm2160+nm2168+nm2176+nm2184+nm2192+nm2200+nm2208+nm2216+nm2224+nm2232+nm2240+nm2248+nm2256+nm2264+nm2272+nm2400+nm2408+nm2416+nm2424+nm2448+nm2464+nm2472+nm2480+nm2496+nm2504+nm2512+nm2528+nm2536+nm2544+nm2552+nm2560+nm2568+nm2576+nm2584+nm2592+nm2600+nm2648+nm2656+nm2664, nirs_filtered.data, really.big=F, nvmax=39, method="exhaustive")
subsets = regsubsets( N ~ 1+nm2076+nm2304+nm2428+nm2352+nm2468+nm2116+nm2244+nm2280+nm2176+nm2368+nm2240+nm2548+nm2532+nm2640+nm2032+nm2052+nm2136+nm2628+nm2268+nm2204+nm1976+nm2252+nm2500+nm2484+nm2284+nm2292+nm2576+nm2288+nm2060+nm2448+nm2480+nm2092+nm2524+nm2248+nm2208+nm2256+nm2308+nm2440+nm2488, nirs_filtered.data, really.big=F, nvmax=39)
optModelId = which.min(summary(subsets)$cp)
summary(subsets)$cp
optModelId
coeffs = as.vector(coef(subsets, optModelId))



# ##########################################################
# ## SIMULATION
# ##########################################################
set.seed(13)
getDesignmatrix = function(subsets, modelId, responsevar, data) {
  X <- summary(subsets)$which
  xvars <- dimnames(X)[[2]][-1]
  id <- X[modelId,]
  form <- reformulate(xvars[which(id[-1])], responsevar, id[1])
  mod.nir <- lm(form, data)
  designMatrix <- model.matrix(mod.nir)
  
  return(designMatrix)
}

getMeans = function(designMatrix, subsets, optModelId) {
  coeffs = as.vector(coef(subsets, optModelId))
  means = designMatrix %*% coeffs
  
  return(means)
}

getSDs = function(designMatrix, means, data.origin) {
  residuals = as.vector(data.origin$N - means)
  sd = sqrt(as.numeric(t(residuals) %*% residuals / (dim(data.origin)[1] - dim(designMatrix)[2])))
  
  return(sd)
}

# Generiert num Pseudodaten für Zielgröße N
# num: Anzahl zu erzeugender Pseudodaten
# TIMES: Anzahl wiederholungen
simulate = function(subsets, modelId, data.origin, num, TIMES) {
  data = data.frame(matrix(nrow = num, ncol=TIMES));
  
  designmatrix = getDesignmatrix(subsets, modelId, "N", data.origin)
  means = getMeans(designmatrix, subsets, modelId)
  sd = getSDs(designmatrix, means, data.origin)
  
  for (i in 1:TIMES) {
    data[i] = rnorm(num, mean=means, sd = sd)  
  }
  
  dataAVG = rowMeans(data)
  return(dataAVG)
}

# simulate 
simdata = nirs_filtered.data;
simdata$N = simulate(subsets, optModelId, nirs.data, nrow(nirs.data), 1000)
plot(nirs.data[,4], nirs.data[,2], type="p", col=1, ylim=c(-0.1, 0.75))
points(nirs.data[,4], simdata$N, type="p", col=2, ylim=c(-0.1, 0.75))

simsets533 = regsubsets( N ~ 1+nm2144+nm2148+nm2152+nm2156+nm2160+nm2164+nm2168+nm2172+nm2176+nm2180+nm2216+nm2220+nm2408+nm2412+nm2416+nm2424+nm2428+nm2476+nm2480+nm2504+nm2508+nm2512+nm2552+nm2556+nm2560+nm2564+nm2568+nm2572+nm2576+nm2580+nm2584+nm2588+nm2592+nm2596+nm2600+nm2656+nm2660+nm2664, simdata, really.big=F, nvmax=39)
optModelId533 = which.min(summary(simsets533)$cp)
print(paste("========= Num Rows: 533 ========="))
print(paste("Optimal Model (ID): ", optModelId533))
print("CP Value")
print(summary(simsets533)$cp)

# Shrink designmatrix
getDataSubset = function(data.origin, numRows) {
  if (numRows < nrow(data.origin)) {
    return(nirs.data[sample(nrow(data.origin), numRows), ])
  }
  
  return(nirs.data)
}

# rum simulation and second model selection for a subset of nirs.data with increasing size
# TODO: repeat simulation and model selection x times
for (numRows in seq(50, 500, 50)) {
  print(paste("========= Num Rows: ", numRows, " ========="))
  nirs.data100 = getDataSubset(nirs.data, numRows)
  simdata100 = subset(nirs.data100, select=c(c("SOC", "N", "pH"), as.character(filteredAvg1stDev[,"nm"])))
  simdata100$N = simulate(subsets, optModelId, nirs.data100, numRows, 1000)
  plot(nirs.data[,4], nirs.data[,2], type="p", col=1, ylim=c(-0.1, 0.75))
  points(nirs.data100[,4], simdata100$N, type="p", col=2, ylim=c(-0.1, 0.75))
  
  simsets100 = regsubsets( N ~ 1+nm2144+nm2148+nm2152+nm2156+nm2160+nm2164+nm2168+nm2172+nm2176+nm2180+nm2216+nm2220+nm2408+nm2412+nm2416+nm2424+nm2428+nm2476+nm2480+nm2504+nm2508+nm2512+nm2552+nm2556+nm2560+nm2564+nm2568+nm2572+nm2576+nm2580+nm2584+nm2588+nm2592+nm2596+nm2600+nm2656+nm2660+nm2664, simdata100, really.big=F, nvmax=39)
  optModelId100 = which.min(summary(simsets100)$cp)
  print(paste("Optimal Model (ID): ", optModelId100))
  print("CP Value")
  print(summary(simsets100)$cp)
}
