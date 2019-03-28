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

nirs.data = read.csv("NIR.csv", sep=";")
head(nirs.data)

# nir data only
nir.data = nirs.data[, 4:322]

# wavelengths
x=seq(1400, 2664, 4)


plot(seq(1400, 2672, 4), nir.data[1,], type='l', col=1, ylim=c(0.1,0.7))
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
subsets = regsubsets( N ~ 1+nm2144+nm2148+nm2152+nm2156+nm2160+nm2164+nm2168+nm2172+nm2176+nm2180+nm2216+nm2220+nm2408+nm2412+nm2416+nm2424+nm2428+nm2476+nm2480+nm2504+nm2508+nm2512+nm2552+nm2556+nm2560+nm2564+nm2568+nm2572+nm2576+nm2580+nm2584+nm2588+nm2592+nm2596+nm2600+nm2656+nm2660+nm2664, nirs_filtered.data, really.big=F, nvmax=39, method="seqrep")
optModelId = which.min(summary(subsets)$cp)
coeffs = as.vector(coef(subsets, optModelId))
optModelId

cp.ground.model = summary(subsets)$cp[20]

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

# Maximales Modell
lm.max =lm(N~1+nm2144+nm2148+nm2152+nm2156+nm2160+nm2164+
             nm2168+nm2172+nm2176+nm2180+nm2216+nm2220+nm2408+
             nm2412+nm2416+nm2424+nm2428+nm2476+nm2480+nm2504+
             nm2508+nm2512+nm2552+nm2556+nm2560+nm2564+nm2568+
             nm2572+nm2576+nm2580+nm2584+nm2588+nm2592+nm2596+
             nm2600+nm2656+nm2660+nm2664, data = nirs_filtered.data)
rss.max = sum(residuals(lm.max)^2)
sigma2.tilde.full = rss.max/(533-40)

spse = (rss.max+ 2*sigma2.tilde.full*40)
spse.2 = cp.ground.model * sigma2.tilde.full
