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
abline(h=0.0015, col=2)


# filter dataset for colunms having an average value (1. deviation) above threshold
df = data.frame(x, colnames(nirs.data[, 4:320]), Mavg)
colnames(df) = c("x", "nm", "y")
filteredAvg1stDev = subset(df, y >= 0.0015, select=c("x", "nm", "y"));
nrow(filteredAvg1stDev)
plot(filteredAvg1stDev$x, filteredAvg1stDev$y, type='p', col=1)

nirs_filtered.data = subset(nirs.data, select=c(c("SOC", "N", "pH"), as.character(filteredAvg1stDev[,"nm"])))

### Mallows CP
require(leaps)

#subsets = regsubsets( N ~ nm1512+nm1520+nm1524+nm1528+nm1532+nm1536+nm1548+nm1552+nm2140+nm2144+nm2148+nm2152+nm2156+nm2160+nm2164+nm2168+nm2172+nm2176+nm2180+nm2184+nm2188+nm2192+nm2196+nm2200+nm2204+nm2208+nm2212+nm2216+nm2220+nm2224+nm2228+nm2232+nm2236+nm2240+nm2244+nm2248+nm2252+nm2256+nm2260+nm2264+nm2268+nm2272+nm2384+nm2388+nm2396+nm2400+nm2404+nm2408+nm2412+nm2416+nm2420+nm2424+nm2428+nm2432+nm2448+nm2452+nm2456+nm2464+nm2468+nm2472+nm2476+nm2480+nm2484+nm2488+nm2492+nm2496+nm2500+nm2504+nm2508+nm2512+nm2528+nm2532+nm2536+nm2540+nm2544+nm2548+nm2552+nm2556+nm2560+nm2564+nm2568+nm2572+nm2576+nm2580+nm2584+nm2588+nm2592+nm2596+nm2600+nm2604+nm2608+nm2636+nm2640+nm2652+nm2656+nm2660+nm2664, nirs_filtered.data, really.big=T, nvmax=97)
subsets = regsubsets( N ~ nm2144+nm2148+nm2152+nm2156+nm2160+nm2164+nm2168+nm2172+nm2176+nm2180+nm2216+nm2220+nm2408+nm2412+nm2416+nm2424+nm2428+nm2476+nm2480+nm2504+nm2508+nm2512+nm2552+nm2556+nm2560+nm2564+nm2568+nm2572+nm2576+nm2580+nm2584+nm2588+nm2592+nm2596+nm2600+nm2656+nm2660+nm2664, nirs_filtered.data, really.big=F, nvmax=38)
summary(subsets)$cp
