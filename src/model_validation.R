nirs.data = read.csv("../NIR.csv", sep=";")

require(leaps)

varHigh3.subsets = regsubsets( N ~ 1+nm1468+nm1472+nm1476+nm1480+nm1484+nm1488+nm1492+nm1496+nm1500+
                                 nm1504+nm1508+nm1512+nm1516+nm1520+nm1524+nm1528+nm1532+nm1536+nm1540+
                                 nm1544+nm1548+nm1552+nm1556+nm1560+nm1564+nm1568+nm1572+nm1576+nm1580+
                                 nm1996+nm2128+nm2132+nm2136+nm2140+nm2144+nm2148+nm2152+nm2156+nm2160+
                                 nm2164+nm2168+nm2172+nm2176+nm2180+nm2184+nm2188+nm2192+nm2196+nm2200+
                                 nm2204+nm2208+nm2212+nm2216+nm2220+nm2224+nm2228+nm2232+nm2236+nm2240+
                                 nm2244+nm2248+nm2252+nm2256+nm2260+nm2264+nm2268+nm2272+nm2276+nm2280+
                                 nm2284+nm2288+nm2292+nm2296+nm2372+nm2376+nm2380+nm2384+nm2388+nm2392+
                                 nm2396+nm2400+nm2404+nm2408+nm2412+nm2416+nm2420+nm2424+nm2428+nm2432+
                                 nm2436+nm2440+nm2444+nm2448+nm2452+nm2456+nm2460+nm2464+nm2468+nm2472+
                                 nm2476+nm2480+nm2484+nm2488+nm2492+nm2496+nm2500+nm2504+nm2508+nm2512+
                                 nm2516+nm2520+nm2524+nm2528+nm2532+nm2536+nm2540+nm2544+nm2548+nm2552+
                                 nm2556+nm2560+nm2564+nm2568+nm2572+nm2576+nm2580+nm2584+nm2588+nm2592+
                                 nm2596+nm2600+nm2604+nm2608+nm2612+nm2616+nm2620+nm2624+nm2628+nm2632+
                                 nm2636+nm2640+nm2644+nm2648+nm2652+nm2656+nm2660+nm2664, nirs.data, really.big=T, nvmax=148, method="backward")

subsets = varHigh3.subsets
optModelId <<- which.min(summary(subsets)$cp)


# #########################################
# Evaluation of model by R^2
# #########################################

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

### can be added to modelselection from here

getR2 = function(data.origin, subsetSize, modelId, TIMES){
  y_exp = simulate(subsets, modelId, data.origin, 533, TIMES)
  y_mean = sum(data.origin[,2])/dim(data.origin)[1]
  y_obs = c(data.origin[,2])
  var_exp = y_exp - y_mean
  var_true = y_obs - y_mean
  R2 = (t(var_exp) %*% var_exp) / (t(var_true) %*% var_true)
  return(R2)
}

R2 = getR2(nirs.data, 533, optModelId, 1000)

####

