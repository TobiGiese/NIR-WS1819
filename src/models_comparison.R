getModelFromSubsets = function(subsets, modelId, responsevar, data) {
  X <- summary(subsets)$which
  xvars <- dimnames(X)[[2]][-1]
  id <- X[modelId,]
  form <- reformulate(xvars[which(id[-1])], responsevar, id[1])
  lm <- lm(form, data)
  
  return(lm)
}

calculateRss = function(model) {
  rss = sum(residuals(model)^2)
}

calcauletSigma2TildeFull = function(fullModel, fullModelSize, sampleSize) {
  sigma2.tilde.full = calculateRss(fullModel) / (sampleSize - fullModelSize)
  
  return(sigma2.tilde.full)
}

calculateTrueSpse = function(model, modelSize, sampleSize) {
  rss = calculateRss(model)
  sigma2.tilde = rss/(sampleSize - modelSize)
  spse = (sampleSize+modelSize)*sigma2.tilde
  
  return(spse)
}

calculateTrueSpse2 = function(model, modelSize, fullModel, sampleSize, fullModelSize) {
  sigma2.tilde.full = calcauletSigma2TildeFull(fullModel, fullModelSize, sampleSize)
  spse = calculateRss(model) + 2*sigma2.tilde.full*modelSize
  
  return(spse)
}

calculateEstimatedSpse = function(cp, fullModel, fullModelSize, sampleSize) {
  sigma_2_tilde_full = calcauletSigma2TildeFull(fullModel, fullModelSize, sampleSize)
  spse = cp * sigma_2_tilde_full + sampleSize * sigma_2_tilde_full#
  
  return(spse)
}


nirs.data = read.csv("../NIR.csv", sep=";")

# model with highest 1. deviation
devHigh.subsets = regsubsets( N ~ 1+nm2144+nm2148+nm2152+I(nm2156^2)+I(nm2160^2)+I(nm2164^2)+I(nm2168^2)+I(nm2172^2)+nm2176+nm2180+nm2216+nm2220+nm2408+nm2412+nm2416+nm2424+nm2428+nm2476+nm2480+nm2504+nm2508+nm2512+nm2552+nm2556+nm2560+nm2564+nm2568+nm2572+nm2576+nm2580+nm2584+nm2588+nm2592+nm2596+nm2600+nm2656+nm2660+nm2664, nirs.data, really.big=F, nvmax=39)
devHigh.lm.full = lm(N ~ 1+nm2144+nm2148+nm2152+I(nm2156^2)+I(nm2160^2)+I(nm2164^2)+I(nm2168^2)+I(nm2172^2)+nm2176+nm2180+nm2216+nm2220+nm2408+nm2412+nm2416+nm2424+nm2428+nm2476+nm2480+nm2504+nm2508+nm2512+nm2552+nm2556+nm2560+nm2564+nm2568+nm2572+nm2576+nm2580+nm2584+nm2588+nm2592+nm2596+nm2600+nm2656+nm2660+nm2664, data=nirs.data)
devHigh.cp = min(summary(devHigh.subsets)$cp)
devHigh.optModelId = which.min(summary(devHigh.subsets)$cp)
devHigh.lm.opt = getModelFromSubsets(devHigh.subsets, devHigh.optModelId, "N", nirs.data)

# model with lowest 1. deviation
devLow.subsets = regsubsets( N ~ 1+nm1440+nm1452+nm1456+nm1460+nm1464+nm1472+nm1476+nm1860+nm1864+nm1868+nm1872+nm1876+nm1880+nm1884+nm1888+nm1892+nm1896+nm1900+nm1904+nm1908+nm1912+nm1916+nm1920+nm1924+nm1928+nm1932+nm1936+nm1940+nm1944+nm1948+nm1952+nm1956+nm1960+nm1964+nm1968+nm1972+nm1976+nm1980+nm1984+nm1988+nm2008+nm2012+nm2016+nm2020+nm2024+nm2028+nm2040+nm2044+nm2048+nm2052+nm2056+nm2060+nm2088+nm2092+nm2108+nm2112+nm2116+nm2124+nm2328+nm2332+nm2336+nm2340+nm2344+nm2348+nm2352+nm2356+nm2360+nm2364+nm2368, nirs.data, really.big=T, nvmax=70, method="backward")
devLow.lm.full = lm(N ~ 1+nm1440+nm1452+nm1456+nm1460+nm1464+nm1472+nm1476+nm1860+nm1864+nm1868+nm1872+nm1876+nm1880+nm1884+nm1888+nm1892+nm1896+nm1900+nm1904+nm1908+nm1912+nm1916+nm1920+nm1924+nm1928+nm1932+nm1936+nm1940+nm1944+nm1948+nm1952+nm1956+nm1960+nm1964+nm1968+nm1972+nm1976+nm1980+nm1984+nm1988+nm2008+nm2012+nm2016+nm2020+nm2024+nm2028+nm2040+nm2044+nm2048+nm2052+nm2056+nm2060+nm2088+nm2092+nm2108+nm2112+nm2116+nm2124+nm2328+nm2332+nm2336+nm2340+nm2344+nm2348+nm2352+nm2356+nm2360+nm2364+nm2368, data=nirs.data)
devLow.cp = min(summary(devLow.subsets)$cp)
devLow.optModelId = which.min(summary(devLow.subsets)$cp)
devLow.lm.opt = getModelFromSubsets(devLow.subsets, devLow.optModelId, "N", nirs.data)

# model with highest variablilty, based on normalized distance between max and min (first dev)
varHigh.subsets = regsubsets( N ~ 1+nm1456+nm1460+nm1464+nm1472+nm1476+nm1908+nm1912+nm1932+nm1936+nm1960+nm1964+nm1968+nm1972+nm1976+nm1980+nm1984+nm1996+nm2020+nm2112+nm2132+nm2348+nm2360+nm2364+nm2392+nm2396+nm2436+nm2456+nm2460+nm2492+nm2496+nm2520+nm2524+nm2616+nm2624+nm2628+nm2632, nirs.data, really.big=T, nvmax=37)
varHigh.lm.full = lm(N ~ 1+nm1456+nm1460+nm1464+nm1472+nm1476+nm1908+nm1912+nm1932+nm1936+nm1960+nm1964+nm1968+nm1972+nm1976+nm1980+nm1984+nm1996+nm2020+nm2112+nm2132+nm2348+nm2360+nm2364+nm2392+nm2396+nm2436+nm2456+nm2460+nm2492+nm2496+nm2520+nm2524+nm2616+nm2624+nm2628+nm2632, data=nirs.data)
varHigh.cp = min(summary(varHigh.subsets)$cp)
varHigh.optModelId = which.min(summary(varHigh.subsets)$cp)
varHigh.lm.opt = getModelFromSubsets(varHigh.subsets, varHigh.optModelId, "N", nirs.data)

#  model with highest 1. deviation (without abs)
devHigh2.subsets = regsubsets( N ~ 1+nm1444+nm1548+nm1552+nm1556+nm1560+nm1564+nm1568+nm1572+nm1592+nm2188+nm2192+nm2196+nm2200+nm2204+nm2208+nm2212+nm2216+nm2220+nm2224+nm2228+nm2232+nm2236+nm2240+nm2244+nm2248+nm2252+nm2256+nm2260+nm2264+nm2268+nm2272+nm2276+nm2280+nm2284+nm2424+nm2428+nm2432+nm2636, nirs.data, really.big=T, nvmax=39)
devHigh2.lm.full = lm(N ~ 1+nm1444+nm1548+nm1552+nm1556+nm1560+nm1564+nm1568+nm1572+nm1592+nm2188+nm2192+nm2196+nm2200+nm2204+nm2208+nm2212+nm2216+nm2220+nm2224+nm2228+nm2232+nm2236+nm2240+nm2244+nm2248+nm2252+nm2256+nm2260+nm2264+nm2268+nm2272+nm2276+nm2280+nm2284+nm2424+nm2428+nm2432+nm2636, data=nirs.data)
devHigh2.cp = min(summary(devHigh2.subsets)$cp)
devHigh2.optModelId = which.min(summary(devHigh2.subsets)$cp)
devHigh2.lm.opt = getModelFromSubsets(devHigh2.subsets, devHigh2.optModelId, "N", nirs.data)

#  model with highest variability, based on normalized distance between max and  min (without abs)
varHigh2.subsets = regsubsets( N ~ 1+nm1456+nm1460+nm1464+nm1472+nm1476+nm1484+nm1936+nm1956+nm1960+nm1964+nm1968+nm1972+nm1976+nm1980+nm1984+nm2112+nm2124+nm2136+nm2184+nm2188+nm2360+nm2436+nm2456+nm2492+nm2496+nm2520+nm2524+nm2616+nm2620+nm2624+nm2628+nm2632+nm2640+nm2644+nm2648, nirs.data, really.big=T, nvmax=36)
varHigh2.lm.full = lm(N ~ 1+nm1456+nm1460+nm1464+nm1472+nm1476+nm1484+nm1936+nm1956+nm1960+nm1964+nm1968+nm1972+nm1976+nm1980+nm1984+nm2112+nm2124+nm2136+nm2184+nm2188+nm2360+nm2436+nm2456+nm2492+nm2496+nm2520+nm2524+nm2616+nm2620+nm2624+nm2628+nm2632+nm2640+nm2644+nm2648, data=nirs.data)
varHigh2.cp = min(summary(varHigh2.subsets)$cp)
varHigh2.optModelId = which.min(summary(varHigh2.subsets)$cp)
varHigh2.lm.opt = getModelFromSubsets(varHigh2.subsets, varHigh2.optModelId, "N", nirs.data)

# model with highest variability, based on distance between max and  min (without abs)
varHigh3.subsets = regsubsets( N ~ 1+nm1500+nm1504+nm1508+nm1512+nm1516+nm1520+nm1524+nm1528+nm1532+nm1536+nm1540+nm1548+nm1552+nm1556+nm2136+nm2140+nm2144+nm2148+nm2152+nm2156+nm2160+nm2164+nm2168+nm2172+nm2176+nm2180+nm2184+nm2188+nm2192+nm2196+nm2200+nm2204+nm2220+nm2380+nm2384+nm2388+nm2392+nm2396+nm2400+nm2404+nm2408+nm2412+nm2416+nm2420+nm2424+nm2428+nm2432+nm2436+nm2440+nm2444+nm2448+nm2452+nm2456+nm2460+nm2464+nm2468+nm2472+nm2476+nm2480+nm2484+nm2488+nm2492+nm2496+nm2500+nm2504+nm2508+nm2512+nm2516+nm2520+nm2524+nm2528+nm2532+nm2552+nm2556+nm2560+nm2564+nm2568+nm2572+nm2576+nm2580+nm2584+nm2588+nm2592+nm2596+nm2600+nm2604+nm2608+nm2612+nm2616+nm2620+nm2624+nm2628+nm2632+nm2636+nm2640+nm2644+nm2648+nm2652+nm2656+nm2660+nm2664, nirs.data, really.big=T, nvmax=102, method="backward")
varHigh3.lm.full = lm(N ~ 1+nm1500+nm1504+nm1508+nm1512+nm1516+nm1520+nm1524+nm1528+nm1532+nm1536+nm1540+nm1548+nm1552+nm1556+nm2136+nm2140+nm2144+nm2148+nm2152+nm2156+nm2160+nm2164+nm2168+nm2172+nm2176+nm2180+nm2184+nm2188+nm2192+nm2196+nm2200+nm2204+nm2220+nm2380+nm2384+nm2388+nm2392+nm2396+nm2400+nm2404+nm2408+nm2412+nm2416+nm2420+nm2424+nm2428+nm2432+nm2436+nm2440+nm2444+nm2448+nm2452+nm2456+nm2460+nm2464+nm2468+nm2472+nm2476+nm2480+nm2484+nm2488+nm2492+nm2496+nm2500+nm2504+nm2508+nm2512+nm2516+nm2520+nm2524+nm2528+nm2532+nm2552+nm2556+nm2560+nm2564+nm2568+nm2572+nm2576+nm2580+nm2584+nm2588+nm2592+nm2596+nm2600+nm2604+nm2608+nm2612+nm2616+nm2620+nm2624+nm2628+nm2632+nm2636+nm2640+nm2644+nm2648+nm2652+nm2656+nm2660+nm2664, data=nirs.data)
varHigh3.cp = min(summary(varHigh3.subsets)$cp)
varHigh3.optModelId = which.min(summary(varHigh3.subsets)$cp)
varHigh3.lm.opt = getModelFromSubsets(varHigh3.subsets, varHigh3.optModelId, "N", nirs.data)

# #########################################
# calculate true SPSE
# #########################################

devHigh.spse.true = calculateTrueSpse(devHigh.lm.opt, devHigh.optModelId, 533)
devLow.spse.true = calculateTrueSpse(devLow.lm.opt, devLow.optModelId, 533)
varHigh.spse.true = calculateTrueSpse(varHigh.lm.opt, varHigh.optModelId, 533)
devHigh2.spse.true = calculateTrueSpse(devHigh2.lm.opt, devHigh2.optModelId, 533)
varHigh2.spse.true = calculateTrueSpse(varHigh2.lm.opt, varHigh2.optModelId, 533)
varHigh3.spse.true = calculateTrueSpse(varHigh3.lm.opt, varHigh3.optModelId, 533)

devHigh.spse.true
devLow.spse.true
varHigh.spse.true
devHigh2.spse.true
varHigh2.spse.true
varHigh3.spse.true


# ###################################################
# calculate true SPSE vith RSS ans sigme.tilde.full
# ##################################################

devHigh.spse.true2 = calculateTrueSpse2(devHigh.lm.opt, devHigh.optModelId, devHigh.lm.full, 533, 39)
devLow.spse.true2 = calculateTrueSpse2(devLow.lm.opt, devLow.optModelId, devLow.lm.full, 533, 70)
varHigh.spse.true2 = calculateTrueSpse2(varHigh.lm.opt, varHigh.optModelId, varHigh.lm.full, 533, 37)
devHigh2.spse.true2 = calculateTrueSpse2(devHigh2.lm.opt, devHigh2.optModelId, devHigh2.lm.full, 533, 39)
varHigh2.spse.true2 = calculateTrueSpse2(varHigh2.lm.opt, varHigh2.optModelId, varHigh2.lm.full, 533, 36)
varHigh3.spse.true2 = calculateTrueSpse2(varHigh3.lm.opt, varHigh3.optModelId, varHigh3.lm.full, 533, 37)

devHigh.spse.true2
devLow.spse.true2
varHigh.spse.true2
devHigh2.spse.true2
varHigh2.spse.true2
varHigh3.spse.true2

# ###################################################
# estimate SPSE based on CP value
# ##################################################

devHigh.spse.est.cp = calculateEstimatedSpse(devHigh.cp, devHigh.lm.full, 39, 533)
devLow..spse.est.cp = calculateEstimatedSpse(devLow.cp, devLow.lm.full, 70, 533)
varHigh.spse.est.cp = calculateEstimatedSpse(varHigh.cp, varHigh.lm.full, 37, 533)
devHigh2.spse.est.cp = calculateEstimatedSpse(devHigh2.cp, devHigh2.lm.full, 39, 533)
varHigh2.spse.est.cp = calculateEstimatedSpse(varHigh2.cp, varHigh2.lm.full, 36, 533)
varHigh3.spse.est.cp = calculateEstimatedSpse(varHigh3.cp, varHigh3.lm.full, 37, 533)

devHigh.spse.est.cp 
devLow..spse.est.cp 
varHigh.spse.est.cp 
devHigh2.spse.est.cp
varHigh2.spse.est.cp
varHigh3.spse.est.cp
