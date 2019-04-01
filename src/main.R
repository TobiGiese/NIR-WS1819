# ###############
# DEPENDENCIES
# ##############
require(leaps)
source("helper.R")
source("modelselect.R")
source("simulation.R")
source("spse.R")


# ###############
# MAIN
# ##############
nirs.data = read.csv("../NIR.csv", sep=";")

set.seed(13)

# general 
COLUMN_ID_N = 2
COLUMN_ID_NM_FROM = 4
COLUMN_ID_NM_TO = ncol(nirs.data)
NUM_ROWS = nrow(nirs.data)
WAVELENGTHS = seq(1400, 2672, 4)
RESPONSEVAR = "N"

# modelselection
THRESHOLD = 0.001

# simulation
SAMPLE_SIZES = c(150, 200, 250, 300, 350, 400, 450, 500, 533)
SIM_ITERATIONS = 10


# ###############################
# 1. Modelselection
# ###############################
nirs.criterion.slopedist = criterionSlopeDist(nirs.data, COLUMN_ID_NM_FROM, COLUMN_ID_NM_TO)
selectedCriterion = nirs.criterion.slopedist

# select wavelengths with highes variablilty
plot(WAVELENGTHS[1:length(WAVELENGTHS)-1], selectedCriterion, type = "l", col=1)
abline(h=THRESHOLD, col=2)
fullModel.features = selectFeatures(nirs.data, WAVELENGTHS, THRESHOLD, selectedCriterion)

# build model from selected wavelength
fullModel.formula = buildFormula(fullModel.features, "nm", RESPONSEVAR)

# modelselection based on Mallow's CP from full model
nirs.subsets = regsubsets(fullModel.formula, nirs.data, really.big=T, nvmax=nrow(fullModel.features)+1, method="backward")
nirs.lm.full = lm(fullModel.formula, data=nirs.data)
nirs.cp = min(summary(nirs.subsets)$cp)
nirs.optModelId = which.min(summary(nirs.subsets)$cp)
nirs.lm.opt = getModelFromSubsets(nirs.subsets, nirs.optModelId, "N", nirs.data)
nirs.spse.true = calculateTrueSpse(nirs.lm.opt, nirs.optModelId, NUM_ROWS)
nirs.spse.true

# plot slope with selected features
plot(WAVELENGTHS[1:length(WAVELENGTHS)-1], selectedCriterion, type = "l", col=1)
selfeat = as.integer(substr(names(nirs.lm.opt$coefficients)[-1],3,6))
sapply(selfeat, addSelectedFeatureToPlot, 3)

# ###############################
# 2. Simulation
# ###############################
R2 = calculateR2(nirs.lm.opt, nirs.data, NUM_ROWS, 1000)
R2

sim.models = data.frame(matrix(nrow=length(SAMPLE_SIZES), ncol=SIM_ITERATIONS))
sim.cp = data.frame(matrix(nrow=length(SAMPLE_SIZES), ncol=SIM_ITERATIONS))
for (run in 1:SIM_ITERATIONS) {
  for (i in 1:length(SAMPLE_SIZES)) {
    simResult = simulateOnDataSubset(nirs.lm.opt, nirs.data, fullModel.formula, fullModel.features, SAMPLE_SIZES[i], 1000)
    sim.models[i, run] = simResult["id"]
    sim.cp[i, run] = simResult["cp"]
  }
  print(paste(run, "/", SIM_ITERATIONS, sep=""))
}

# ############################
# calculate SPSEs and compare
# ############################
sigma2.tilde.full = calcauletSigma2TildeFull(nirs.lm.full, nrow(fullModel.features)+1, NUM_ROWS)

spse_true = calculateTrueSpse(nirs.lm.opt, nirs.optModelId, nirs.lm.full, NUM_ROWS, nrow(fullModel.features)+1)

sim.spse.est = (sim.cp + SAMPLE_SIZES) * sigma2.tilde.full
plot(x=150, y=mean(as.vector(t(sim.spse.est[1,]))), xlim=c(0,550), ylim=c(0,2), pch=16, col=1)
for (i in 2:length(SAMPLE_SIZES)) {
  points(x=SAMPLE_SIZES[i], y=mean(as.vector(t(sim.spse.est[i,]))), pch=16, col=1)
  print(mean(as.vector(t(sim.spse.est[i,]))))
}
points(x=NUM_ROWS, y=spse_true, col=2, pch=16)

# ###########################
# plots
# ##########################
# 1. correlation plot
simN = simulate(nirs.lm.opt, nirs.data, 533, 1000)
plot(nirs.data$N, simN, xlim=c(0,0.8), ylim=c(0,0.8))
abline(a=0,b=1, col=2)
