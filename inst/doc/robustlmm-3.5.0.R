### R code from vignette source 'robustlmm-3.5.0.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(prompt = "R> ", continue = "+  ", width = 76,
        useFancyQuotes = FALSE, digits = 4)
suppressPackageStartupMessages(library("robustlmm"))
set.seed(20260623)
## Numbers quoted from the simulation studies are precomputed into this
## RDS by vignettes/precompute-3.5.0.R and pulled in with \Sexpr{} below,
## so they always match the studies and are never hand-typed.
nums <- readRDS("robustlmm-3.5.0-numbers.rds")
fmt <- function(x, digits = 2) formatC(x, format = "f", digits = digits)
## Sweave does not echo R warnings into the vignette the way an
## interactive session does; this helper prints them next to the call
## that triggers them.
echoWarnings <- function(expr)
    withCallingHandlers(expr, warning = function(w) {
        writeLines(strwrap(paste("Warning:", conditionMessage(w)),
                           width = 74))
        invokeRestart("muffleWarning")
    })


###################################################
### code chunk number 2: base-fit
###################################################
fm <- rlmer(Reaction ~ Days + (Days | Subject), sleepstudy,
            method = "DASvar")


###################################################
### code chunk number 3: vcov-sandwich
###################################################
vcov(fm, type = "default")
echoWarnings(vcov(fm, type = "sandwich"))


###################################################
### code chunk number 4: confint-wald
###################################################
confint(fm, method = "Wald")
suppressWarnings(  # silence the small-J sandwich note shown above
    confint(fm, method = "Wald", vcov_type = "sandwich"))


###################################################
### code chunk number 5: confint-df
###################################################
confint(fm, method = "Wald", df = "satterthwaite")


###################################################
### code chunk number 6: confint-boot (eval = FALSE)
###################################################
## ## parametric BCa bootstrap (needs the 'confintROB' package)
## confint(fm, method = "BCa",  nsim = 1000, seed = 20260601)
## ## wild bootstrap
## confint(fm, method = "boot", boot.type = "wild", nsim = 1000,
##         seed = 20260601)


###################################################
### code chunk number 7: anova-single
###################################################
anova(fm)


###################################################
### code chunk number 8: anova-two
###################################################
fm0 <- rlmer(Reaction ~ 1 + (Days | Subject), sleepstudy,
             method = "DASvar")
anova(fm0, fm)


###################################################
### code chunk number 9: anova-ddf
###################################################
anova(fm, ddf = "satterthwaite")


###################################################
### code chunk number 10: anova-boot (eval = FALSE)
###################################################
## anova(fm0, fm, test = "boot", nsim = 1000, seed = 20260601)


###################################################
### code chunk number 11: summary-df
###################################################
summary(fm, df = "satterthwaite")$coefficients


###################################################
### code chunk number 12: emmeans-df
###################################################
emmeans::emtrends(fm, ~ 1, var = "Days")


###################################################
### code chunk number 13: predict-interval
###################################################
nd <- data.frame(Days = 0:2, Subject = "308")
suppressWarnings(predict(fm, newdata = nd, interval = "confidence"))
suppressWarnings(predict(fm, newdata = nd, interval = "prediction"))


###################################################
### code chunk number 14: cooks-obs
###################################################
fpen <- rlmer(diameter ~ (1 | plate) + (1 | sample), Penicillin,
              method = "DASvar")
cd <- cooks.distance(fpen)
sort(cd, decreasing = TRUE)[1:5]


###################################################
### code chunk number 15: cooks-cluster
###################################################
cdc <- cooks.distance(fm, groups = "Subject")
sort(cdc, decreasing = TRUE)[1:5]


###################################################
### code chunk number 16: hatvalues
###################################################
range(hatvalues(fm))


###################################################
### code chunk number 17: cluster-fig
###################################################
op <- par(mar = c(4, 4, 1, 1))
dotchart(sort(cdc), xlab = "Cluster-level Cook's distance",
         ylab = "Subject", pch = 19, cex = 0.7)
par(op)


###################################################
### code chunk number 18: ransac-init
###################################################
fr <- rlmer(Reaction ~ Days + (Days | Subject), sleepstudy,
            method = "DASvar", init = "ransac")
fixef(fr)


###################################################
### code chunk number 19: ransac-lme4
###################################################
rs <- ransac_lme4(Reaction ~ Days + (Days | Subject), sleepstudy,
                  K = 20L, seed = 20260601)
rs$K_used


###################################################
### code chunk number 20: ransac-consensus
###################################################
fc <- rlmer_ransac(Reaction ~ Days + (Days | Subject), sleepstudy,
                   K = 20L, n_starts = 3L, seed = 20260601)
attr(fc, "consensus")


###################################################
### code chunk number 21: design-weights
###################################################
dat <- sleepstudy
dat$Days[1] <- 40                 # an extreme, high-leverage design point
fw <- rlmer(Reaction ~ Days + (Days | Subject), dat,
            method = "DASvar", design.weights = "mcd")
sum(getME(fw, "design.weights") < 1)


###################################################
### code chunk number 22: structured-data
###################################################
set.seed(20260623)
datasets <- generateRepeatedMeasuresDatasets(
    1, numberOfSubjects = 40, numberOfVisits = 4, numberOfReplicates = 3,
    structure = "ar1", marginalSd = 1.5, correlation = 0.6,
    trueBeta = 10, trueSigma = 0.7)
d <- datasets$generateData(1)


###################################################
### code chunk number 23: structured-fits
###################################################
rcs  <- rlmer(y ~ 1 + cs(0 + visit | subject),  d, method = "DASvar")
rar1 <- rlmer(y ~ 1 + ar1(0 + visit | subject), d, method = "DASvar")
VarCorr(rcs)
VarCorr(rar1)


###################################################
### code chunk number 24: sessionInfo
###################################################
## (session info intentionally not printed to keep the vignette stable)
