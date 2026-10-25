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
## numbers of the MMRM section, written by vignettes/precompute-mmrm.R
mm <- readRDS("mmrm-sensitivity-numbers.rds")
fmt <- function(x, digits = 2) formatC(x, format = "f", digits = digits)
## The slowest chunks are not evaluated (eval=FALSE, echo=FALSE); their
## output is precomputed by vignettes/precompute-chunks.R (to keep the
## check time down) and showChunk() writes code and output as Sweave
## would. It stops if a chunk's code no longer matches the code the
## output was computed from.
chunkOutput <- if (file.exists("robustlmm-3.5.0-chunks.rds"))
    readRDS("robustlmm-3.5.0-chunks.rds")
showChunk <- function(label) {
    src <- readLines("robustlmm-3.5.0.Rnw")
    s <- grep(paste0("^<<", label, "[,>]"), src)
    e <- s + match(TRUE, grepl("^@", src[-seq_len(s)]))
    x <- chunkOutput[[label]]
    if (!identical(src[(s + 1):(e - 1)], x$code))
        stop("the stored output of chunk '", label, "' is out of date; ",
             "re-run vignettes/precompute-chunks.R")
    out <- "\\begin{Schunk}"
    input <- FALSE
    for (st in x$steps) {
        if (!input) out <- c(out, "\\begin{Sinput}")
        out <- c(out, paste0(c(getOption("prompt"),
                               rep(getOption("continue"),
                                   length(st$src) - 1L)), st$src))
        input <- TRUE
        if (length(st$output)) {
            out <- c(out, "\\end{Sinput}", "\\begin{Soutput}", st$output,
                     "\\end{Soutput}")
            input <- FALSE
        }
    }
    if (input) out <- c(out, "\\end{Sinput}")
    writeLines(c(out, "\\end{Schunk}"))
}
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
confint(fm, method = "Wald", vcov_type = "sandwich")


###################################################
### code chunk number 5: confint-df (eval = FALSE)
###################################################
## confint(fm, method = "Wald", df = "satterthwaite")


###################################################
### code chunk number 6: confint-df-output
###################################################
showChunk("confint-df")


###################################################
### code chunk number 7: confint-boot (eval = FALSE)
###################################################
## ## parametric BCa bootstrap (needs the 'confintROB' package)
## confint(fm, method = "BCa",  nsim = 1000, seed = 20260601)
## ## wild bootstrap
## confint(fm, method = "boot", boot.type = "wild", nsim = 1000,
##         seed = 20260601)


###################################################
### code chunk number 8: anova-single
###################################################
anova(fm)


###################################################
### code chunk number 9: anova-two
###################################################
fm0 <- rlmer(Reaction ~ 1 + (Days | Subject), sleepstudy,
             method = "DASvar")
anova(fm0, fm)


###################################################
### code chunk number 10: anova-ddf
###################################################
anova(fm, ddf = "satterthwaite")


###################################################
### code chunk number 11: anova-boot (eval = FALSE)
###################################################
## anova(fm0, fm, test = "boot", nsim = 1000, seed = 20260601)


###################################################
### code chunk number 12: summary-df
###################################################
summary(fm, df = "satterthwaite")$coefficients


###################################################
### code chunk number 13: emmeans-df
###################################################
emmeans::emtrends(fm, ~ 1, var = "Days")


###################################################
### code chunk number 14: predict-interval
###################################################
nd <- data.frame(Days = 0:2, Subject = "308")
suppressWarnings(predict(fm, newdata = nd, interval = "confidence"))
suppressWarnings(predict(fm, newdata = nd, interval = "prediction"))


###################################################
### code chunk number 15: cooks-obs (eval = FALSE)
###################################################
## fpen <- rlmer(diameter ~ (1 | plate) + (1 | sample), Penicillin,
##               method = "DASvar")
## cd <- cooks.distance(fpen)
## sort(cd, decreasing = TRUE)[1:5]


###################################################
### code chunk number 16: cooks-obs-output
###################################################
showChunk("cooks-obs")


###################################################
### code chunk number 17: cooks-cluster
###################################################
cdc <- cooks.distance(fm, groups = "Subject")
sort(cdc, decreasing = TRUE)[1:5]


###################################################
### code chunk number 18: hatvalues
###################################################
range(hatvalues(fm))


###################################################
### code chunk number 19: cluster-fig
###################################################
op <- par(mar = c(4, 4, 1, 1))
dotchart(sort(cdc), xlab = "Cluster-level Cook's distance",
         ylab = "Subject", pch = 19, cex = 0.7)
par(op)


###################################################
### code chunk number 20: ransac-init (eval = FALSE)
###################################################
## fr <- rlmer(Reaction ~ Days + (Days | Subject), sleepstudy,
##             method = "DASvar", init = "ransac")
## fixef(fr)


###################################################
### code chunk number 21: ransac-init-output
###################################################
showChunk("ransac-init")


###################################################
### code chunk number 22: ransac-lme4
###################################################
rs <- ransac_lme4(Reaction ~ Days + (Days | Subject), sleepstudy,
                  K = 20L, seed = 20260601)
rs$K_used


###################################################
### code chunk number 23: ransac-consensus (eval = FALSE)
###################################################
## fc <- rlmer_ransac(Reaction ~ Days + (Days | Subject), sleepstudy,
##                    K = 20L, n_starts = 3L, seed = 20260601)
## attr(fc, "consensus")


###################################################
### code chunk number 24: ransac-consensus-output
###################################################
showChunk("ransac-consensus")


###################################################
### code chunk number 25: design-weights
###################################################
dat <- sleepstudy
dat$Days[1] <- 40                 # an extreme, high-leverage design point
fw <- rlmer(Reaction ~ Days + (Days | Subject), dat,
            method = "DASvar", design.weights = "mcd")
sum(getME(fw, "design.weights") < 1)


###################################################
### code chunk number 26: structured-data
###################################################
set.seed(20260623)
datasets <- generateRepeatedMeasuresDatasets(
    1, numberOfSubjects = 40, numberOfVisits = 4, numberOfReplicates = 3,
    structure = "ar1", marginalSd = 1.5, correlation = 0.6,
    trueBeta = 10, trueSigma = 0.7)
d <- datasets$generateData(1)


###################################################
### code chunk number 27: structured-fits (eval = FALSE)
###################################################
## rcs  <- rlmer(y ~ 1 + cs(0 + visit | subject),  d, method = "DASvar")
## rar1 <- rlmer(y ~ 1 + ar1(0 + visit | subject), d, method = "DASvar")
## VarCorr(rcs)
## VarCorr(rar1)


###################################################
### code chunk number 28: structured-fits-output
###################################################
showChunk("structured-fits")


###################################################
### code chunk number 29: mmrm-fit (eval = FALSE)
###################################################
## library("mmrm")
## fitMmrm <- mmrm(FEV1 ~ RACE + SEX + ARMCD * AVISIT + us(AVISIT | USUBJID),
##                 data = fev_data)


###################################################
### code chunk number 30: mmrm-rlmer (eval = FALSE)
###################################################
## fitRob <- rlmer(FEV1 ~ RACE + SEX + ARMCD * AVISIT + (1 | USUBJID),
##                 data = fev_data)
## fev_data$visit <- as.numeric(fev_data$AVISIT) - 2.5
## fitRobSlope <- rlmer(FEV1 ~ RACE + SEX + ARMCD * AVISIT +
##                          (visit | USUBJID),
##                      data = fev_data, method = "DASvar")


###################################################
### code chunk number 31: mmrm-emmeans (eval = FALSE)
###################################################
## library("emmeans")
## contrast(emmeans(fitMmrm, ~ ARMCD | AVISIT), "revpairwise")
## contrast(emmeans(fitRob, ~ ARMCD | AVISIT), "revpairwise")


###################################################
### code chunk number 32: mmrm-contrast (eval = FALSE)
###################################################
## L <- c("ARMCDTRT" = 1, "ARMCDTRT:AVISITVIS4" = 1)
## contrast_test(fitRob, L)
## contrast_test(fitRob, L, vcov_type = "sandwich")


###################################################
### code chunk number 33: mmrm-table
###################################################
ctab <- mm$contrasts
lab <- c(mmrm = "MMRM, unstructured covariance",
         lmer = "\\code{lmer}, random intercept",
         rlmer = "\\code{rlmer}, random intercept",
         rlmerSandwich = "\\quad with the cluster sandwich",
         slope = "\\code{rlmer}, random intercept and slope",
         slopeSandwich = "\\quad with the cluster sandwich")
cat("\\begin{table}[htb]\n\\centering\n",
    "\\begin{tabular}{lrrrr}\n\\toprule \n",
    "Fit & Estimate & SE & df & $p$ \\\\\n\\midrule \n", sep = "")
for (r in rownames(ctab))
    cat(lab[[r]], " & ", fmt(ctab[r, "estimate"]), " & ",
        fmt(ctab[r, "se"]), " & ",
        if (is.finite(ctab[r, "df"])) fmt(ctab[r, "df"], 0) else "--",
        " & ", if (ctab[r, "p"] < 0.001) "$<$0.001" else fmt(ctab[r, "p"], 3),
        " \\\\\n", sep = "")
cat("\\bottomrule \n\\end{tabular}\n",
    "\\caption[Treatment difference at the last visit]{Treatment ",
    "difference in FEV1 at the last visit. The standard errors without ",
    "the sandwich are model-based; ``--'' is the normal reference.}\n",
    "\\label{tab:mmrm}\n\\end{table}\n", sep = "")


###################################################
### code chunk number 34: mmrm-augment (eval = FALSE)
###################################################
## aug <- generics::augment(fitRob)
## aug[order(aug$.wt_e)[1:6], ]
## table(aug$AVISIT[aug$.wt_e < 0.8])


###################################################
### code chunk number 35: mmrm-cooks (eval = FALSE)
###################################################
## cd <- cooks.distance(fitRobSlope, groups = TRUE)
## head(sort(cd, decreasing = TRUE))


###################################################
### code chunk number 36: sessionInfo
###################################################
## (session info intentionally not printed to keep the vignette stable)
