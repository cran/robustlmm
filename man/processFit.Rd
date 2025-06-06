% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processFit.R
\name{processFit}
\alias{processFit}
\alias{processFit.lmerMod}
\alias{processFit.rlmerMod}
\alias{processFit.heavyLme}
\alias{processFit.lqmm}
\alias{processFit.rlme}
\alias{processFit.varComprob}
\title{Process Fitted Objects}
\usage{
processFit(
  obj,
  all = FALSE,
  coefs = TRUE,
  stdErrors = all,
  tValues = all,
  sigma = TRUE,
  thetas = TRUE,
  b = all,
  meanB = all,
  meanAbsB = all,
  residuals = all,
  converged = TRUE,
  numWarnings = all,
  procTime = all,
  ...
)

\method{processFit}{lmerMod}(
  obj,
  all = FALSE,
  coefs = TRUE,
  stdErrors = all,
  tValues = all,
  sigma = TRUE,
  thetas = TRUE,
  b = all,
  meanB = all,
  meanAbsB = all,
  residuals = all,
  converged = TRUE,
  numWarnings = all,
  procTime = all,
  ...
)

\method{processFit}{rlmerMod}(
  obj,
  all = FALSE,
  coefs = TRUE,
  stdErrors = all,
  tValues = all,
  sigma = TRUE,
  thetas = TRUE,
  b = all,
  meanB = all,
  meanAbsB = all,
  residuals = all,
  converged = TRUE,
  numWarnings = all,
  procTime = all,
  ...
)

\method{processFit}{heavyLme}(
  obj,
  all = FALSE,
  coefs = TRUE,
  stdErrors = all,
  tValues = all,
  sigma = TRUE,
  thetas = TRUE,
  b = all,
  meanB = all,
  meanAbsB = all,
  residuals = all,
  converged = TRUE,
  numWarnings = all,
  procTime = all,
  ...
)

\method{processFit}{lqmm}(
  obj,
  all = FALSE,
  coefs = TRUE,
  stdErrors = all,
  tValues = all,
  sigma = TRUE,
  thetas = TRUE,
  b = all,
  meanB = all,
  meanAbsB = all,
  residuals = all,
  converged = TRUE,
  numWarnings = all,
  procTime = all,
  ...
)

\method{processFit}{rlme}(
  obj,
  all = FALSE,
  coefs = TRUE,
  stdErrors = all,
  tValues = all,
  sigma = TRUE,
  thetas = TRUE,
  b = all,
  meanB = all,
  meanAbsB = all,
  residuals = all,
  converged = TRUE,
  numWarnings = all,
  procTime = all,
  ...
)

\method{processFit}{varComprob}(
  obj,
  all = FALSE,
  coefs = TRUE,
  stdErrors = all,
  tValues = all,
  sigma = TRUE,
  thetas = TRUE,
  b = all,
  meanB = all,
  meanAbsB = all,
  residuals = all,
  converged = TRUE,
  numWarnings = all,
  procTime = all,
  isInterceptCorrelationSlopeModel,
  ...
)
}
\arguments{
\item{obj}{object returned by the fitting method.}

\item{all}{logical, shorthand to enable all exports.}

\item{coefs}{logical, if true coefficients are added to export.}

\item{stdErrors}{logical, if true, standard errors are added to export.}

\item{tValues}{logical, if true, t-values are added to export.}

\item{sigma}{logical, if true, sigma is added to export.}

\item{thetas}{logical, if true, thetas are added to export.}

\item{b}{scalar logical or index vector, if true, all random effects are
added to export. If an index vector is given, then only the corresponding
random effects are added to the export. The same order as in \code{lmer}
is used for all methods.}

\item{meanB}{logical, if true, the mean of the random effects is added to
the export.}

\item{meanAbsB}{logical, if true, the mean of the absolute value of the
random effects is added to the export.}

\item{residuals}{scalar logical or index vector, similar to argument
\code{b}, just returning the residuals.}

\item{converged}{logical, if true, convergence code is added to export.}

\item{numWarnings}{logical, if true, the number of warnings generated during
the fitting process is added to export.}

\item{procTime}{logical, if true, time needed to fit object is added to
export.}

\item{...}{optional parameters used for some implementations.}

\item{isInterceptCorrelationSlopeModel}{optional logical, can be used to
override the assumption that a model with three variance components can be
interpreted as having intercept, correlation and slope.}
}
\value{
List with extracted values, most items can be suppressed
  to save disk space.
\item{\code{label}: }{Name of fitting method used to create the fit}
\item{\code{datasetIndex}: }{Index of the dataset in the dataset list}
\item{\code{coefficients}: }{Vector of estimated fixed-effects coefficients of the fitted model}
\item{\code{standardErrors}: }{Vector of estimated standard errors of the fixed-effects coefficients}
\item{\code{tValues}: }{Vector of t-Values (or z-Values depending on fitting method)
  of the fixed-effects coefficients}
\item{\code{sigma}: }{Estimated residual standard error}
\item{\code{thetas}: }{Vector of random-effects parameter estimates. As parameterized as by
  \code{\link[lme4]{lmer}} and \code{\link{rlmer}}.}
\item{\code{b}: }{Vector of requested predicted random-effects.}
\item{\code{meanB}: }{Vector of means of the predicted random-effects.}
\item{\code{meanAbsB}: }{Vector of means of the absolute values of the  predicted random-effects.}
\item{\code{residuals}: }{Vector of requested residuals.}
\item{\code{converged}: }{Convergence status as reported by the fitting method. \code{0} means converged.
  If not available, \code{NA} is used. Other values are to be interpreted carefully as codes
  vary from method to method.}
\item{\code{numberOfWarnings}: }{the number of warnings generated during the fitting process.}
\item{\code{proc.time}: }{Vector of times (user, system, elapsed) as reported by \code{\link[base]{proc.time}}
  required to fit the model.}
}
\description{
Methods to process fitted objects and convert into a data structure that is
useful in post-processing.
}
\details{
Warning. \code{processFit.varComprob} uses simplistic logic to
  convert from the parameterisation used in the robustvarComp package to
  \code{theta} as used in \code{\link[lme4]{lmer}} and \code{\link{rlmer}}. If
  there are three variance components, the code assumes that they are
  intercept, correlation and slope. Otherwise the code assumes that the
  variance components are independent. Exports \code{b} and \code{residuals}
  are not supported.
}
\examples{
  set.seed(1)
  oneWay <- generateAnovaDatasets(1, 1, 10, 4,
                                  lmeFormula = y ~ 1,
                                  heavyLmeRandom = ~ 1,
                                  heavyLmeGroups = ~ Var2,
                                  lqmmRandom = ~ 1,
                                  lqmmGroup = "Var2",
                                  groups = cbind(rep(1:4, each = 10), rep(1:10, 4)),
                                  varcov = matrix(1, 4, 4),
                                  lower = 0)
  processFit(fitDatasets_lmer(oneWay)[[1]], all = TRUE)
  processFit(fitDatasets_rlmer_DASvar(oneWay)[[1]], all = TRUE)
  \dontrun{
    processFit(fitDatasets_heavyLme(oneWay)[[1]], all = TRUE)
  }
  if (require(lqmm)) {
    processFit(fitDatasets_lqmm(oneWay)[[1]], all = TRUE)
  }
  \dontrun{
    processFit(fitDatasets_varComprob_compositeTau(oneWay)[[1]], all = TRUE)
  }
}
