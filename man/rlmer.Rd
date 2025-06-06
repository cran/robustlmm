% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rlmer.R, R/lmer.R
\name{rlmer}
\alias{rlmer}
\alias{lmerNoFit}
\title{Robust Scoring Equations Estimator for Linear Mixed Models}
\usage{
rlmer(
  formula,
  data,
  ...,
  method = c("DAStau", "DASvar"),
  setting,
  rho.e,
  rho.b,
  rho.sigma.e,
  rho.sigma.b,
  rel.tol = 1e-08,
  max.iter = 40 * (r + 1)^2,
  verbose = 0,
  doFit = TRUE,
  init
)

lmerNoFit(formula, data = NULL, ..., initTheta)
}
\arguments{
\item{formula}{a two-sided linear formula object describing the
fixed-effects part of the model, with the response on the left of a
\code{~} operator and the terms, separated by \code{+} operators, on the
right.  The vertical bar character \code{"|"} separates an expression for
a model matrix and a grouping factor.}

\item{data}{an optional data frame containing the variables named in
\code{formula}.  By default the variables are taken from the environment
from which \code{lmer} is called.}

\item{...}{Additional parameters passed to lmer to find the initial
estimates. See \code{\link[lme4]{lmer}}.}

\item{method}{method to be used for estimation of theta and sigma, see
Details.}

\item{setting}{a string specifying suggested choices for the arguments
\code{rho.e}, \code{rho.sigma.e}, \code{rho.b} and \code{rho.sigma.b}.
Use \code{"RSEn"} (the default) or \code{"RSEa"}. Both use
\code{\link{smoothPsi}} for all the \dQuote{rho} arguments. For
\code{rho.sigma.e}, squared robustness weights are used (see
\code{\link{psi2propII}}). \code{"RSEn"} uses the same tuning parameter as
for \code{rho.e}, which leads to higher robustness but lower efficiency.
\code{"RSEa"} adjusts the tuning parameter for higher asymptotic efficiency
which results in lower robustness (\code{k = 2.28} for default \code{rho.e}).
For diagonal random effects covariance matrices, \code{rho.sigma.b} is
treated exactly as \code{rho.sigma.e}. For block diagonal random effects
covariance matrices (with correlation terms), regular robustness weights
are used for \code{rho.sigma.b}, not squared ones, as they're not needed.
But the tuning parameters are adjusted for both \code{rho.b} and
\code{rho.sigma.b} according to the dimensions of the blocks (for both
\code{"RSEn"} or \code{"RSEa"}). For a block of dimension 2 (e.g.,
correlated random intercept and slope) \code{k = 5.14} is used.}

\item{rho.e}{object of class psi_func, specifying the functions to use for
the huberization of the residuals.}

\item{rho.b}{object of class psi_func or list of such objects (see Details),
specifying the functions to use for the huberization of the random
effects.}

\item{rho.sigma.e}{object of class psi_func, specifying the weight functions
to use for the huberization of the residuals when estimating the variance
components, use the \code{\link{psi2propII}} function to specify squared
weights and custom tuning parameters.}

\item{rho.sigma.b}{(optional) object of class psi_func or list of such
objects, specifying the weight functions to use for the huberization of
the random effects when estimating the variance components (see Details).
Use \code{\link{psi2propII}} to specify squared weights and custom tuning
parameters or \code{\link{chgDefaults}} for regular weights for variance
components including correlation parameters.}

\item{rel.tol}{relative tolerance used as criteria in the fitting process.}

\item{max.iter}{maximum number of iterations allowed.}

\item{verbose}{verbosity of output. Ranges from 0 (none) to 3 (a lot of
output)}

\item{doFit}{logical scalar. When \code{doFit = FALSE} the model is not fit
but instead a structure with the model matrices for the random-effects
terms is returned (used to speed up tests). When \code{doFit = TRUE}, the
default, the model is fit immediately.}

\item{init}{optional lmerMod- or rlmerMod-object to use for starting values,
a list with elements \sQuote{fixef}, \sQuote{u}, \sQuote{sigma},
\sQuote{theta}, or a function producing an lmerMod object.}

\item{initTheta}{parameter to initialize theta with (optional)}
}
\value{
object of class rlmerMod.
}
\description{
Robust estimation of linear mixed effects models, for hierarchical nested
and non-nested, e.g., crossed, datasets.
}
\details{
\describe{ \item{Overview:}{

This function implements the Robust Scoring Equations estimator for linear
mixed effect models. It can be used much like the function
\code{\link[lme4]{lmer}} in the package \code{lme4}. The supported models
are the same as for \code{\link[lme4]{lmer}} (gaussian family only). The
robust approach used is based on the robustification of the scoring
equations and an application of the Design Adaptive Scale approach.

Example analyses and theoretical details on the method are available in the
vignette (see \code{vignette("rlmer")}).

Models are specified using the \code{formula} argument, using the same
syntax as for \code{\link[lme4]{lmer}}. Additionally, one also needs to
specify what robust scoring or weight functions are to be used (arguments
starting with \code{rho.}). By default a smoothed version of the Huber
function is used. Furthermore, the \code{method} argument can be used to
speed up computations at the expense of accuracy of the results. }

\item{Computation methods:}{

Currently, there are two different methods available for fitting models.
They only differ in how the consistency factors for the Design Adaptive
Scale estimates are computed. Available fitting methods for theta and
sigma.e: \itemize{

\item \code{DAStau} (default): For this method, the consistency factors are
computed using numerical quadrature. This is slower but yields more accurate
results. This is the direct analogue to the DAS-estimate in robust linear
regression.

\item \code{DASvar}: This method computes the consistency factors using a
direct approximation which is faster but less accurate. For complex models
with correlated random effects with more than one correlation term, this is
the only method available.

} }

\item{Weight functions:}{

The tuning parameters of the weight functions \dQuote{rho} can be used to
adjust robustness and efficiency of the resulting estimates (arguments
\code{rho.e}, \code{rho.b}, \code{rho.sigma.e} and \code{rho.sigma.b}).
Better robustness will lead to a decrease of the efficiency. With the default
setting, \code{setting = "RSEn"}, the tuning parameters are set to yield
estimates with approximately 95\% efficiency for the fixed effects. The
variance components are estimated with a lower efficiency but better
robustness properties.

One has to use different weight functions and tuning parameters for simple
variance components and for such including correlation parameters. By
default, they are chosen appropriately to the model at hand. However, when
using the \code{rho.sigma.e} and \code{rho.sigma.b} arguments, it is up to
the user to specify the appropriate function. See
\code{\link{asymptoticEfficiency}} for methods to find tuning parameters
that yield a given asymptotic efficiency. \itemize{

\item For simple variance components and the residual error scale use the
function \code{\link{psi2propII}} to change the tuning parameters. This is
similar to Proposal 2 in the location-scale problem (i.e., using the
squared robustness weights of the location estimate for the scale estimate;
otherwise the scale estimate is not robust).

\item For multi-dimensional blocks of random effects modeled, e.g.,
a model with correlated random intercept and slope, (referred to as
block diagonal case below), use the \code{\link{chgDefaults}} function to
change the tuning parameters. The parameter estimation problem is
multivariate, unlike the case without correlation where the problem was
univariate. For the employed estimator, this amounts to switching from
simple scale estimates to estimating correlation matrices. Therefore
different weight functions have to be used. Squaring of the weights (using
the function \code{\link{psi2propII}}) is no longer necessary. To yield
estimates with the same efficiency, the tuning parameters for the
block diagonal are larger than for the simple case. Tables of tuning parameters
are given in Table 2 and 3 of the vignette (\code{vignette("rlmer")}).

} }

\item{Recommended tuning parameters:}{

For a more robust estimate, use \code{setting = "RSEn"} (the default). For
higher efficiency, use \code{setting = "RSEa"}. The settings described in
the following paragraph are used when \code{setting = "RSEa"} is specified.

For the smoothed Huber function the tuning parameters to get approximately
95\% efficiency are \eqn{k=1.345}{k=1.345} for \code{rho.e} and
\eqn{k=2.28}{k=2.28} for \code{rho.sigma.e} (using the squared version). For
simple variance components, the same can be used for \code{rho.b} and
\code{rho.sigma.b}. For variance components including correlation
parameters, use \eqn{k=5.14}{k=5.14} for both \code{rho.b} and
\code{rho.sigma.b}. Tables of tuning parameter are given in Table 2 and 3 of
the vignette (\code{vignette("rlmer")}). }

\item{Specifying (multiple) weight functions:}{

If custom weight functions are specified using the argument \code{rho.b}
(\code{rho.e}) but the argument \code{rho.sigma.b} (\code{rho.sigma.e}) is
missing, then the squared weights are used for simple variance components
and the regular weights are used for variance components including
correlation parameters. The same tuning parameters will be used when
\code{setting = "RSEn"} is used. To get
higher efficiency either use \code{setting = "RSEa"} (and only set arguments
\code{rho.e} and \code{rho.b}). Or specify the tuning parameters by hand
using the \code{\link{psi2propII}} and \code{\link{chgDefaults}} functions.

To specify separate weight functions \code{rho.b} and \code{rho.sigma.b} for
different variance components, it is possible to pass a list instead of a
psi_func object. The list entries correspond to the groups as shown by
\code{VarCorr(.)} when applied to the model fitted with \code{lmer}. A set
of correlated random effects count as just one group. }

\item{\code{lmerNoFit}:}{

The \code{lmerNoFit} function can be used to get trivial starting values.
This is mainly used to verify the algorithms to reproduce the fit by
\code{\link[lme4]{lmer}} when starting from trivial initial values. } }
}
\examples{
## dropping of VC
system.time(print(rlmer(Yield ~ (1|Batch), Dyestuff2, method="DASvar")))

\dontrun{
  ## Default method "DAStau"
  system.time(rfm.DAStau <- rlmer(Yield ~ (1|Batch), Dyestuff))
  summary(rfm.DAStau)
  ## DASvar method (faster, less accurate)
  system.time(rfm.DASvar <- rlmer(Yield ~ (1|Batch), Dyestuff,
                                  method="DASvar"))
  ## compare the two
  compare(rfm.DAStau, rfm.DASvar)

  ## Fit variance components with higher efficiency
  ## psi2propII yields squared weights to get robust estimates
  ## this is the same as using rlmer's argument `setting = "RSEa"`
  rlmer(diameter ~ 1 + (1|plate) + (1|sample), Penicillin,
        rho.sigma.e = psi2propII(smoothPsi, k = 2.28),
        rho.sigma.b = psi2propII(smoothPsi, k = 2.28))

  ## use chgDefaults for variance components including
  ## correlation terms (regular, non squared weights suffice)
  ## this is the same as using rlmer's argument `setting = "RSEa"`
  rlmer(Reaction ~ Days + (Days|Subject), sleepstudy,
        rho.sigma.e = psi2propII(smoothPsi, k = 2.28),
        rho.b = chgDefaults(smoothPsi, k = 5.14, s=10),
        rho.sigma.b = chgDefaults(smoothPsi, k = 5.14, s=10))
}

\dontrun{
  ## start from lmer's initial estimate, not its fit
  rlmer(Yield ~ (1|Batch), Dyestuff, init = lmerNoFit)
}
}
\seealso{
\code{\link[lme4]{lmer}}, \code{vignette("rlmer")}
}
\author{
Manuel Koller, with thanks to Vanda Lourenço for improvements.
}
\keyword{models}
