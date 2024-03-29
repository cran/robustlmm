% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot.rlmerMod}
\alias{plot.rlmerMod}
\alias{print.rlmerMod_plots}
\title{Plot Method for "rlmerMod" objects.}
\usage{
\method{plot}{rlmerMod}(
  x,
  y = NULL,
  which = 1:4,
  title = c("Fitted Values vs. Residuals", "Normal Q-Q vs. Residuals",
    "Normal Q-Q vs. Random Effects", "Scatterplot of Random Effects for Group \\"\%s\\""),
  multiply.weights = FALSE,
  add.line = c("above", "below", "none"),
  ...
)

\method{print}{rlmerMod_plots}(x, ask = interactive() & length(x) > 1, ...)
}
\arguments{
\item{x}{an object as created by \code{rlmer} or \code{rlmer}; or an object
as created by \code{plot.rlmerMod}}

\item{y}{currently ignored.}

\item{which}{integer number between 1 and 4 to specify which plot is
desired.}

\item{title}{Titles for the different plots. The fourth item can be a format
string passed to \code{sprintf} to add the name of the current group.}

\item{multiply.weights}{multiply the residuals / random effects with the
robustness weights when producing the Q-Q plots.}

\item{add.line}{add reference line to plots, use \code{"above"} or
\code{"below"} to show the line above or below the points. Hide the line
with \code{"none"}.}

\item{...}{passed on to \code{\link[ggplot2]{geom_hline}} and
\code{\link[ggplot2]{geom_qq_line}}, to customize how the line is drawn.}

\item{ask}{waits for user input before displaying each plot.}
}
\value{
a list of plots of class \code{\link[ggplot2]{ggplot}} that can be
  used for further modification before plotting (using \code{print}).
}
\description{
Diagnostic plots for objects of class \code{rlmerMod} and \code{lmerMod}.
}
\details{
The robustness weights for estimating the fixed and random effects are used
in the plots, e.g., the ones returned by \code{getME(object, "w_e")} and
\code{getME(object, "w_b")}.
}
\examples{
\dontrun{
  rfm <- rlmer(Yield ~ (1|Batch), Dyestuff)
  plot(rfm)
  fm <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
  plot.rlmerMod(fm)
}
}
\seealso{
\code{\link{getME}}, \code{\link[ggplot2]{ggplot}}
}
