## plot functions for rlmerMod objects

globalVariables("theoretical", add=TRUE)

## internal functions, called via plot.rlmerMod
## TA-plots inclusive coloring for weights
ta <- function(obj, title="", add.line="none", ...) {
    data <- data.frame(fitted = fitted(obj),
                       resid = resid(obj),
                       weights = if (is(obj, "rlmerMod")) getME(obj, "w_e") else 1)
    plt <- ggplot2::ggplot(data, ggplot2::aes(fitted, resid))
    if (add.line == "below")
        plt <- plt + ggplot2::geom_hline(yintercept = 0, ...)
    if (title != "") plt <- plt + ggplot2::ggtitle(title)
    plt <-
        if (is(obj, "rlmerMod"))
            plt + ggplot2::geom_point(ggplot2::aes(color = weights))
        else plt + ggplot2::geom_point()
    if (add.line == "above")
        plt <- plt + ggplot2::geom_hline(yintercept = 0, ...)
    plt
}
## QQ-plots inclusive coloring for weights
qq <- function(obj, type = c("resid", "ranef"), title="",
               multiply.weights=FALSE, add.line="none", ...) {
    type <- match.arg(type)
    val0 <- switch(type,
                  resid = list(resid = data.frame(resid=resid(obj))),
                  ranef = ranef(obj))
    ord0 <- numeric(0)
    for (level in names(val0)) {
        val1 <- val0[[level]]
        for (col in colnames(val1)) {
            val <- val1[[col]]
            ord <- order(val)
            name <- if (ncol(val1) == 1) level else paste(level, col, sep="/")
            data <- data.frame(level = name,
                               sample = val[ord],
                               theoretical = qnorm(ppoints(length(val))))
            data0 <- if (exists("data0")) rbind(data0, data) else data
            ord0 <- c(ord0, length(ord0) + ord)
        }
    }
    data0$weights <-
        if (is(obj, "rlmerMod")) {
            switch(type, resid = getME(obj, "w_e"),
                   ranef = getME(obj, "w_b_vector"))[ord0]
        } else 1
    if (multiply.weights) data0$sample <- data0$sample * data0$weights
    plt <- ggplot2::ggplot(data0, ggplot2::aes(theoretical, sample))
    if (add.line == "below")
        plt <- plt + ggplot2::geom_qq_line(ggplot2::aes(sample = sample), ...)
    if (title != "") plt <- plt + ggplot2::ggtitle(title)
    plt <-
        if (is(obj, "rlmerMod"))
            plt + ggplot2::geom_point(ggplot2::aes(color = weights))
        else plt + ggplot2::geom_point()
    if (add.line == "above")
        plt <- plt + ggplot2::geom_qq_line(ggplot2::aes(sample = sample), ...)
    if (length(levels(data0$level)) > 1) plt + ggplot2::facet_wrap(~ level) else plt
}
## scatterplots for correlated random effects
rsc <- function(obj, title="") {
    r <- ranef(obj)
    if (is(obj, "rlmerMod")) w <- getME(obj, "w_b")
    qn <- function(n) paste("`", n, "`", sep="")
    plots <- list()
    for (g in names(r)) {
        df <- r[[g]]
        nc <- ncol(df)
        if (nc < 2) next
        if (is(obj, "rlmerMod")) df <- cbind(df, weights=w[[g]][,1])
        ## create a plot for all combinations
        for (i in 1L:(nc-1L)) {
            for (j in (i+1L):nc) {
                lplt <- ggplot2::ggplot(df, ggplot2::aes_string(x=qn(colnames(df)[i]),
                                                      y=qn(colnames(df)[j])))
                if (title != "") lplt <- lplt + ggplot2::ggtitle(sprintf(title, g))
                lplt <- if (is(obj, "rlmerMod"))
                    lplt + ggplot2::geom_point(ggplot2::aes(color = weights))
                else lplt + ggplot2::geom_point()
                plots <- c(plots, list(lplt))
            }
        }
    }
    plots
}

##' Diagnostic plots for objects of class \code{rlmerMod} and \code{lmerMod}.
##'
##' The robustness weights for estimating the fixed and random effects are used
##' in the plots, e.g., the ones returned by \code{getME(object, "w_e")} and
##' \code{getME(object, "w_b")}.
##'
##' @title Plot Method for "rlmerMod" objects.
##' @param x an object as created by \code{rlmer} or \code{rlmer}; or an object
##'   as created by \code{plot.rlmerMod}
##' @param y currently ignored.
##' @param which integer number between 1 and 4 to specify which plot is
##'   desired.
##' @param title Titles for the different plots. The fourth item can be a format
##'   string passed to \code{sprintf} to add the name of the current group.
##' @param multiply.weights multiply the residuals / random effects with the
##'   robustness weights when producing the Q-Q plots.
##' @param ask waits for user input before displaying each plot.
##' @param add.line add reference line to plots, use \code{"above"} or
##'   \code{"below"} to show the line above or below the points. Hide the line
##'   with \code{"none"}.
##' @param ... passed on to \code{\link[ggplot2]{geom_hline}} and
##'   \code{\link[ggplot2]{geom_qq_line}}, to customize how the line is drawn.
##' @return a list of plots of class \code{\link[ggplot2]{ggplot}} that can be
##'   used for further modification before plotting (using \code{print}).
##' @seealso \code{\link{getME}}, \code{\link[ggplot2]{ggplot}}
##' @examples
##' \dontrun{
##'   rfm <- rlmer(Yield ~ (1|Batch), Dyestuff)
##'   plot(rfm)
##'   fm <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
##'   plot.rlmerMod(fm)
##' }
##' @export plot.rlmerMod
##' @method plot rlmerMod
##' @export
plot.rlmerMod <- function(x, y=NULL, which=1:4,
                          title = c("Fitted Values vs. Residuals",
                              "Normal Q-Q vs. Residuals",
                              "Normal Q-Q vs. Random Effects",
                              "Scatterplot of Random Effects for Group \"%s\""),
                          multiply.weights=FALSE,
                          add.line=c("above", "below", "none"),
                          ...) {
    if (!inherits(x, "rlmerMod") & !inherits(x, "lmerMod"))
        stop("Use only with 'rlmerMod' and 'lmerMod' objects")
    if (!isPackageInstalled("ggplot2"))
        stop("Package 'ggplot2' is required to use 'plot.rlmerMod'.",
             "Please install it using 'install.packages(\"ggplot2\")'.")
    show <- rep.int(FALSE, 4)
    if (!is.numeric(which) || any(which < 1) || any(which > 4))
        stop("'which' must be in 1:4")
    show[which] <- TRUE
    add.line <- match.arg(add.line)
    plots <- list()

    if (show[1])
        plots[[1]] <- ta(x, title=title[1], add.line=add.line, ...)
    if (show[2])
        plots <- c(plots, list(qq(x, type="resid", title=title[2],
                                  multiply.weights=multiply.weights,
                                  add.line=add.line, ...)))
    if (show[3])
        plots <- c(plots, list(qq(x, type="ranef", title=title[3],
                                  multiply.weights=multiply.weights,
                                  add.line=add.line, ...)))
    if (show[4])
        plots <- c(plots, rsc(x, title=title[4]))

    class(plots) <- "rlmerMod_plots"
    plots
}

##' @rdname plot.rlmerMod
##' @method print rlmerMod_plots
##' @export
print.rlmerMod_plots <- function(x, ask=interactive() & length(x) > 1, ...) {
    if (ask) {
        oldAsk <- devAskNewPage(ask=TRUE)
        on.exit(devAskNewPage(ask=oldAsk))
    }
    invisible(lapply(x, print))
}

### --- FIXME: These plot functions are just little hacks, should redo
### --- them in ggplot at some point.

##' @importFrom lattice dotplot
##' @export
dotplot.ranef.rlmerMod <- getS3method("dotplot", "ranef.mer")

##' @importFrom graphics plot
##' @export
plot.ranef.rlmerMod <- getS3method("plot", "ranef.mer")

##' @importFrom lattice qqmath
##' @export
qqmath.ranef.rlmerMod <- getS3method("qqmath", "ranef.mer")

##' @importFrom graphics plot
##' @export
plot.coef.rlmerMod <- function(x, y, ...) {
    ## remove non-varying columns from frames
    reduced <- lapply(x, function(el)
		      el[, !sapply(el, function(cc) all(cc == cc[1]))])
    plot.ranef.rlmerMod(reduced, ...)
}

##' @importFrom lattice dotplot
##' @export
dotplot.coef.rlmerMod <- function(x, data, ...) {
    mc <- match.call()
    mc[[1]] <- as.name("dotplot.ranef.rlmerMod")
    eval(mc)
}
