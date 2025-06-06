#######################################################
## Basic accessor functions                          ##
#######################################################

## Get Z
## (slot Zt in reModule might be changed to Z, so
## don't write it explicitly)
##
## @title Get Z from reModule
## @param object merMod object
## @param t transpose or not
.Zt <- function(object) object@pp$Zt
getZ <- function(object, t = FALSE) {
    if (t) .Zt(object) else t(.Zt(object))
}

## Get X
##
## @title Get X from predModule
## @param object merMod object
## @param t transpose or not
.X <- function(object) object@pp$X
getX <- function(object, t = FALSE) {
    if (t) t(.X(object)) else .X(object)
}

## Get rho-function used for residuals
##
## @title Get rho_e
## @param object merMod object
## @param which add \dQuote{sigma} for rho.sigma.e
## @export
rho.e <- function(object, which = "default") {
    switch(which,
           sigma=object@rho.sigma.e,
           default=,object@rho.e)
}

## Get rho-function for used random effects
##
## @title Get rho_b
## @param object merMod object
## @param which add \dQuote{sigma} for rho.sigma.e
## @export
rho.b <- function(object, which = "default") {
    ret <- switch(which,
                  sigma=object@rho.sigma.b,
                  default=,object@rho.b)
    ## backwards compatibility:
    if (inherits(ret, "psi_func"))
        ret <- rep.int(list(ret), length(object@blocks))
    ## add names
    names(ret) <- names(object@cnms)
    ret
}

## Get Lambda
##
## @title Get Lambda
## @param object merMod object
Lambda <- function(object) {
    ## FIXME: which theta?
    if (class(object)[1] == "lmerMod") t(object@pp$Lambdat) else t(object@pp$Lambdat())
}

## Get U_b
##
## @title Get U_b
## @param object merMod object
U_b <- function(object) {
   if (class(object)[1] == "lmerMod") t(object@pp$Lambdat) else object@pp$U_b
}

## Get Lind
##
## @title Get Lind
## @param object merMod object
Lind <- function(object) {
  object@pp$Lind
}

## Get lower
##
## @title Get lower
## @param object merMod object
lower <- function(object) {
    object@lower
}

## Get indices of the random effects corresponding to zero variance components
##
## @title Get indices of r.e. corresp. to zero v.c.
## @param object merMod object
getZeroU <- function(object) .zeroB(object)

## Get various numbers of parameters / lengths of vectors.
##
## \itemize{
## \item b, u: length of random effects vector
## \item beta, coef: number of fixed effects
## \item theta: length of vector theta
## \item r, e: number of observations
## }
##
## @title length
## @param x merMod object
## @param what length is requested
len <- function(x, what) switch(what,
                                u=,
                                b=nrow(.Zt(x)),
                                coef=,
                                beta=length(x@beta),
                                theta=length(x@theta),
                                r=,
                                e=length(x@resp$y),
                                stop("unknown length"))

.nobsLmerMod <- function(object, ...) len(object, "e")
##' @importFrom stats nobs
##' @export
nobs.rlmerMod <- .nobsLmerMod

### Get REML (so that we are not coercing all the time)
##' @importFrom lme4 isREML
.isREML <- function(x, ...) {
    as.logical(x@devcomp$dims["REML"])
}

##' The per-observation residuals are returned, i.e., the difference of the
##' observation and the fitted value including random effects. With type one can
##' specify whether the weights should be used or not.
##'
##' @title Get residuals
##' @param object rlmerMod object
##' @param type type of residuals
##' @param scaled scale residuals by residual standard deviation (=scale
##'   parameter)?
##' @param ... ignored
##' @examples
##' \dontrun{
##'   fm <- rlmer(Yield ~ (1|Batch), Dyestuff)
##'   stopifnot(all.equal(resid(fm, type="weighted"),
##'                       resid(fm) * getME(fm, "w_e")))
##' }
##' @importFrom stats residuals resid
##' @export
residuals.rlmerMod <- function(object, type = c("response", "weighted"),
                               scaled=FALSE, ...) {
    type <- match.arg(type)
    r <- switch(type,
                ## FIXME: really??
                response = object@resp$wtres,
                weighted = wgt.e(object) * object@resp$wtres,
                stop("unknown type of residual"))
    if (is.null(nm <- rownames(model.frame(object)))) nm <- seq_along(r)
    names(r) <- nm
    if (scaled) r <- r/sigma(object)
    if (!is.null(na.action <- attr(model.frame(object),"na.action")))
        r <- naresid(na.action,r)
    r
}

### Get sigma (so that we are not coercing all the time)
.sigma <- function(object, ...)
  object@pp$sigma
##' @importFrom lme4 sigma
##' @export
sigma.rlmerMod <- .sigma


### Get deviance
.deviance <- function(object, ...)
    stop("Deviance is not defined for rlmerMod objects")
##' @export
deviance.rlmerMod <- .deviance

.mu <- function(object)
{
    ## Purpose: calculate mu of respModule
    ## ----------------------------------------------------------------------
    ## Arguments: object: lmerMod object
    ## ----------------------------------------------------------------------
    ## Author: Manuel Koller, Date: 11 Apr 2011, 11:40

    ## FIXME: ?? offset will be added in updateMu
    crossprod(.Zt(object), .b(object))@x + (.X(object) %*% .fixef(object))[,1]
}

### Get fixed effects
.fixef <- function(object)
  object@pp$beta

### Get u
b.s <- .u <- function(object, ...)
  object@pp$b.s
u.rlmerMod <- function(object, ...) {
    ret <- b.s(object)
    names(ret) <- dimnames(getZ(object))[[2]]
    ret
}

u.lmerMod <- function(object, ...) object@u

### Get b
.b <- function(object, ...)
  object@pp$b.r
b.rlmerMod <- function(object, ...) {
    ret <- .b(object)
    names(ret) <- dimnames(getZ(object))[[2]]
    ret
}
b.lmerMod <- function(object, ...) {
    ret <- crossprod(getME(object, "Lambdat"), getME(object, "u"))
    names(ret) <- dimnames(getME(object, "Zt"))[[1]]
    ret
}

### Get ranef
##' @importFrom nlme ranef
##' @export
ranef.rlmerMod <- function(object, ...) {
    ## FIXME: add postVar, drop and whichel arguments
    b <- b(object)
    ret <- uArrangedNames(object, b.s = b)
    class(ret) <- "ranef.rlmerMod"
    ret
}

## return u as list arranged like ranef group ranefs with the same subject
uArrangedNames <- function(object, b.s = b.s(object)) {
    ret <- list()
    for (id in unique(names(object@cnms))) {
        lid <- id == names(object@cnms)
        lr <- lapply(object@idx[lid], function(bidx) {
            lret <- b.s[bidx]
            dim(lret) <- dim(bidx)
            ## add rownames
            colnames(lret) <- names(b.s)[bidx[1,]]
            as.data.frame(t(lret)) })
        lr <- do.call(cbind, lr)
        colnames(lr) <- unlist(object@cnms[lid])
        ret <- c(ret, list(lr))
        names(ret)[length(ret)] <- id
    }
    ret
}
## same as uArrangedNames, but do not set names and do not group by id name
uArranged <- function(object, b.s = b.s(object), idx = object@idx) {
    ret <- lapply(idx, function(bidx) {
        lret <- b.s[bidx]
        dim(lret) <- dim(bidx)
        t(lret)
    })
    ret
}

## Construct names of individual theta/sd:cor components
##
## @param object a fixed model
## @param diag.only include only diagonal elements?
## @param old (logical) give backward-compatible results?
## @param prefix a character vector with two elements giving the prefix for
##   diagonal (e.g. \dQuote{sd}) and off-diagonal (e.g. \dQuote{cor}) elements
##   ## @export
tnames <- function(object,diag.only=FALSE,old=TRUE,prefix=NULL) {
    if (old) {
        nc <- c(unlist(mapply(function(g,e) {
            mm <- outer(e,e,paste,sep=".")
            diag(mm) <- e
            mm <- if (diag.only) diag(mm) else mm[lower.tri(mm,diag=TRUE)]
            paste(g,mm,sep=".")
        },
        names(object@cnms),object@cnms)))
        return(nc)
    } else {
        pfun <- function(g,e) {
            mm <- outer(e,e,paste,sep=".")
            mm[] <- paste(mm,g,sep="|")
            if (!is.null(prefix)) mm[] <- paste(prefix[2],mm,sep="_")
            diag(mm) <- paste(e,g,sep="|")
            if (!is.null(prefix))  diag(mm) <- paste(prefix[1],diag(mm),sep="_")
            mm <- if (diag.only) diag(mm) else mm[lower.tri(mm,diag=TRUE)]
        }
        nc <- c(unlist(mapply(pfun,names(object@cnms),object@cnms)))
        return(nc)
    }
}

##' Extract or Get Generalize Components from a Fitted Mixed Effects Model
##'
##' Extract (or \dQuote{get}) \dQuote{components} -- in a generalized
##' sense -- from a fitted mixed-effects model, i.e. from an object
##' of class \code{\linkS4class{rlmerMod}} or \code{\link[lme4:merMod-class]{merMod}}.
##'
##' The function \code{theta} is short for \code{getME(, "theta")}.
##'
##' The goal is to provide \dQuote{everything a user may want} from a fitted
##' \code{rlmerMod} object \emph{as far} as it is not available by methods, such
##' as \code{\link[lme4:fixef.merMod]{fixef}},
##' \code{\link[lme4:ranef.merMod]{ranef}}, \code{\link[lme4:vcov.merMod]{vcov}}, etc.
##'
##' @param object a fitted mixed-effects model of class
##' \code{\linkS4class{rlmerMod}}, i.e. typically the result of
##' \code{\link{rlmer}()}.
##' @param name a character string specifying the name of the
##' \dQuote{component}.  Possible values are:\cr
##' \describe{
##'     \item{\code{"X"}:}{fixed-effects model matrix}
##'     \item{\code{"Z"}:}{random-effects model matrix}
##'     \item{\code{"Zt"}:}{transpose of random-effects model matrix}
##'     \item{\code{"Ztlist"}:}{list of components of the transpose of the random-effects model matrix,
##'              separated by individual variance component}
##'     \item{\code{"mmList"}:}{list of raw model matrices associated with random effects terms}
##'     \item{\code{"y"}:}{response vector}
##'     \item{\code{"mu"}:}{conditional mean of the response}
##'     \item{\code{"u"}:}{conditional mode of the \dQuote{spherical} random effects variable}
##'     \item{\code{"b.s"}:}{synonym for \dQuote{u}}
##'     \item{\code{"b"}:}{conditional mode of the random effects variable}
##'     \item{\code{"Gp"}:}{groups pointer vector.  A pointer to the beginning of each group
##'               of random effects corresponding to the random-effects terms.}
##'     \item{\code{"Tp"}:}{theta pointer vector.  A pointer to the beginning
##'               of the theta sub-vectors corresponding to the
##'               random-effects terms, beginning with 0 and including
##'               a final element giving the total number of random effects}
##'     \item{\code{"Lambda"}:}{relative covariance factor of the random effects.}
##'     \item{\code{"U_b"}:}{synonym for \dQuote{Lambda}}
##'     \item{\code{"Lambdat"}:}{transpose of the relative covariance factor of the random effects.}
##'     \item{\code{"Lind"}:}{index vector for inserting elements of \eqn{\theta}{theta} into the
##'                 nonzeros of \eqn{\Lambda}{Lambda}}
##'     \item{\code{"A"}:}{Scaled sparse model matrix (class
##'      \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}}) for
##'      the unit, orthogonal random effects, \eqn{U},
##'       equal to \code{getME(.,"Zt") \%*\% getME(.,"Lambdat")}}
##'     \item{\code{"sigma"}:}{residual standard error}
##'     \item{\code{"flist"}:}{a list of the grouping variables (factors) involved in the random effect terms}
##'     \item{\code{"fixef"}:}{fixed-effects parameter estimates}
##'     \item{\code{"beta"}:}{fixed-effects parameter estimates (identical to the result of \code{\link[lme4:fixef.merMod]{fixef}}, but without names)}
##'     \item{\code{"theta"}:}{random-effects parameter estimates: these are parameterized as the relative Cholesky factors of each random effect term}
##'     \item{\code{"ST"}:}{A list of S and T factors in the TSST' Cholesky
##'               factorization of the relative variance matrices of the random
##'               effects associated with each random-effects term.  The unit lower
##'               triangular matrix, \eqn{T}, and the diagonal matrix, \eqn{S}, for
##'               each term are stored as a single matrix with diagonal elements
##'               from \eqn{S} and off-diagonal elements from \eqn{T}.}
##'     \item{\code{"is_REML"}:}{returns \code{TRUE} for rlmerMod-objects (for compatibility with lme4)}
##'     \item{\code{"n_rtrms"}:}{number of random-effects terms}
##'     \item{\code{"n_rfacs"}:}{number of distinct random-effects grouping factors}
##'     \item{\code{"N"}:}{number of rows of \code{X}}
##'     \item{\code{"n"}:}{length of the response vector, \code{y}}
##'     \item{\code{"p"}:}{number of columns of the fixed effects model matrix, \code{X}}
##'     \item{\code{"q"}:}{number of columns of the random effects model matrix, \code{Z}}
##'     \item{\code{"p_i"}:}{numbers of columns of the raw model matrices, \code{mmList}}
##'     \item{\code{"l_i"}:}{numbers of levels of the grouping factors}
##'     \item{\code{"q_i"}:}{numbers of columns of the term-wise model matrices, \code{ZtList}}
##'     \item{\code{"k"}:}{number of random effects terms}
##'     \item{\code{"m_i"}:}{numbers of covariance parameters in each term}
##'     \item{\code{"m"}:}{total number of covariance parameters, i.e., the
##'                        same as \code{dim@nth} below.}
##'     \item{\code{"cnms"}:}{the \dQuote{component names}, a \sQuote{list}.}
##'     \item{\code{"devcomp"}:}{a list consisting of a named numeric vector,
##'         \code{cmp}, and a named integer vector, \code{dims}, describing
##'         the fitted model.  The elements of \code{cmp} are:\cr
##'         \describe{
##'             \item{ldL2}{always NA, for consistency with lme4 output}
##'             \item{ldRX2}{always NA, for consistency with lme4 output}
##'             \item{wrss}{always NA, for consistency with lme4 output}
##'             \item{ussq}{always NA, for consistency with lme4 output}
##'             \item{pwrss}{always NA, for consistency with lme4 output}
##'             \item{drsum}{always NA, for consistency with lme4 output}
##'             \item{REML}{always NA, for consistency with lme4 output}
##'             \item{dev}{always NA, for consistency with lme4 output}
##'             \item{sigmaML}{always NA, for consistency with lme4 output}
##'             \item{sigmaREML}{REML estimate of residual standard deviation}
##'         } The elements of \code{dims} are:\cr
##'         \describe{
##'             \item{N}{number of rows of \code{X}}
##'             \item{n}{length of \code{y}}
##'             \item{p}{number of columns of \code{X}}
##'             \item{nmp}{\code{n-p}}
##'             \item{nth}{length of \code{theta}}
##'             \item{q}{number of columns of \code{Z}}
##'             \item{nAGQ}{see \code{\link[lme4]{glmer}}}
##'             \item{compDev}{see \code{\link[lme4]{glmerControl}}}
##'             \item{useSc}{\code{TRUE} if model has a scale parameter}
##'             \item{reTrms}{number of random effects terms}
##'             \item{REML}{\code{0} indicates the model was fitted by maximum
##'                 likelihood, any other positive integer indicates fitting by
##'                 restricted maximum likelihood}
##'             \item{GLMM}{\code{TRUE} if a GLMM}
##'             \item{NLMM}{\code{TRUE} if an NLMM}
##'         }
##'     }
##'     \item{\code{"offset"}:}{model offset}
##'     \item{\code{"lower"}:}{lower bounds on random-effects model
##'         parameters (i.e, "theta" parameters). In order to constrain
##'         random effects covariance matrices to be semi-positive-definite,
##'         this vector is equal to 0 for elements of
##'         the \code{theta} vector corresponding to diagonal elements of
##'         the Cholesky factor, \code{-Inf}
##'         otherwise. (\code{getME(.,"lower")==0} can be used as a test to
##'         identify diagonal elements, as in \code{isSingular}.)
##'     }
##'     \item{\code{"rho_e"}:}{rho function used for the residuals}
##'     \item{\code{"rho_b"}:}{list of rho functions used for the random effects}
##'     \item{\code{"rho_sigma_e"}:}{rho function used for the residuals when estimating sigma}
##'     \item{\code{"rho_sigma_b"}:}{list of rho functions used for the random effects when estimating the covariance parameters}
##'     \item{\code{"M"}:}{list of matrices, blocks of the Henderson's equations and the matrices used for computing the linear approximations of the estimates of beta and spherical random effects.}
##'     \item{\code{"w_e"}:}{robustness weights associated with the observations}
##'     \item{\code{"w_b"}:}{robustness weights associated with the spherical random effects, returned in the same format as \code{\link[lme4:ranef.merMod]{ranef}()}}
##'     \item{\code{"w_b_vector"}:}{robustness weights associated with the spherical random effects, returned as one long vector}
##'     \item{\code{"w_sigma_e"}:}{robustness weights associated with the observations when estimating sigma}
##'     \item{\code{"w_sigma_b"}:}{robustness weights associated with the spherical random effects when estimating the covariance parameters, returned in the same format as \code{\link[lme4:ranef.merMod]{ranef}()}}
##'     \item{\code{"w_sigma_b_vector"}:}{robustness weights associated with the spherical random effects when estimating the covariance parameters, returned as one long vector}
##'
##'      %% -- keep at the very end:
##'      \item{\code{"ALL"}:}{get all of the above as a \code{\link[base]{list}}.}
##' }
##' @param ... potentially further arguments; not here.
##' @return Unspecified, as very much depending on the \code{\link[base]{name}}.
##' @seealso \code{\link[stats]{getCall}()};
##' more standard methods for rlmerMod objects, such as \code{\link[lme4:ranef.merMod]{ranef}},
##' \code{\link[lme4:fixef.merMod]{fixef}}, \code{\link[lme4:vcov.merMod]{vcov}}, etc.:
##' see \code{methods(class="rlmerMod")}
##' @keywords utilities
##' @examples
##' ## shows many methods you should consider *before* using getME():
##' methods(class = "rlmerMod")
##'
##' ## doFit = FALSE to speed up example
##' (fm1 <- rlmer(Reaction ~ Days + (Days|Subject), sleepstudy,
##'               method="DASvar", doFit=FALSE))
##' Z <- getME(fm1, "Z")
##' stopifnot(is(Z, "CsparseMatrix"),
##'           c(180,36) == dim(Z),
##' 	  all.equal(fixef(fm1), b1 <- getME(fm1, "beta"),
##' 		    check.attributes=FALSE, tolerance = 0))
##'
##' ## A way to get *all* getME()s :
##' ## internal consistency check ensuring that all work:
##' parts <- getME(fm1, "ALL")
##' str(parts, max=2)
##' stopifnot(identical(Z,  parts $ Z),
##'           identical(b1, parts $ beta))
## % S3 generic now imported:
##' @importFrom lme4 getME
##' @rdname getME
##' @export
##' @method getME rlmerMod
getME.rlmerMod <-
    function(object,
             name = c("X", "Z", "Zt", "Ztlist", "mmList", "y", "mu",
                      "u", "b.s", "b",
                      "Gp", "Tp", "Lambda", "Lambdat", "Tlist",
                      "A", "U_b", "Lind", "sigma", "flist",
                      "fixef", "beta", "theta", "ST", "is_REML",
                      "n_rtrms", "n_rfacs", "N", "n", "p", "q",
                      "p_i", "l_i", "q_i", "k", "m_i", "m",
                      "cnms", "devcomp", "offset", "lower",
                      "rho_e", "rho_b", "rho_sigma_e", "rho_sigma_b",
                      "M", "w_e", "w_b", "w_b_vector", "w_sigma_e",
                      "w_sigma_b", "w_sigma_b_vector"), ...)
{
    if(missing(name)) stop("'name' must not be missing")
    ## Deal with multiple names -- "FIXME" is inefficiently redoing things
    if (length(name <- as.character(name)) > 1) {
        names(name) <- name
        return(lapply(name, getME, object = object))
    }
    if(name == "ALL") ## recursively get all provided components
        return(sapply(eval(formals()$name),
                      getME.rlmerMod, object=object, simplify=FALSE))
    stopifnot(is(object, "rlmerMod"))
    name <- match.arg(name)
    rsp  <- object@resp
    PR   <- object@pp
    dc   <- object@devcomp
    cmp  <- dc $ cmp
    cnms <- object@cnms
    dims <- dc $ dims
    Tpfun <- function(cnms) {
	ltsize <- function(n) n*(n+1)/2 # lower triangle size
	cLen <- cumsum(ltsize(vapply(cnms,length, 1L)))
	setNames(c(0, cLen),
		 c("beg__", names(cnms))) ## such that diff( Tp ) is well-named
    }
    if (any(name == c("p_i", "q_i", "m_i")))
        p_i <- vapply(mmList.rlmerMod(object), ncol, 1L)
    if (any(name == c("l_i", "q_i")))
        l_i <- vapply(object@flist, nlevels, 1L)
    switch(name,
	   "X" = getX(object, t=FALSE),
	   "Z" = getZ(object, t=FALSE),
	   "Zt"= getZ(object, t=TRUE),
           "Ztlist" =
       {
           getInds <- function(i) {
               n <- diff(object@Gp)[i]      ## number of elements in this block
               nt <- length(cnms[[i]]) ## number of REs
               inds <- lapply(seq(nt),seq,to=n,by=nt)  ## pull out individual RE indices
               inds <- lapply(inds,function(x) x + object@Gp[i])  ## add group offset
           }
           inds <- do.call(c,lapply(seq_along(cnms),getInds))
           setNames(lapply(inds,function(i) PR$Zt[i,]),
                    tnames(object,diag.only=TRUE))
       },
	   "mmList" = mmList.rlmerMod(object),
           "y" = rsp$y,
           "mu"= rsp$mu,
           "u" =,
           "b.s" = b.s(object),
           "b" = b(object),
	       "L" = PR$L(),
	   "Lambda"= ,
           "U_b" = PR$U_b,
	   "Lambdat"= PR$Lambdat(),
           "A"= PR$Lambdat() %*% PR$Zt,
           "Lind" = PR$Lind,
           "sigma" = sigma(object),
           "Gp" = object@Gp,
           "Tp" = Tpfun(cnms) ,
           "flist" = object@flist,
	       "fixef" = fixef(object),
	   "beta" = object@beta,
           "theta"= theta(object),
	   "ST" = setNames(lme4::vec2STlist(object@theta, n = lengths(cnms)),
	                 names(cnms)), Tlist = {
	                     nc <- lengths(cnms)
	                     nt <- length(nc)
	                     ans <- vector("list", nt)
	                     names(ans) <- names(nc)
	                     pos <- 0L
	                     th <- object@theta
	                     for (i in 1:nt) {
	                         nci <- nc[i]
	                         tt <- matrix(0, nci, nci)
	                         inds <- lower.tri(tt, diag = TRUE)
	                         nthi <- sum(inds)
	                         tt[inds] <- th[pos + seq_len(nthi)]
	                         pos <- pos + nthi
	                         ans[[i]] <- tt
	                     }
	                     ans
	                 },
	   "n_rtrms" = length(object@flist),
           "n_rfacs" = length(object@flist),
	   "N" = dims[["N"]], "n" = dims[["n"]], "p" = dims[["p"]],
	   "q" = dims[["q"]], "p_i" = p_i, "l_i" = l_i,
	   "q_i" = p_i * l_i, "k" = length(cnms),
	   "m_i" = choose(p_i + 1, 2), "m" = dims[["nth"]],
           "cnms" = cnms,
           "devcomp" = dc,
           "offset" = rsp$offset,
           "lower" = object@lower,
           "rho_e" = rho.e(object),
           "rho_b" = rho.b(object),
           "rho_sigma_e" = rho.e(object, "sigma"),
           "rho_sigma_b" = rho.b(object, "sigma"),
           "M" = PR$ M(),
           "w_e" = wgt.e(object),
           "w_b" = uArrangedNames(object, wgt.b(object)),
           "w_b_vector" = wgt.b(object),
           "w_sigma_e" = wgt.e(object, use.rho.sigma=TRUE),
           "w_sigma_b" = uArrangedNames(object, wgt.b(object, use.rho.sigma=TRUE)),
           "w_sigma_b_vector" = wgt.b(object, use.rho.sigma=TRUE),
           "is_REML" = TRUE,
	   "..foo.." =# placeholder!
	   stop(gettextf("'%s' is not implemented yet",
			 sprintf("getME(*, \"%s\")", name))),
	   ## otherwise
	   stop(sprintf("Mixed-Effects extraction of '%s' is not available for class \"%s\"",
			name, class(object))))
}## {getME}

##' The function \code{theta} is short for \code{getME(, "theta")}.
##'
##' @rdname getME
##' @examples
##' stopifnot(all.equal(theta(fm1), getME(fm1, "theta")))
##' @export
theta <- function(object) {
    if (is(object, "rlmerMod")) {
        ## add names like lme4
        tt <- object@pp$theta
        nc <- c(unlist(mapply(function(g,e) {
            mm <- outer(e,e,paste,sep=".")
            diag(mm) <- e
            mm <- mm[lower.tri(mm,diag=TRUE)]
            paste(g,mm,sep=".")
        }, names(object@cnms),object@cnms)))
        names(tt) <- nc
        tt
    } else getME(object, "theta")
}

.zeroB <- function(object, pp = object@pp)
  pp$zeroB

.U_eZU_b <- function(object)
  object@pp$U_eZU_b

.U_eX <- function(object)
  object@pp$U_eX

.U_e <- function(object)
  object@pp$U_e

.U_b <- function(object)
  object@pp$U_b

.Lambda_b <- function(object)
  object@pp$Lambda_b

.kappa_b <- function(object)
  object@pp$kappa_b

.kappa_e <- function(object)
  object@pp$kappa_e
