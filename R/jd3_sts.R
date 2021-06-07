#' @include jd3_ts.R jd3_rslts.R protobuf.R
#' @checkmate
NULL

#' Title
#'
#' @param y 
#' @param X 
#' @param level 
#' @param slope 
#' @param cycle 
#' @param noise 
#' @param seasonal 
#' @param diffuse.regs 
#' @param tol 
#'
#' @return
#' @export
#'
#' @examples
std<-rjd3modelling::td(12, c(1992,1), 228, c(1,1,1,1,2,3,0))
ts<-function(y, X=NULL, level=1, slope=1, cycle=-1, noise=1
              , seasonal=c("Trigonometric", "Dummy", "Crude", "HarrisonStevens", "Fixed", "Unused"), diffuse.regs=T, tol=1e-9){
  
  if (!is.ts(y)){
    stop("y must be a time series")
  }
  seasonal<-match.arg(seasonal)
  jts<-ts_r2jd(y)
  jx<-matrix_r2jd(X)
  jsts<-.jcall("demetra/sts/r/Bsm", "Ldemetra/sts/BasicStructuralModel;", "process", jts, jx,
              as.integer(level), as.integer(slope), as.integer(cycle), as.integer(noise), seasonal, as.logical(diffuse.regs), tol)
  buffer<-.jcall("demetra/sts/r/Bsm", "[B", "toBuffer", jsts)
  p<-RProtoBuf::read(sts.Bsm, buffer)
  return (p2r_sts_rslts(p))
}

#' Title
#'
#' @param y Series
#' @param model Model for calendar effects
#' \itemize{
#'   \item{td2: }{leap year + week days (week-end derived)}
#'   \item{td3: }{leap year + week days + saturdays (sundays derived)}
#'   \item{td7: }{leap year + all days (sundays derived)} 
#'   \item{full: }{td3 + easter effect}
#'   \item{none: }{no calendar effect}
#'   }
#' @param nf 
#'
#' @return
#' @export
#'
#' @examples
sts.forecast<-function(y, model=c("none", "td2", "td3", "td7", "full"), nf=12){
  model<-match.arg(model)
  if (!is.ts(y)){
    stop("y must be a time series")
  }
  jf<-.jcall("demetra/sts/r/Bsm", "Ldemetra/math/matrices/MatrixType;", "forecast", ts_r2jd(y), model, as.integer((nf)))
  return (matrix_jd2r(jf))
  
}

p2r_sts_rslts<-function(p){
  
  return (structure(list(
    description=p2r_sts_description(p$description),
    estimation=p2r_sts_estimation(p$estimation),
    decomposition=p2r_sts_components(p$components)),
    class="JD3STS")
  )
}

p2r_sts_estimation<-function(p){
  return (list(
    y=p$y,
    X=p2r_matrix(p$x),
    parameters=p2r_parameters_estimation(p$parameters),
    b=p$b,
    bvar=p2r_matrix(p$bcovariance),
    likelihood=p2r_diffuselikelihood(p$likelihood),
    res=p$residuals))
}

p2r_sts_description<-function(p){
  return (list(
    log=p$log,
    preadjustment = enum_extract(modelling.LengthOfPeriod, p$preadjustment),
    bsm=p2r_spec_bsm(p$bsm),
    variables=p2r_variables(p$variables)))
}

p2r_sts_components<-function(p){
  return (list(
    level=p2r_sts_component(p$level),
    slope=p2r_sts_component(p$slope),
    cycle=p2r_sts_component(p$cycle),
    seasonal=p2r_sts_component(p$seasonal),
    noise=p2r_sts_component(p$noise)
  ))
}

p2r_sts_component<-function(p){
  if (is.null(p)) return (NULL) else return (p$as.list())
}




p2r_spec_bsm<-function(p){
  return (list(
    level=p2r_parameter(p$level),
    slope=p2r_parameter(p$slope),
    seas=p2r_parameter(p$seas),
    seasmodel=enum_extract(sts.SeasonalModel, p$seasonal_model),
    noise=p2r_parameter(p$noise),
    cycle=p2r_parameter(p$cycle),
    cyclelength=p2r_parameter(p$cycle_period),
    cyclefactor=p2r_parameter(p$cycle_factor)
  ))
  
}

#' Title
#'
#' @param m 
#'
#' @return
#' @export
#'
#' @examples
print.JD3STS<-function(m){
  cat("Structural time series", "\n\n")
  cat("Variances:\n")
  s<-m$description$bsm$level
  if (! is.null(s)) cat("level: ", format(round(s$value, 6), scientific = FALSE), "\n")
  s<-m$description$bsm$slope
  if (! is.null(s)) cat("slope: ", format(round(s$value, 6), scientific = FALSE), "\n")
  s<-m$description$bsm$seas
  if (! is.null(s)) cat("seasonal: ", format(round(s$value, 6), scientific = FALSE), "\n")
  s<-m$description$bsm$noise
  if (! is.null(s)) cat("noise: ", format(round(s$value, 6), scientific = FALSE), "\n\n")
  s<-m$description$bsm$cycle
  if (! is.null(s)) {
    cat("cycle: ", format(round(s$value, 6), scientific = FALSE), "\n\n")
  }
  
  s<-m$estimation$likelihood$ll
  cat("LogLikelihood: ", format(round(s, 5), scientific = FALSE), "\n")
  s<-m$estimation$likelihood$aic
  cat("AIC: ", format(round(s, 5), scientific = FALSE), "\n\n")
  
  cat("Regression:\n")
  regs<-do.call("rbind", lapply(m$description$variables, function(z){z$coeff}))
  xregs<-cbind(regs, stde=NA, t=NA, pvalue=NA)
  stde<-sqrt(diag(m$estimation$bvar))
  sel<-xregs$type=='ESTIMATED'
  t<-xregs$value[sel]/stde
  ndf<-m$estimation$likelihood$nobs-m$estimation$likelihood$ndiffuse-m$estimation$likelihood$nparams+1
  pval<-2*pt(abs(t), ndf, lower.tail = F)
  xregs$stde[sel]<-stde
  xregs$t[sel]<-t
  xregs$pvalue[sel]<-pval
  print(xregs[-2])
  
}