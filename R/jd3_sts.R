#' @include jd3_ts.R jd3_rslts.R protobuf.R
#' @checkmate
NULL

#' Structural time series
#'
#' @param y 
#' @param level 
#' @param slope 
#' @param cycle 
#' @param noise 
#' @param seasonal 
#'
#' @return
#' @export
#'
#' @examples
sts<-function(y, X=NULL, level=1, slope=1, cycle=-1, noise=1
              , seasonal=c("Trigonometric", "Dummy", "Crude", "HarrisonStevens", "Fixed", "Unused"), tol=1e-12){
  
  if (!is.ts(y)){
    stop("y must be a time series")
  }
  seasonal<-match.arg(seasonal)
  jsts<-.jcall("demetra/sts/r/Bsm", "Ldemetra/sts/BasicStructuralModel;", "process", ts_r2jd(y), matrix_r2jd(X),
              as.integer(level), as.integer(slope), as.integer(cycle), as.integer(noise), seasonal, tol)
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
    estimation=p2r_sts_estimation(p$estimation)),
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
print.JDSTS<-function(m){
  cat("Structural time series", "\n\n")
  cat("Variances:\n")
  s<-m$model$level
  if (! is.na(s) && s >=0) cat("level: ", format(round(s, 6), scientific = FALSE), "\n")
  s<-m$model$slope
  if (! is.na(s) && s >=0) cat("slope: ", format(round(s, 6), scientific = FALSE), "\n")
  s<-m$model$seas
  if (! is.na(s) && s >=0) cat("seas: ", format(round(s, 6), scientific = FALSE), "\n")
  s<-m$model$n
  if (! is.na(s) && s >=0) cat("noise: ", format(round(s, 6), scientific = FALSE), "\n\n")
  s<-m$likelihood$ll
  cat("LogLikelihood: ", format(round(s, 5), scientific = FALSE), "\n")
  s<-m$estimation$score
  cat("Scores: ", format(round(s, 5), scientific = FALSE), "\n")
  
  #      ll<-proc_numeric(object@internal,"likelihood.ll")
  #      cat("Log likelihood = ", format(round(ll, 4), scientific = FALSE), "\n")
}