#' @include utils.R
NULL

#' Title
#'
#' @param y 
#' @param period 
#' @param level -1 = no level, 0 = fixed level, 1 = sotchastic level
#' @param slope 
#' @param noise 
#' @param seasonal Seasonal model
#' @param X Regression variables (same length as y) or NULL
#' @param X.td Specification for trading days clustering. 
#' Contains thr group id fo Mondays... Sundays 
#' (for instance (1,1,1,1,1,0,0) for week days or (1,1,1,1,1,2,0) for week.saturdays/Sundays variables).
#' Contrasts are used. Can be NULL 
#'
#' @return
#' @export
#'
#' @examples
seasonalbreaks<-function(y, period=NA, level=1, slope=1, noise=1, seasonal=c("HarrisonStevens", "Trigonometric", "Dummy", "Crude", "Fixed", "Unused"),
                       X=NULL,X.td=NULL){
  
  data<-as.numeric(y)
  if (is.ts(y)){
    period<-frequency(y)
  }else{
    if (! is.null(X.td)){
      stop("y must be a time series when X.td is used")
    }
    if (is.na(period)){
      stop("y must be a time series or period must be specified")
    }
  }
  seasonal<-match.arg(seasonal)
  if (! is.null(X.td)){
    td<-rjd3modelling::td.forTs(y, X.td)
    X<-cbind(X, td)
  }

  so<-.jcall("demetra/sts/r/StsOutliersDetection", "[D", "seasonalBreaks", data, as.integer(period), 
               as.integer(level), as.integer(slope), as.integer(noise), seasonal, rjd3toolkit:::matrix_r2jd(X))
  
  if (is.ts(y)){
    return (ts(so, frequency = period, start=start(y)))
  }else{
    return (so)
  }
}
  
