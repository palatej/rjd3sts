p2r_diffuselikelihood<-function(p){
  return (structure(list(nobs=p$nobs, ndiffuse=p$ndiffuse, nparams=p$nparams, ndf=p$degrees_of_freedom,
                         ll=p$log_likelihood, adjll=p$adjusted_log_likelihood,
                         aic=p$aic, aicc=p$aicc, bic=p$bic, ssq=p$ssq, ldet=p$ldet, dcorr=p$dcorrection),
                    class = "JD3DIFFUSELIKELIHOOD"))
}

proc_diffuselikelihood<-function(jrslt, prefix){
  return (list(
    ll=.JD3_ENV$proc_numeric(jrslt, paste(prefix,"ll", sep="")),
    adjustedll=.JD3_ENV$proc_numeric(jrslt, paste(prefix,"adjustedll", sep="")),
    ssq=.JD3_ENV$proc_numeric(jrslt, paste(prefix,"ssqerr", sep="")),
    nobs=.JD3_ENV$proc_int(jrslt, paste(prefix,"nobs", sep="")),
    ndiffuse=.JD3_ENV$proc_int(jrslt, paste(prefix,"ndiffuse", sep="")),
    nparams=.JD3_ENV$proc_int(jrslt, paste(prefix,"nparams", sep="")),
    df=.JD3_ENV$proc_int(jrslt, paste(prefix,"df", sep="")))
  )
}

