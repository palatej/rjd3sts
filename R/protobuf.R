enum_extract<-function(type, p){
  name<-type$value(number=p)$name()
  return (substring(name, regexpr("_", name)+1))
}

enum_of<-function(type, code, prefix){
  if (is.null(code)){
    return (as.integer(0))
  }
  i<-type$value(name=paste(prefix, code, sep='_'))$number()
}

p2r_diffuselikelihood<-function(p){
  return (structure(list(nobs=p$nobs, ndiffuse=p$ndiffuse, nparams=p$nparams, ndf=p$degrees_of_freedom,
                         ll=p$log_likelihood, adjll=p$adjusted_log_likelihood,
                         aic=p$aic, aicc=p$aicc, bic=p$bic, ssq=p$ssq, ldet=p$ldet, dcorr=p$dcorrection),
                    class = "JD3DIFFUSELIKELIHOOD"))
}

p2r_matrix<-function(p){
  m<-matrix(data=p$values, nrow = p$nrows, ncol = p$ncols)
  `attr<-`(m, "name", p$name)
  return (m)
}

p2r_ts<-function(p){
  if (length(p$values) == 0)
    return (NULL)
  s<-ts(data=p$values, frequency = p$annual_frequency, start = c(p$start_year, p$start_period))
  `attr<-`(s, "name", p$name)
  return (s)
}

p2r_test<-function(p){
  return (rjd3toolkit::statisticaltest(p$value, p$pvalue, p$description))
}

p2r_parameter<-function(p){
  if (! p$has("type")) return (NULL)
  return (list(value = p$value, type=enum_extract(jd3.ParameterType, p$type)))
}

p2r_parameters_rslt<-function(p){
  if (is.null(p))
    return (NULL)
  if (length(p) == 0)
    return (NULL)
  value<-sapply(p, function(z){z$value})
  type<-sapply(p, function(z){enum_extract(jd3.ParameterType, z$type)})
   return (data.frame(value=value, type=type))
}

p2r_parameters_rsltx<-function(p){
  if (is.null(p))
    return (NULL)
  if (length(p) == 0)
    return (NULL)
  value<-sapply(p, function(z){z$value})
  type<-sapply(p, function(z){enum_extract(jd3.ParameterType, z$type)})
  description<-sapply(p, function(z){z$description})
  
  rslt<-data.frame(value=value, type=type)
  row.names(rslt)<-description
  
  return (rslt)
}

p2r_parameters_estimation<-function(p){
  if (is.null(p))
    return (NULL)
  return (list(val=p$value, score=p$score, cov=p2r_matrix(p$covariance), description=p$description))
}

p2r_variables<-function(p){
  return (lapply(p, function(v){p2r_variable(v)}))
}

p2r_variable<-function(p){
  name<-p$name
  type<-enum_extract(modelling.VariableType, p$var_type)
  coeff<-p2r_parameters_rsltx(p$coefficients)
  return (list(name=name, type=type, coeff=coeff))
}
