#' @title Augmented Backward Elimination
#'
#' @description  Function \code{abe} performs Augmented backward elimination where variable selection is based on the change-in-estimate and significance or information criteria.
#' It can also make a backward-selection based on significance or information criteria only by turning off the change-in-estimate criterion.
#'
#'
#' @param fit An object of a class \code{"lm"}, \code{"glm"} or \code{"coxph"} representing the fit.
#' Note, the functions should be fitted with argument \code{x=TRUE} and \code{y=TRUE}.
#' @param data data frame used when fitting the object \code{fit}.
#' @param include a vector containing the names of variables that will be included in the final model. These variables are used as only passive variables during modeling. These variables might be exposure variables of interest or known confounders.
#' They will never be dropped from the working model in the selection process,
#' but they will be used passively in evaluating change-in-estimate criteria of other variables.
#' Note, variables which are not specified as include or active in the model fit are assumed to be active and passive variables.
#' @param active a vector containing the names of active variables. These less important explanatory variables will only be used as active,
#' but not as passive variables when evaluating the change-in-estimate criterion.
#' @param tau  Value that specifies the threshold of the relative change-in-estimate criterion. Default is set to 0.05.
#' @param exp.beta Logical specifying if exponent is used in formula to standardize the criterion. Default is set to TRUE.
#' @param exact Logical, specifies if the method will use exact change-in-estimate or its approximation. Default is set to FALSE, which means that the method will use approximation proposed by Dunkler et al.
#' Note, setting to TRUE can severely slow down the algorithm, but setting to FALSE can in some cases lead to a poor approximation of the change-in-estimate criterion.
#' @param criterion String that specifies the strategy to select variables for the black list.
#' Currently supported options are significance level \code{'alpha'}, Akaike information criterion \code{'AIC'} and Bayesian information criterion \code{'BIC'}.
#' If you are using significance level, in that case you have to specify the value of 'alpha' (see parameter \code{alpha}) and type of the test statistic (see parameter \code{type.test}). Default is set to \code{"alpha"}.
#' @param alpha Value that specifies the level of significance as explained above. Default is set to 0.2.
#' @param type.test String that specifies which test should be performed in case the \code{criterion = "alpha"}.
#' Possible values are \code{"F"} and \code{"Chisq"} (default) for class \code{"lm"}, \code{"Rao"}, \code{"LRT"}, \code{"Chisq"} (default), \code{"F"} for class \code{"glm"} and \code{"Chisq"} for class \code{"coxph"}. See also \code{\link{drop1}}.
#' @param type.factor String that specifies how to treat factors, see details, possible values are \code{"factor"} and \code{"individual"}.
#' @param verbose Logical that specifies if the variable selection process should be printed. Note: this can severely slow down the algorithm.
#'
#' @details
#' Using the default settings ABE will perform augmented backward elimination based on significance.
#' The level of significance will be set to 0.2. All variables will be treated as "passive or active".
#' Approximated change-in-estimate will be used. Threshold of the relative change-in-estimate criterion will be 0.05.
#' Setting tau to a very large number (e.g. \code{Inf}) turns off the change-in-estimate criterion, and ABE will only perform backward elimination.
#' Specifying \code{"alpha" = 0} will include variables only because of the change-in-estimate criterion,
#' as then variables are not safe from exclusion because of their p-values.
#' Specifying \code{"alpha" = 1} will always include all variables.
#'
#' When using \code{type.factor="individual"} each dummy variable of a factor is treated as an individual explanatory variable, hence only this dummy variable can be removed from the model (warning: use sensible coding for the reference group).
#' Using \code{type.factor="factor"} will look at the significance of removing all dummy variables of the factor and can drop the entire variable from the model.
#'
#' @return An object of class \code{"lm"}, \code{"glm"} or \code{"coxph"} representing the model chosen by abe method.
#' @references Daniela Dunkler, Max Plischke, Karen Lefondre, and Georg Heinze. Augmented backward elimination: a pragmatic and purposeful way to develop statistical models. PloS one, 9(11):e113677, 2014.
#'
#' @seealso \code{\link{abe.resampling}}, \code{\link{lm}}, \code{\link{glm}} and \code{\link{coxph}}
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @author Sladana Babic
#' @export
#' @examples
#' # simulate some data:
#'
#' set.seed(1)
#' n=100
#' x1<-runif(n)
#' x2<-runif(n)
#' x3<-runif(n)
#' y<--5+5*x1+5*x2+ rnorm(n,sd=5)
#' dd<-data.frame(y,x1,x2,x3)
#'
#' # fit a simple model containing only numeric covariates
#' fit<-lm(y~x1+x2+x3,x=TRUE,y=TRUE,data=dd)
#'
#' # perform ABE with "x1" as only passive and "x2" as only active
#' # using the exact change in the estimate of 5% and significance
#' # using 0.2 as a threshold
#' abe.fit<-abe(fit,data=dd,include="x1",active="x2",
#' tau=0.05,exp.beta=FALSE,exact=TRUE,criterion="alpha",alpha=0.2,
#' type.test="Chisq",verbose=TRUE)
#'
#' summary(abe.fit)
#'
#' # similar example, but turn off the change-in-estimate and perform
#' # only backward elimination
#'
#' abe.fit<-abe(fit,data=dd,include="x1",active="x2",
#' tau=Inf,exp.beta=FALSE,exact=TRUE,criterion="alpha",alpha=0.2,
#' type.test="Chisq",verbose=TRUE)
#'
#' summary(abe.fit)
#'
#' # an example with the model containing categorical covariates:
#' dd$x3<-rbinom(n,size=3,prob=1/3)
#' dd$y1<--5+5*x1+5*x2+ rnorm(n,sd=5)
#' fit<-lm(y1~x1+x2+factor(x3),x=TRUE,y=TRUE,data=dd)
#'
#' # treat "x3" as a single covariate:
#'
#' abe.fit.fact<-abe(fit,data=dd,include="x1",active="x2",
#' tau=0.05,exp.beta=FALSE,exact=TRUE,criterion="alpha",alpha=0.2,
#' type.test="Chisq",verbose=TRUE,type.factor="factor")
#'
#' summary(abe.fit.fact)
#'
#' # treat each dummy of "x3" as a separate covariate:
#'
#' abe.fit.ind<-abe(fit,data=dd,include="x1",active="x2",
#' tau=0.05,exp.beta=FALSE,exact=TRUE,criterion="alpha",alpha=0.2,
#' type.test="Chisq",verbose=TRUE,type.factor="individual")
#'
#' summary(abe.fit.ind)


abe<-function(fit,data=NULL,include=NULL,active=NULL,tau=0.05,exp.beta=TRUE,exact=FALSE,criterion="alpha",alpha=0.2,type.test="Chisq",type.factor=NULL,verbose=T){
  if (is.null(data)) stop("Supply the data which were used when fitting the full model.")
assign(as.character(substitute(data)),data)

if (!"x"%in%names(fit)) stop("the model should be fitted with: x=T")
if (nrow(fit$x)!=nrow(data)) stop("Data object contains missing values. Remove all the missing values and refit the model.")

if (class(fit)[1]=="lm") if (!"y"%in%names(fit)) stop("the model should be fitted with: y=T")

if (class(fit)[1]=="coxph"&exp.beta==F) stop("setting exp.beta=F for the cox model is not supported")

if (!class(fit)[1]%in%c("lm","glm","coxph","brglmFit")) stop("this model is not supported")

if (class(fit)[1]=="lm"&exp.beta==T) warning("using change in estimate for exp(b) with linear model, try to use exp.beta=F")

if (sum(unlist(lapply(strsplit(colnames(model.matrix(fit)),split=":"),function(x) length(x)!=1)))!=0) stop("interaction effects are not supported")

if (length(criterion)!=1) stop("you need to specify a single criterion")

if (colnames(model.matrix(fit))[1]=="(Intercept)") xm<-as.matrix(fit$x)[,-1] else xm<-as.matrix(fit$x)

if (!is.matrix(xm) ) stop("performing variable selection with a single variable in the model is meaningless")
if (sum(criterion%in%c("alpha","AIC","BIC"))==0) stop("valid criteria are alpha, AIC and BIC")
if (criterion=="alpha") if (alpha<0|alpha>1) stop("specify alpha between zero and one")
if (is.null(tau)) stop("Specify tau.")
if (tau<0) stop("Tau has to be >=0.")
if (sum(my_grepl("offset",names(attributes(fit$terms)$dataClasses)))!=0){

  warning("offset variables are in the model treating them as only passive")

  offset.var<-names(attributes(fit$terms)$dataClasses)[my_grepl("offset",names(attributes(fit$terms)$dataClasses))]

  if (!is.null(active)){
    for (ii in 1:length(active)){
      if (sum(my_grepl(  active[ii],offset.var )!=0)) active[ii]<-NA
    }

  active=active[!is.na(active)]
  if (length(active)==0) active=NULL
  }
  if (!is.null(include)){
    for (ii in 1:length(include)){
      if (sum(my_grepl(  include[ii],offset.var )!=0)) include[ii]<-NA
    }

  include=include[!is.na(include)]
  if (length(include)==0) include=NULL
  }

}

if (class(fit)[1]=="coxph"){
  if (sum(my_grepl("strata",attributes(fit$terms)$term.labels))!=0){
    strata.vars<-attributes(fit$terms)$term.labels[my_grepl("strata",attributes(fit$terms)$term.labels)]
    if ( is.null(include)&is.null(active) ) {
      if (criterion=="alpha") stop("The model includes strata and stratification variables are either only active or passive and active, using alpha as a criterion is inappropriate, use either AIC or BIC.")
      warning("The model includes strata and stratification variables are either only active or passive and active, using approximate change-in-estimate is inappropriate, using exact change-in-estimate.")
      exact=TRUE
      }
    if (!is.null(active)&is.null(include)){

        if (criterion=="alpha") stop("The model includes strata and stratification variables are either only active or passive and active, using alpha as a criterion is inappropriate, use either AIC or BIC.")
        warning("The model includes strata and stratification variables are either only active or passive and active, using approximate change-in-estimate is inappropriate, using exact change-in-estimate.")
      exact=TRUE


    }
    if (!is.null(active)&!is.null(include)){
      log<-rep(NA,length(active))
      for (ii in 1:length(active)){
        log[ii]<-sum(my_grepl(active[ii],strata.vars))
      }
      log1<-rep(NA,length(include))
      for (ii in 1:length(include)){
        log1[ii]<-sum(my_grepl(include[ii],strata.vars))
      }
      if (sum(log)!=0&sum(log1)!=0) stop("stratification variable cannot be specified as include and active")
      if (sum(log1)==0){
        if (criterion=="alpha") stop("The model includes strata and stratification variables are either only active or passive and active, using alpha as a criterion is inappropriate, use either AIC or BIC.")
        warning("The model includes strata and stratification variables are either only active or passive and active, using approximate change-in-estimate is inappropriate, using exact change-in-estimate.")
        exact=TRUE

      }
    }

    if (!is.null(include)&is.null(active)){
      log1<-rep(NA,length(include))
      for (ii in 1:length(include)){
        log1[ii]<-sum(my_grepl(include[ii],strata.vars))
      }
      if (sum(log1)==0){
        if (criterion=="alpha") stop("The model includes strata and stratification variables are either only active or passive and active, using alpha as a criterion is inappropriate, use either AIC or BIC.")
        warning("The model includes strata and stratification variables are either only active or passive and active, using approximate change-in-estimate is inappropriate, using exact change-in-estimate.")
        exact=TRUE

      }
    }

        }

}

if (length( my_grep("matrix",attributes(fit$terms)$dataClasses[-1]))!=0){ warning("The model contains a covariate which was fitted as matrix (splines etc.), treating it as a single covariate. Using type.factor=factor for factors. Using approximate change-in-estimate is inappropriate; using exact change in estimate instead.");bt<-abe.fact1(fit,data,include,active,tau,exp.beta,exact=TRUE,criterion,alpha,type.test,verbose)} else {

if (sum(attributes(fit$terms)$dataClasses[!my_grepl("strata",names(attributes(fit$terms)$dataClasses))]=="factor")>0)  {

if (is.null(type.factor)) {type.factor="factor"; warning("There are factors in the model but type.factor is not specified, using type.factor=factor")}

	if (type.factor=="factor") bt<-abe.fact1(fit,data,include,active,tau,exp.beta,exact,criterion,alpha,type.test,verbose) else bt<-abe.fact2(fit,data,include,active,tau,exp.beta,exact,criterion,alpha,type.test,verbose)


	} else  bt<-abe.num(fit,data,include,active,tau,exp.beta,exact,criterion,alpha,type.test,verbose)
}

bt
}







#' Resampled Augmented Backward Elimination
#'
#' Performs Augmented backward elimination on re-sampled datasets using different bootstrap and re-sampling techniques.
#'
#'
#' @param fit An object of a class \code{"lm"}, \code{"glm"} or \code{"coxph"} representing the fit.
#' Note, the functions should be fitted with argument \code{x=TRUE} and \code{y=TRUE}.
#' @param data data frame used when fitting the object \code{fit}.
#' @param include a vector containing the names of variables that will be included in the final model. These variables are used as passive variables during modeling. These variables might be exposure variables of interest or known confounders.
#' They will never be dropped from the working model in the selection process,
#' but they will be used passively in evaluating change-in-estimate criteria of other variables.
#' Note, variables which are not specified as include or active in the model fit are assumed to be active and passive variables.
#' @param active a vector containing the names of active variables. These less important explanatory variables will only be used as active,
#' but not as passive variables when evaluating the change-in-estimate criterion.
#' @param tau  Value that specifies the threshold of the relative change-in-estimate criterion. Default is set to 0.05.
#' @param exp.beta Logical specifying if exponent is used in formula to standardize the criterion. Default is set to TRUE.
#' @param exact Logical, specifies if the method will use exact change-in-estimate or approximated. Default is set to FALSE, which means that the method will use approximation proposed by Dunkler et al.
#' Note, setting to TRUE can severely slow down the algorithm, but setting to FALSE can in some cases lead to a poor approximation of the change-in-estimate criterion.
#' @param criterion String that specifies the strategy to select variables for the blacklist.
#' Currently supported options are significance level \code{'alpha'}, Akaike information criterion \code{'AIC'} and Bayesian information criterion \code{'BIC'}.
#' If you are using significance level, in that case you have to specify the value of 'alpha' (see parameter \code{alpha}). Default is set to \code{"alpha"}.
#' @param alpha Value that specifies the level of significance as explained above. Default is set to 0.2.
#' @param type.test String that specifies which test should be performed in case the \code{criterion = "alpha"}.
#' Possible values are \code{"F"} and \code{"Chisq"} (default) for class \code{"lm"}, \code{"Rao"}, \code{"LRT"}, \code{"Chisq"} (default), \code{"F"} for class \code{"glm"} and \code{"Chisq"} for class \code{"coxph"}. See also \code{\link{drop1}}.
#' @param type.factor String that specifies how to treat factors, see details, possible values are \code{"factor"} and \code{"individual"}.
#' @param num.boot number of bootstrap re-samples
#' @param type.boot String that specifies the type of bootstrap. Possible values are \code{"bootstrap"}, \code{"mn.bootstrap"}, \code{"subsampling"},  see details
#' @param prop.sampling Sampling proportion. Only applicable for \code{type.boot="mn.bootstrap"} and \code{type.boot="subsampling"}, defaults to 0.5. See details.
#' @return an object of class \code{abe} for which \code{summary}, \code{plot} and \code{pie.abe} functions are available.
#' A list with the following elements:
#'
#' \code{models} the final models obtained after performing ABE on re-sampled datasets, each object in the list is of the same class as \code{fit}
#'
#' \code{alpha} the vector of significance levels used
#'
#' \code{tau} the vector of threshold values for the change-in-estimate
#'
#' \code{num.boot} number of re-sampled datasets
#'
#' \code{criterion} criterion used when constructing the black-list
#'
#' \code{all.vars} a list of variables used when estimating \code{fit}
#'
#' \code{fit.or} the initial model
#'
#' \code{misc} the parameters of the call to \code{abe.resampling}
#'
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @author Sladana Babic
#' @details \code{type.boot} can be \code{bootstrap} (n observations drawn from the original data with replacement), \code{mn.bootstrap} (m out of n observations drawn from the original data with replacement), \code{subsampling} (m out of n observations drawn from the original data without replacement), where m is [prop.sampling*n].
#' @references Daniela Dunkler, Max Plischke, Karen Lefondre, and Georg Heinze. Augmented backward elimination: a pragmatic and purposeful way to develop statistical models. PloS one, 9(11):e113677, 2014.
#' @references Riccardo De Bin, Silke Janitza, Willi Sauerbrei and Anne-Laure Boulesteix. Subsampling versus Bootstrapping in Resampling-Based Model Selection for Multivariable Regression. Biometrics 72, 272-280, 2016.
#' @seealso \code{\link{abe}}, \code{\link{summary.abe}}, \code{\link{print.abe}}, \code{\link{print.abe.default}}, \code{\link{plot.abe}}, \code{\link{pie.abe}}
#' @export
#' @examples
#' # simulate some data and fit a model
#'
#' set.seed(1)
#' n=100
#' x1<-runif(n)
#' x2<-runif(n)
#' x3<-runif(n)
#' y<--5+5*x1+5*x2+ rnorm(n,sd=5)
#' dd<-data.frame(y=y,x1=x1,x2=x2,x3=x3)
#' fit<-lm(y~x1+x2+x3,x=TRUE,y=TRUE,data=dd)
#'
#' # use ABE on 50 bootstrap re-samples considering different
#' # change-in-estimate thresholds and significance levels
#'
#' fit.boot<-abe.resampling(fit,data=dd,include="x1",active="x2",
#' tau=c(0.05,0.1),exp.beta=FALSE,exact=TRUE,
#' criterion="alpha",alpha=c(0.2,0.05),type.test="Chisq",
#' num.boot=50,type.boot="bootstrap")
#'
#' summary(fit.boot)
#'
#' # use ABE on 50 subsamples randomly selecting 50% of subjects
#' # considering different change-in-estimate thresholds and
#' # significance levels
#'
#' fit.resample<-abe.resampling(fit,data=dd,include="x1",active="x2",
#' tau=c(0.05,0.1),exp.beta=FALSE,exact=TRUE,
#' criterion="alpha",alpha=c(0.2,0.05),type.test="Chisq",
#' num.boot=50,type.boot="subsampling",prop.sampling=0.5)
#'
#' summary(fit.resample)

abe.resampling<-function(fit,data=NULL,include=NULL,active=NULL,tau=0.05,exp.beta=TRUE,exact=FALSE,criterion="alpha",alpha=0.2,type.test="Chisq",type.factor=NULL,num.boot=100,type.boot=c("bootstrap","mn.bootstrap","subsampling"),prop.sampling=0.5){

  if (is.null(data)) stop("Supply the data which were used when fitting the full model.")

  if (!"x"%in%names(fit)) stop("the model should be fitted with: x=T")
  if (nrow(fit$x)!=nrow(data)) stop("Data contains missing values. Remove all the missing values and refit the model.")

  if (class(fit)[1]=="lm") if (!"y"%in%names(fit)) stop("the model should be fitted with: y=T")

  if (class(fit)[1]=="coxph"&exp.beta==F) stop("setting exp.beta=F for the cox model is not supported")

  if (!class(fit)[1]%in%c("lm","glm","coxph","brglmFit")) stop("this model is not supported")

  if (class(fit)[1]=="lm"&exp.beta==T) warning("using change in estimate for exp(b) with linear model, try to use exp.beta=F")

  if (sum(unlist(lapply(strsplit(colnames(model.matrix(fit)),split=":"),function(x) length(x)!=1)))!=0) stop("interaction effects are not supported")

  if (length(criterion)!=1) stop("you need to specify a single criterion")
  if (sum(criterion%in%c("alpha","AIC","BIC"))==0) stop("valid criteria are alpha, AIC and BIC")
  if (criterion=="alpha") if (sum(alpha<0)!=0|sum(alpha>1)!=0) stop("specify alphas between zero and one")
  if (is.null(tau)) stop("Specify tau.")
  if (sum(tau<0)!=0) stop("Taus has to be >=0.")

  if (length(type.boot)!=1) stop("you need to specify a single resampling method")

  if (criterion!="alpha") alpha=NULL

  if (colnames(model.matrix(fit))[1]=="(Intercept)") xm<-as.matrix(fit$x)[,-1] else xm<-as.matrix(fit$x)

  if (!is.matrix(xm) ) stop("performing variable selection with a single variable in the model is meaningless")

  if (sum(my_grepl("offset",names(attributes(fit$terms)$dataClasses)))!=0){

    warning("offset variables are in the model treating them as only passive")

    offset.var<-names(attributes(fit$terms)$dataClasses)[my_grepl("offset",names(attributes(fit$terms)$dataClasses))]

    if (!is.null(active)){
      for (ii in 1:length(active)){
        if (sum(my_grepl(  active[ii],offset.var )!=0)) active[ii]<-NA
      }

      active=active[!is.na(active)]
      if (length(active)==0) active=NULL
    }
    if (!is.null(include)){
      for (ii in 1:length(include)){
        if (sum(my_grepl(  include[ii],offset.var )!=0)) include[ii]<-NA
      }

      include=include[!is.na(include)]
      if (length(include)==0) include=NULL
    }

  }
  if (class(fit)[1]=="coxph"){
    if (sum(my_grepl("strata",attributes(fit$terms)$term.labels))!=0){
      strata.vars<-attributes(fit$terms)$term.labels[my_grepl("strata",attributes(fit$terms)$term.labels)]
      if ( is.null(include)&is.null(active) ) {
        if (criterion=="alpha") stop("The model includes strata and stratification variables are either only active or passive and active, using alpha as a criterion is inappropriate, use either AIC or BIC.")
        warning("The model includes strata and stratification variables are either only active or passive and active, using approximate change-in-estimate is inappropriate, using exact change-in-estimate.")
        exact=TRUE
      }
      if (!is.null(active)&is.null(include)){

        if (criterion=="alpha") stop("The model includes strata and stratification variables are either only active or passive and active, using alpha as a criterion is inappropriate, use either AIC or BIC.")
        warning("The model includes strata and stratification variables are either only active or passive and active, using approximate change-in-estimate is inappropriate, using exact change-in-estimate.")
        exact=TRUE


      }
      if (!is.null(active)&!is.null(include)){
        log<-rep(NA,length(active))
        for (ii in 1:length(active)){
          log[ii]<-sum(my_grepl(active[ii],strata.vars))
        }
        log1<-rep(NA,length(include))
        for (ii in 1:length(include)){
          log1[ii]<-sum(my_grepl(include[ii],strata.vars))
        }
        if (sum(log)!=0&sum(log1)!=0) stop("stratification variable cannot be specified as include and active")
        if (sum(log1)==0){
          if (criterion=="alpha") stop("The model includes strata and stratification variables are either only active or passive and active, using alpha as a criterion is inappropriate, use either AIC or BIC.")
          warning("The model includes strata and stratification variables are either only active or passive and active, using approximate change-in-estimate is inappropriate, using exact change-in-estimate.")
          exact=TRUE

        }
      }

      if (!is.null(include)&is.null(active)){
        log1<-rep(NA,length(include))
        for (ii in 1:length(include)){
          log1[ii]<-sum(my_grepl(include[ii],strata.vars))
        }
        if (sum(log1)==0){
          if (criterion=="alpha") stop("The model includes strata and stratification variables are either only active or passive and active, using alpha as a criterion is inappropriate, use either AIC or BIC.")
          warning("The model includes strata and stratification variables are either only active or passive and active, using approximate change-in-estimate is inappropriate, using exact change-in-estimate.")
          exact=TRUE

        }
      }

    }

  }
  if (length( my_grep("matrix",attributes(fit$terms)$dataClasses[-1]))!=0){warning("The model contains a covariate which was fitted as a matrix (splines etc.), treating it as a single covariate. Using type.factor=factor for factors. Using approximate change-in-estimate is inappropriate; using exact change in estimate instead.")} else {


  if (sum(attributes(fit$terms)$dataClasses[!my_grepl("strata",names(attributes(fit$terms)$dataClasses))]=="factor")>0&is.null(type.factor)) {type.factor="factor"; warning("There are factors in the model but type.factor is not specified, using type.factor=factor")}

  if (sum(attributes(fit$terms)$dataClasses[!my_grepl("strata",names(attributes(fit$terms)$dataClasses))]=="factor")>0) if(type.factor=="factor"&exact==FALSE) {warning("there are factors in the model, using approximate change-in-estimate with this type.factor is inappropriate; using exact change in estimate instead"); exact=T}
  }

  cnms<-attributes(fit$terms)$term.labels

  if (!is.null(include)) {
    include.l<-list()
    for (i in 1:length(include)) {
      if (sum(my_grepl(include[i],cnms))==0) stop("at least one include variable is not in the model")
      include.l[[i]]<-cnms[my_grep(include[i],cnms)]
    }
    include2<-unlist(include.l)
  } else include2<-include

  if (!is.null(active))  {
    active.l<-list()
    for (i in 1:length(active)) {
      if (sum(my_grepl(active[i],cnms))==0) stop("at least one active variable is not in the model")
      active.l[[i]]<-cnms[my_grep(active[i],cnms)]
    }
    active2<-unlist(active.l)
  } else active2<-active


   if (sum(include2%in%active2)!=0) stop("at least one include variable is also specified as active")



  if ( sum(cnms%in%include2)==length(cnms) ) stop("all variables are specified as pasive, cannot perform variable selection")
  if ( sum(cnms%in%active2)==length(cnms) ) {
    include=active=NULL
    warning("all variables are specified as active, treating all variables as active or pasive")
  }

  if (criterion[1]=="alpha"&is.null(alpha)) stop("specify alpha")


  if (criterion[1]=="alpha") k<-NULL
  if (criterion[1]=="AIC") k<-2

  n<-nrow(fit$x)
  if (type.boot!="bootstrap") m<-round(n*prop.sampling,0) else m<-n
  if (criterion[1]=="BIC") k<-log(m)



  boot<-list()

  i=0

  if (!is.null(alpha)&!is.null(tau)){

  for (a in alpha){

    if (criterion[1]=="alpha") k<-qchisq(1-a,df=1)

  for (t in tau){
    for (ii in 1:num.boot){
     i=i+1

     if (type.boot=="bootstrap") ids<-sample(1:n,n,replace=T)
     if (type.boot=="mn.bootstrap") ids<-sample(1:n,m,replace=T)
     if (type.boot=="subsampling") ids<-sample(1:n,m,replace=F)

    data.boot<-data[ids,]

     fit.i<-my_update_boot(fit,data=data.boot)
    if (length( my_grep("matrix",attributes(fit$terms)$dataClasses[-1]))!=0){boot[[i]]<-abe.fact1.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=a,type.test,k) } else {
    if (sum(attributes(fit$terms)$dataClasses[!my_grepl("strata",names(attributes(fit$terms)$dataClasses))]=="factor")>0)  {

    if (type.factor=="factor") {


                boot[[i]]<-abe.fact1.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=a,type.test,k) } else  {

                boot[[i]]<-abe.fact2.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=a,type.test,k)#,data=data.boot)


                      }


    } else  boot[[i]]<-abe.num.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=a,type.test,k)

    }}
   }
  }}


  if (is.null(alpha)&!is.null(tau)){


      for (t in tau){
        for (ii in 1:num.boot){
          i=i+1
          if (type.boot=="bootstrap") ids<-sample(1:n,n,replace=T)
          if (type.boot=="mn.bootstrap") ids<-sample(1:n,m,replace=T)
          if (type.boot=="subsampling") ids<-sample(1:n,m,replace=F)

          data.boot<-data[ids,]

            fit.i<-my_update_boot(fit,data=data.boot)

          if (length( my_grep("matrix",attributes(fit$terms)$dataClasses[-1]))!=0){boot[[i]]<-abe.fact1.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=alpha,type.test,k) } else {

          if (sum(attributes(fit$terms)$dataClasses[!my_grepl("strata",names(attributes(fit$terms)$dataClasses))]=="factor")>0)  {

            if (type.factor=="factor") {


              boot[[i]]<-abe.fact1.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=alpha,type.test,k) } else  {

                boot[[i]]<-abe.fact2.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=alpha,type.test,k)#,data=data.boot)
                      }


          } else  boot[[i]]<-abe.num.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=alpha,type.test,k)

        }}

    }}

  if (!is.null(alpha)&is.null(tau)){


    for (a in alpha){

      if (criterion[1]=="alpha") k<-qchisq(1-a,df=1)

      for (ii in 1:num.boot){
        i=i+1
        if (type.boot=="bootstrap") ids<-sample(1:n,n,replace=T)
        if (type.boot=="mn.bootstrap") ids<-sample(1:n,m,replace=T)
        if (type.boot=="subsampling") ids<-sample(1:n,m,replace=F)

        data.boot<-data[ids,]

         fit.i<-my_update_boot(fit,data=data.boot)

        if (length( my_grep("matrix",attributes(fit$terms)$dataClasses[-1]))!=0){boot[[i]]<-abe.fact1.boot(fit.i,data.boot,include,active,tau=tau,exp.beta,exact,criterion,alpha=a,type.test,k) } else {

        if (sum(attributes(fit$terms)$dataClasses[!my_grepl("strata",names(attributes(fit$terms)$dataClasses))]=="factor")>0)  {

          if (type.factor=="factor") {


            boot[[i]]<-abe.fact1.boot(fit.i,data.boot,include,active,tau=tau,exp.beta,exact,criterion,alpha=a,type.test,k) } else  {

              boot[[i]]<-abe.fact2.boot(fit.i,data.boot,include,active,tau=tau,exp.beta,exact,criterion,alpha=a,type.test,k)#,data=data.boot)

                       }


        } else  boot[[i]]<-abe.num.boot(fit.i,data.boot,include,active,tau=tau,exp.beta,exact,criterion,alpha=a,type.test,k)

      }}

    }}

  fit.or<-fit
if (length( my_grep("matrix",attributes(fit$terms)$dataClasses[-1]))==0){
  if (sum(attributes(fit$terms)$dataClasses[!my_grepl("strata",names(attributes(fit$terms)$dataClasses))]=="factor")>0){
    if (type.factor=="individual") {
      df<-as.data.frame(model.matrix(fit))

      names(df)<-gsub("factor", replacement="", names(df), ignore.case = FALSE, perl = FALSE,
                      fixed = FALSE, useBytes = FALSE)

      names(df)<-gsub(")", replacement=".", names(df), ignore.case = FALSE, perl = FALSE,
                      fixed = FALSE, useBytes = FALSE)

      names(df)<-gsub("\\(", replacement="", names(df), ignore.case = FALSE, perl = FALSE,
                      fixed = FALSE, useBytes = FALSE)

      check.names<-names(df)
      var.mod<-attributes(fit$terms)$term.labels
      if (sum(my_grepl("strata",var.mod))!=0)  check.names<-c(check.names, var.mod[my_grepl("strata",var.mod)] )



      if (!is.null(include)) {
        include.l<-list()
        for (i in 1:length(include)) {
          if (sum(my_grepl(include[i],check.names))==0) stop("at least one include variable is not in the model")
          include.l[[i]]<-check.names[my_grep(include[i],check.names)]
        }
        include<-unlist(include.l)
      }
      if (!is.null(active))  {
        active.l<-list()
        for (i in 1:length(active)) {
          if (sum(my_grepl(active[i],check.names))==0) stop("at least one active variable is not in the model")
          active.l[[i]]<-check.names[my_grep(active[i],check.names)]
        }
        active<-unlist(active.l)
      }

        if ( colnames(model.matrix(fit))[1]=="(Intercept)" )   updt.f<-as.formula(paste("~",paste(check.names[-1],collapse="+"))) else  updt.f<-as.formula(paste("~",paste(check.names,collapse="+"),"-1"))

      df<-cbind(df,model.frame(fit),data)

      fit.or<-my_update(fit, updt.f   ,data=df)

      }
  }
}

misc<-list(tau=tau,criterion=criterion,alpha=alpha,type.boot=type.boot,prop.sampling=prop.sampling)

  res<-list(models=boot,alpha=alpha,tau=tau,num.boot=num.boot,criterion=criterion,all.vars=names(coef(fit.or)),fit.or=fit.or,misc=misc,call=match.call())

  class(res)<-"abe"
   res
}



#' Print Function
#'
#' Prints a summary table of a resampled version of ABE using the latest guideliness.
#' The table displays the coefficient estimates and standard errors from the initial model (model with all covariates),
#' the relative inclusion frequencies of the covariates from the initial model (using \code{type.boot="resampling"} and \code{prop.sampling =0.5}),
#' resampled median and percentiles for the estimates of the regression coefficients for each variable from the initial model,
#' root mean squared difference ratio (RMSD) and relative bias conditional on selection (RBCS) all using \code{type.boot="bootstrap"}.
#' While not required, it makes sense to call \code{\link{abe.resampling}} with \code{type.boot="bootstrap"}. If it is not specified in this way, the print function will override the argument \code{type.boot} and refit \code{\link{abe.resampling}} (twice) which greatly increases the computing time.
#'
#' @param x an object of class \code{"abe"}, an object returned by a call to \code{\link{abe.resampling}}, preferably using \code{type.boot="bootstrap"} or \code{type.boot="subsampling"} with \code{prop.sampling=0.5}.
#' @param conf.level the confidence level, defaults to 0.95
#' @param alpha the alpha value for which the output is to be printed, defaults to \code{NULL}
#' @param tau the tau value for which the output is to be printed, defaults to \code{NULL}
#' @param ... additional arguments affecting the summary produced.
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @author Sladana Babic
#' @seealso \code{\link{print.abe}}, \code{\link{abe.resampling}}, \code{\link{summary.abe}}, \code{\link{plot.abe}}, \code{\link{pie.abe}}
#' @export
#' @examples
#' set.seed(1)
#' n=100
#' x1<-runif(n)
#' x2<-runif(n)
#' x3<-runif(n)
#' y<--5+5*x1+5*x2+ rnorm(n,sd=5)
#' dd<-data.frame(y=y,x1=x1,x2=x2,x3=x3)
#' fit<-lm(y~x1+x2+x3,x=TRUE,y=TRUE,data=dd)
#'
#' fit.resampling<-abe.resampling(fit,data=dd,include="x1",active="x2",
#' tau=c(0.05,0.1),exp.beta=FALSE,exact=TRUE,
#' criterion="alpha",alpha=c(0.2,0.05),type.test="Chisq",
#' num.boot=50,type.boot="bootstrap")
#'
#' print.abe.default(fit.resampling,conf.level=0.95,alpha=0.2,tau=0.05)


print.abe.default<-function(x,conf.level=0.95,alpha=NULL,tau=NULL,...){
  cat("Printing the output of a call to abe.resampling using the latest guideliness. Be patient, as this might take a while.")
  object<-x
  type.boot<-object$misc$type.boot

  if (type.boot=="subsampling") {
      if  (object$misc$prop.sampling==0.5 ) subs<-object else subs<-update(object,prop.sampling =0.5 )
      boot<-update(object,type.boot = "bootstrap" )
      }
  if (type.boot=="bootstrap"){
    boot<-object
    subs<-update(object,type.boot ="subsampling",prop.sampling =0.5 )
  }
  if (type.boot=="mn.bootstrap"){
    subs<-update(object,type.boot ="subsampling",prop.sampling =0.5 )
    boot<-update(object,type.boot = "bootstrap" )
  }
  if (conf.level<0|conf.level>1) stop("Confidence level out of range")
  if (!is.null(tau)) if(tau%in%object$misc$tau==FALSE) stop("This value of tau was not used in a call to abe.resampling.")
  if (object$misc$criterion=="alpha") if(!is.null(alpha)) if(alpha%in%object$misc$alpha==FALSE) stop("This value of alpha was not used in a call to abe.resampling.")


  if (object$misc$criterion=="alpha") if(is.null(alpha)) alpha=object$misc$alpha[1]

  if (is.null(tau)) tau=object$misc$tau[1]

  cat("Printing results of a call to abe.resampling for:\n  tau=",tau,"\n  criterion=\"",object$misc$criterion,"\"",sep="")
  if (object$misc$criterion=="alpha" ) cat("  alpha=",alpha,sep="")

  set<-paste("tau=",tau,sep="")
  sea<-paste("alpha=",alpha,sep="")

  if (object$misc$criterion=="alpha") rs<-paste(set,sea,sep=".") else rs<-set
  ss.susb<-summary(subs,conf.level = conf.level)
  ss.boot<-summary(boot,conf.level = conf.level)

  ss<-summary(object,conf.level = conf.level)

  ss$var.rel.frequencies<-ss.susb$var.rel.frequencies
  ss$model.rel.frequencies<-ss.susb$model.rel.frequencies
  ss$var.coefs<-ss.boot$var.coefs

  vars.num<-c(object$all.vars,attributes(object$fit.or$terms)$term.labels[my_grepl("strata",attributes(object$fit.or$terms)$term.labels)])

  idv<-which(vars.num%in%object$all.vars==TRUE)
  idg<-which(vars.num%in%object$all.vars==FALSE)


  if (length(idv)!=length(vars.num)) {
    pt<-paste(paste(names(ss$var.rel.frequencies[rs,])[idg],": ",ss$var.rel.frequencies[rs,idg],sep=""),collapse=",")
    cat("\n\n Inclusion relative frequencies of the offset/stratification variable(s):",pt)

  }
  mat <- cbind(coef(object$fit.or),
               sqrt(diag(vcov(object$fit.or))),
               c(ss$var.rel.frequencies[rs,idv]),
               t(ss$var.coefs[[rs]][,idv]))
  rownames(mat)<-names(ss$var.rel.frequencies[rs,])[idv]

  colnames(mat)[1:3] <- c("Estimate init.", "Std. Error init.",
                          "Incl. Freq.")
  colnames(mat)[4:6]<-c("Estimate, 50%", paste("Estimate ",(1-conf.level)/2*100,"%",sep=""),paste("Estimate ",100-(1-conf.level)/2*100,"%",sep="") )
  colnames(mat)[7:8]<-c("RMSD ratio","RBCS")
  cat("\n\n")

  print(mat)
}



#' Summary Function
#'
#' makes a summary of a resampled version of ABE
#'
#' @param object an object of class \code{"abe"}, an object returned by a call to \code{\link{abe.resampling}}
#' @param conf.level the confidence level, defaults to 0.95
#' @param ... additional arguments affecting the summary produced.

#' @return a list with the following elements:
#'
#' \code{var.rel.frequencies}: inclusion relative frequencies for all variables from the initial model
#'
#' \code{model.rel.frequencies}: relative frequencies of the final models
#'
#' \code{var.coefs}: bootstrap/resampled medians and percentiles for the estimates of the regression coefficients for each variable from the initial model
#'
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @author Sladana Babic
#' @seealso \code{\link{abe.resampling}}, \code{\link{print.abe.default}}, \code{\link{print.abe}}, \code{\link{plot.abe}}, \code{\link{pie.abe}}
#' @export
#' @examples
#' set.seed(1)
#' n=100
#' x1<-runif(n)
#' x2<-runif(n)
#' x3<-runif(n)
#' y<--5+5*x1+5*x2+ rnorm(n,sd=5)
#' dd<-data.frame(y=y,x1=x1,x2=x2,x3=x3)
#' fit<-lm(y~x1+x2+x3,x=TRUE,y=TRUE,data=dd)
#'
#' fit.boot<-abe.resampling(fit,data=dd,include="x1",active="x2",
#' tau=c(0.05,0.1),exp.beta=FALSE,exact=TRUE,
#' criterion="alpha",alpha=c(0.2,0.05),type.test="Chisq",
#' num.boot=50,type.boot="bootstrap")
#'
#' summary(fit.boot)$var.rel.frequencies


summary.abe<-function(object,conf.level=0.95,...){

if (object$criterion!="alpha") alphas<-NULL else {

if (!is.null(object$tau)) alphas<-paste("alpha=",rep(object$alpha,each=length(object$tau)*object$num.boot),sep="")
if (is.null(object$tau)) alphas<-paste("alpha=",rep(object$alpha,each=1*object$num.boot),sep="")

}


if (is.null(object$tau)) taus<-NULL else {

if (is.null(object$alpha)) {taus<- paste("tau=",rep(rep( object$tau,each=object$num.boot  ),1),sep="")} else {

taus<- paste("tau=",rep(rep( object$tau,each=object$num.boot  ),length(object$alpha)),sep="")
}
}

if (!is.null(object$tau)&!is.null(object$alpha)) boot.iter<- rep(1:object$num.boot,length(object$alpha)*length(object$tau))

if (is.null(object$tau)&!is.null(object$alpha)) boot.iter<- rep(1:object$num.boot,length(object$alpha))

if (!is.null(object$tau)&is.null(object$alpha)) boot.iter<- rep(1:object$num.boot,length(object$tau))


 vars.model<-lapply(object$models,function(x) if (is.numeric(x)) "empty model" else c(names(coef(x)),attributes(x$terms)$term.labels[my_grepl("strata",attributes(x$terms)$term.labels)]))
vars.model.cf<-lapply(object$models,function(x) if (is.numeric(x)) "empty model" else names(coef(x)))


ggff<-function(x){ tbx<-table(x); {for (jj in which((tbx>1)==T)) {x[x==names(tbx[jj])]<-paste(x[x==names(tbx[jj])],1:sum(x==names(tbx[jj])),sep="")  }} ;x  }


vars.model<-lapply(vars.model,ggff  )
vars.model.cf<-lapply(vars.model.cf,ggff   )

coefs.model<-lapply(object$models,function(x) if (is.numeric(x)) NULL else  coef(x))


object$all.vars<-ggff(object$all.vars)


mm<-lapply(1:length(vars.model.cf),function(i,x,y,z) {zz<-list();zz[[i]]<-rep(0,length(y));zz[[i]][y%in%x[[i]]==T]<-z[[i]];names(zz[[i]])<-y;zz[[i]]},  vars.model.cf,object$all.vars,coefs.model  )


ff<-function(x){

gg<-matrix(unlist(x),nrow=length(x),ncol=length(object$all.vars),byrow=T)

sum<-apply(gg,2,function(x)  c(median(x,na.rm=T),quantile(x,na.rm=T,probs=c((1-conf.level)/2,conf.level+(1-conf.level)/2)))  )

colnames(sum)<-object$all.vars

rownames(sum)<-c("median","CI lower","CI upper")
sum

}



rmsd<-lapply(1:length(mm),function(i,x,y) (x[[i]]-y)**2/diag(vcov(object$fit.or))  ,mm,  coef(object$fit.or)) #here 0 is imputed when var is not in the model!
rel.bias<-lapply(1:length(mm),function(i,x,y) (x[[i]]/(y)  )  ,mm,  coef(object$fit.or)) #this inputs 0 if var j is not in the model (different than what Georg proposes!), note that later with funk I divide by VIF and substract 1!


#rmsd<-lapply(1:length(mm),function(i,x,y) (ifelse(x[[i]]==0,NA,x[[i]])-y)**2/diag(vcov(object$fit.or))  ,mm,  coef(object$fit.or)) #here the sum is only within subsets where j is included
#rel.bias<-lapply(1:length(mm),function(i,x,y) (ifelse(x[[i]]==0,NA,x[[i]])/(y)  )  ,mm,  coef(object$fit.or)) #here the sum is only within subsets where j is included



ff1<-function(x){

  if (length(x)!=1) gg<-matrix(unlist(x),nrow=length(x),ncol=length(object$all.vars),byrow=T) else gg<-matrix(unlist(x),nrow=length(unlist(x))/length(object$all.vars),ncol=length(object$all.vars),byrow=T)

  sum<-apply(gg,2,function(x)  sqrt(mean(x,na.rm=T)))

names(sum)<-object$all.vars

sum

}


ff.relb<-function(x){

  if (length(x)!=1) gg<-matrix(unlist(x),nrow=length(x),ncol=length(object$all.vars),byrow=T) else gg<-matrix(unlist(x),nrow=length(unlist(x))/length(object$all.vars),ncol=length(object$all.vars),byrow=T)

  sum<-apply(gg,2,function(x)  mean(x,na.rm=T))

  names(sum)<-object$all.vars

  sum

}

vars.model.col<-lapply(vars.model,function(x) paste(x,collapse="+"))
vars.num<-c(object$all.vars,attributes(object$fit.or$terms)$term.labels[my_grepl("strata",attributes(object$fit.or$terms)$term.labels)])

funk<-function(i,x,y) {(x[[i]]/y[i,]-1)*100}

if (object$criterion=="alpha"&!is.null(object$tau)) {
    ss<-lapply(split(vars.model,list(taus,alphas)),function(x) table(unlist(x))[vars.num]/object$num.boot)
    ss.col<-lapply(split(vars.model.col,list(taus,alphas)),function(x) table(unlist(x))[order(-table(unlist(x)))]/object$num.boot)

    gg<-lapply(ss,function(x) {z<-ifelse(vars.num%in%names(x)==F,0,x);names(z)<-vars.num;z} )

    gg1<-matrix(unlist(gg),byrow=T,ncol=length(vars.num))
    colnames(gg1)<-vars.num
    rownames(gg1)<-names(ss)

    ss1<-lapply(split(mm,list(taus,alphas)),ff)
    ss2<-lapply(split(rmsd,list(taus,alphas)),ff1)

    ss1<-lapply(1:length(ss1),function(i,x,y) {rr<-rbind(x[[i]],y[[i]]);colnames(rr)<-colnames(x[[i]]); rownames(rr)<-c(rownames(x[[i]]),"RMSD ratio") ;rr   }, ss1,ss2)
    names(ss1)<-names(ss2)

    ss3<-lapply(split(rel.bias,list(taus,alphas)),ff.relb)
    ss3<-lapply(1:length(ss3),funk,ss3,matrix(gg1[,!my_grepl("strata",colnames(gg1))],ncol=sum(!my_grepl("strata",colnames(gg1)))))

    ss1<-lapply(1:length(ss1),function(i,x,y) {rr<-rbind(x[[i]],y[[i]]);colnames(rr)<-colnames(x[[i]]); rownames(rr)<-c(rownames(x[[i]]),"Relative Cond. Bias") ;rr   }, ss1,ss3)
    names(ss1)<-names(ss2)

    }

if (object$criterion=="alpha"&is.null(object$tau)) {
    ss<-lapply(split(vars.model,list(alphas)),function(x) table(unlist(x))[vars.num]/object$num.boot)
    ss.col<-lapply(split(vars.model.col,list(alphas)),function(x) table(unlist(x))[order(-table(unlist(x)))]/object$num.boot)

    gg<-lapply(ss,function(x) {z<-ifelse(vars.num%in%names(x)==F,0,x);names(z)<-vars.num;z} )

    gg1<-matrix(unlist(gg),byrow=T,ncol=length(vars.num))
    colnames(gg1)<-vars.num
    rownames(gg1)<-names(ss)

    ss1<-lapply(split(mm,list(alphas)),ff)
    ss2<-lapply(split(rmsd,list(alphas)),ff1)

    ss1<-lapply(1:length(ss1),function(i,x,y) {rr<-rbind(x[[i]],y[[i]]);colnames(rr)<-colnames(x[[i]]); rownames(rr)<-c(rownames(x[[i]]),"RMSD ratio") ;rr   }, ss1,ss2)
    names(ss1)<-names(ss2)


    ss3<-lapply(split(rel.bias,list(alphas)),ff.relb)
    ss3<-lapply(1:length(ss3),funk,ss3,matrix(gg1[,!my_grepl("strata",colnames(gg1))],ncol=sum(!my_grepl("strata",colnames(gg1)))))

    ss1<-lapply(1:length(ss1),function(i,x,y) {rr<-rbind(x[[i]],y[[i]]);colnames(rr)<-colnames(x[[i]]); rownames(rr)<-c(rownames(x[[i]]),"Relative Cond. Bias") ;rr   }, ss1,ss3)
    names(ss1)<-names(ss2)

    }

if (object$criterion!="alpha"&!is.null(object$tau)){
    ss<-lapply(split(vars.model,list(taus)),function(x) table(unlist(x))[vars.num]/object$num.boot)
    ss.col<-lapply(split(vars.model.col,list(taus)),function(x) table(unlist(x))[order(-table(unlist(x)))]/object$num.boot)

    gg<-lapply(ss,function(x) {z<-ifelse(vars.num%in%names(x)==F,0,x);names(z)<-vars.num;z} )

    gg1<-matrix(unlist(gg),byrow=T,ncol=length(vars.num))
    colnames(gg1)<-vars.num
    rownames(gg1)<-names(ss)
    ss1<-lapply(split(mm,list(taus)),ff)
    ss2<-lapply(split(rmsd,list(taus)),ff1)

    ss1<-lapply(1:length(ss1),function(i,x,y) {rr<-rbind(x[[i]],y[[i]]);colnames(rr)<-colnames(x[[i]]); rownames(rr)<-c(rownames(x[[i]]),"RMSD ratio") ;rr   }, ss1,ss2)
    names(ss1)<-names(ss2)

    ss3<-lapply(split(rel.bias,list(taus)),ff.relb)
    ss3<-lapply(1:length(ss3),funk,ss3,matrix(gg1[,!my_grepl("strata",colnames(gg1))],ncol=sum(!my_grepl("strata",colnames(gg1)))))

    ss1<-lapply(1:length(ss1),function(i,x,y) {rr<-rbind(x[[i]],y[[i]]);colnames(rr)<-colnames(x[[i]]); rownames(rr)<-c(rownames(x[[i]]),"Relative Cond. Bias") ;rr   }, ss1,ss3)
    names(ss1)<-names(ss2)

    }


if (object$criterion!="alpha"&is.null(object$tau)) {

    ss<-lapply(vars.model,function(x) table(unlist(x))[vars.num]/object$num.boot)
    ss.col<-lapply(vars.model.col,function(x) table(unlist(x))[order(-table(unlist(x)))]/object$num.boot)

    gg<-lapply(ss,function(x) {z<-ifelse(vars.num%in%names(x)==F,0,x);names(z)<-vars.num;z} )

    gg1<-matrix(unlist(gg),byrow=T,ncol=length(vars.num))
    rownames(gg1)<-names(ss)
    colnames(gg1)<-vars.num
    ss1<-lapply(mm,ff)
    ss1<-lapply(mm,ff)
    ss2<-lapply(rmsd,ff1)

    ss1<-lapply(1:length(ss1),function(i,x,y) {rr<-rbind(x[[i]],y[[i]]);colnames(rr)<-colnames(x[[i]]); rownames(rr)<-c(rownames(x[[i]]),"RMSD ratio") ;rr   }, ss1,ss2)
    names(ss1)<-names(ss2)

    ss3<-lapply( rel.bias,ff.relb)
    ss3<-lapply(1:length(ss3),funk,ss3,matrix(gg1[,!my_grepl("strata",colnames(gg1))],ncol=sum(!my_grepl("strata",colnames(gg1)))))

    ss1<-lapply(1:length(ss1),function(i,x,y) {rr<-rbind(x[[i]],y[[i]]);colnames(rr)<-colnames(x[[i]]); rownames(rr)<-c(rownames(x[[i]]),"Relative Cond. Bias") ;rr   }, ss1,ss3)
    names(ss1)<-names(ss2)

    }

list(var.rel.frequencies=gg1,model.rel.frequencies=ss.col,var.coefs=ss1)

}


#' Print Function
#'
#' Prints a summary table of a bootstrapped/resampled version of ABE.
#' The table displays the coefficient estimates and standard errors from the initial model (model with all covariates),
#' the relative inclusion frequencies of the covariates from the initial model,
#' resampled median and percentiles for the estimates of the regression coefficients for each variable from the initial model,
#' root mean squared difference ratio (RMSD) and relative bias conditional on selection (RBCS).
#' In cotrast to funciton \code{\link{print.abe.default}} which calculates these measures using either bootstrap or resampling (accordng to the latest guideliness), the function \code{\link{print.abe}} calculates the measures using only the option specified in the argument \code{type.boot} in a call to \code{\link{abe.resampling}}.
#'
#' @param x an object of class \code{"abe"}, an object returned by a call to \code{\link{abe.resampling}}
#' @param conf.level the confidence level, defaults to 0.95
#' @param alpha the alpha value for which the output is to be printed, defaults to \code{NULL}
#' @param tau the tau value for which the output is to be printed, defaults to \code{NULL}
#' @param ... additional arguments affecting the summary produced.
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @author Sladana Babic
#' @seealso \code{\link{abe.resampling}}, \code{\link{summary.abe}}, \code{\link{print.abe.default}}, \code{\link{plot.abe}}, \code{\link{pie.abe}}
#' @export
#' @examples
#' set.seed(1)
#' n=100
#' x1<-runif(n)
#' x2<-runif(n)
#' x3<-runif(n)
#' y<--5+5*x1+5*x2+ rnorm(n,sd=5)
#' dd<-data.frame(y=y,x1=x1,x2=x2,x3=x3)
#' fit<-lm(y~x1+x2+x3,x=TRUE,y=TRUE,data=dd)
#'
#' fit.boot<-abe.resampling(fit,data=dd,include="x1",active="x2",
#' tau=c(0.05,0.1),exp.beta=FALSE,exact=TRUE,
#' criterion="alpha",alpha=c(0.2,0.05),type.test="Chisq",
#' num.boot=50,type.boot="bootstrap")
#'
#' print(fit.boot,conf.level=0.95,alpha=0.2,tau=0.05)


print.abe<-function(x,conf.level=0.95,alpha=NULL,tau=NULL,...){
  object<-x
  if (conf.level<0|conf.level>1) stop("Confidence level out of range")
  if (!is.null(tau)) if(tau%in%object$misc$tau==FALSE) stop("This value of tau was not used in a call to abe.resampling.")
  if (object$misc$criterion=="alpha") if(!is.null(alpha)) if(alpha%in%object$misc$alpha==FALSE) stop("This value of alpha was not used in a call to abe.resampling.")


  if (object$misc$criterion=="alpha") if(is.null(alpha)) alpha=object$misc$alpha[1]

  if (is.null(tau)) tau=object$misc$tau[1]

  cat("Printing results of a call to abe.resampling for:\n  tau=",tau,"\n  criterion=\"",object$misc$criterion,"\"",sep="")
  if (object$misc$criterion=="alpha" ) cat("  alpha=",alpha,sep="")

  set<-paste("tau=",tau,sep="")
  sea<-paste("alpha=",alpha,sep="")

  if (object$misc$criterion=="alpha") rs<-paste(set,sea,sep=".") else rs<-set
ss<-summary(object,conf.level = conf.level)

vars.num<-c(object$all.vars,attributes(object$fit.or$terms)$term.labels[my_grepl("strata",attributes(object$fit.or$terms)$term.labels)])

idv<-which(vars.num%in%object$all.vars==TRUE)
idg<-which(vars.num%in%object$all.vars==FALSE)


if (length(idv)!=length(vars.num)) {
   pt<-paste(paste(names(ss$var.rel.frequencies[rs,])[idg],": ",ss$var.rel.frequencies[rs,idg],sep=""),collapse=",")
   cat("\n\n Inclusion relative frequencies of the offset/stratification variable(s):",pt)

}
 mat <- cbind(coef(object$fit.or),
                                sqrt(diag(vcov(object$fit.or))),
                                  c(ss$var.rel.frequencies[rs,idv]),
                                  t(ss$var.coefs[[rs]][,idv]))
rownames(mat)<-names(ss$var.rel.frequencies[rs,])[idv]

 colnames(mat)[1:3] <- c("Estimate init.", "Std. Error init.",
                        "Incl. Freq.")
colnames(mat)[4:6]<-c("Estimate, 50%", paste("Estimate ",(1-conf.level)/2*100,"%",sep=""),paste("Estimate ",100-(1-conf.level)/2*100,"%",sep="") )
colnames(mat)[7:8]<-c("RMSD ratio","RBCS")
cat("\n\n")

print(mat)
}

#' Plot Function
#'
#' Plot function for the resampled/bootstrapped version of ABE.
#'
#' @param x an object of class \code{"abe"}, an object returned by a call to \code{\link{abe.resampling}}
#' @param type.plot string which specifies the type of the plot. See details.
#' @param alpha values of alpha for which the plot is to be made (can be a vector of length >1)
#' @param tau values of tau for which the plot is to be made (can be a vector of length >1)
#' @param variable variables for which the plot is to be made (can be a vector of length >1)
#' @param ... Arguments to be passed to methods, such as graphical parameters (see \code{\link{barplot}}, \code{\link{hist}}).
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @author Sladana Babic
#' @details when using \code{type.plot="coefficients"} the function plots a histogram of the estimated regression coefficients for the specified variables, alpha(s) and tau(s) obtained from different re-sampled datasets.
#' When the variable is not included in the final model, its regression coefficient is set to zero.
#' When using \code{type.plot="variables"} the function plots a barplot of the relative inclusion frequencies of the specified variables, for the specified values of alpha and tau.
#' When using \code{type.plot="models"} the function plots a barplot of the relative frequencies of the final models for specified alpha(s) and tau(s).
#' @export
#' @seealso \code{\link{abe.resampling}}, \code{\link{summary.abe}}, \code{\link{pie.abe}}
#' @examples
#' set.seed(1)
#' n=100
#' x1<-runif(n)
#' x2<-runif(n)
#' x3<-runif(n)
#' y<--5+5*x1+5*x2+ rnorm(n,sd=5)
#' dd<-data.frame(y=y,x1=x1,x2=x2,x3=x3)
#' fit<-lm(y~x1+x2+x3,x=TRUE,y=TRUE,data=dd)
#'
#' fit.boot<-abe.resampling(fit,data=dd,include="x1",active="x2",
#' tau=c(0.05,0.1),exp.beta=FALSE,exact=TRUE,
#' criterion="alpha",alpha=c(0.2,0.05),type.test="Chisq",
#' num.boot=50,type.boot="bootstrap")
#'
#' plot(fit.boot,type.plot="coefficients",
#' alpha=0.2,tau=0.1,variable=c("x1","x3"),
#' col="light blue")
#'
#' plot(fit.boot,type.plot="variables",
#' alpha=0.2,tau=0.1,variable=c("x1","x2","x3"),
#' col="light blue",horiz=TRUE,las=1)
#'
#' par(mar=c(4,6,4,2))
#' plot(fit.boot,type.plot="models",
#' alpha=0.2,tau=0.1,col="light blue",horiz=TRUE,las=1)


plot.abe<-function(x,type.plot=c("coefficients","models","variables"),alpha=NULL,tau=NULL,variable=NULL,...){
object<-x

ggff<-function(x){ tbx<-table(x); {for (jj in which((tbx>1)==T)) {x[x==names(tbx[jj])]<-paste(x[x==names(tbx[jj])],1:sum(x==names(tbx[jj])),sep="")  }} ;x  }


object$all.vars<-ggff(object$all.vars)


  if (!is.null(variable)){
      for (i in 1:length(variable)) if(sum(my_grepl(variable[i],object$all.vars))==0) stop("At least one specified variable was not included in the initial model.")

       variable.l<-list()
      for (i in 1:length(variable)) {
        variable.l[[i]]<-object$all.vars[my_grep(variable[i],object$all.vars)]
      }
      variable<-unlist(variable.l)
     }



  if (object$criterion!="alpha") alphas<-NULL else {

    if (!is.null(object$tau)) alphas<-paste("alpha=",rep(object$alpha,each=length(object$tau)*object$num.boot),sep="")
    if (is.null(object$tau)) alphas<-paste("alpha=",rep(object$alpha,each=1*object$num.boot),sep="")
  }


  if (is.null(object$tau)) taus<-NULL else {

    if (is.null(object$alpha)) {taus<- paste("tau=",rep(rep( object$tau,each=object$num.boot  ),1),sep="")} else {

      taus<- paste("tau=",rep(rep( object$tau,each=object$num.boot  ),length(object$alpha)),sep="")
    }
  }


  if (!is.null(alpha)) if (sum(paste("alpha=",alpha,sep="")%in%unique(alphas))!=length(alpha)) stop("This value of alpha was not considered when using abe.resampling.")
  if (!is.null(tau)) if (sum(paste("tau=",tau,sep="")%in%unique(taus))!=length(tau)) stop("This value of tau was not considered when using abe.resampling.")




  if (!is.null(object$tau)&!is.null(object$alpha)) boot.iter<- rep(1:object$num.boot,length(object$alpha)*length(object$tau))

  if (is.null(object$tau)&!is.null(object$alpha)) boot.iter<- rep(1:object$num.boot,length(object$alpha))

  if (!is.null(object$tau)&is.null(object$alpha)) boot.iter<- rep(1:object$num.boot,length(object$tau))


vars.model<-lapply(object$models,function(x) if (is.numeric(x)) "empty model" else c(names(coef(x)),attributes(x$terms)$term.labels[my_grepl("strata",attributes(x$terms)$term.labels)]))
vars.model.cf<-lapply(object$models,function(x) if (is.numeric(x)) "empty model" else names(coef(x)))


ggff<-function(x){ tbx<-table(x); {for (jj in which((tbx>1)==T)) {x[x==names(tbx[jj])]<-paste(x[x==names(tbx[jj])],1:sum(x==names(tbx[jj])),sep="")  }} ;x  }


vars.model<-lapply(vars.model,ggff  )
vars.model.cf<-lapply(vars.model.cf,ggff   )

coefs.model<-lapply(object$models,function(x) if (is.numeric(x)) NULL else  coef(x))


object$all.vars<-ggff(object$all.vars)


coefs.model<-lapply(1:length(vars.model.cf),function(i,x,y,z) {zz<-list();zz[[i]]<-rep(0,length(y));zz[[i]][y%in%x[[i]]==T]<-z[[i]];names(zz[[i]])<-y;zz[[i]]},  vars.model.cf,object$all.vars,coefs.model  )





if (!is.null(alpha)|!is.null(tau)){

  if ( !is.null(alpha)&!is.null(tau)){
    tau.i<-paste("tau=",tau,sep="")
    alpha.i<-paste("alpha=",alpha,sep="")
    coefs.model<-coefs.model[alphas%in%alpha.i&taus%in%tau.i]
    alphas.ii<-alphas[alphas%in%alpha.i&taus%in%tau.i]
    taus<-taus[alphas%in%alpha.i&taus%in%tau.i]
    alphas<-alphas.ii
  }

  if ( !is.null(alpha)&is.null(tau)){
    alpha.i<-paste("alpha=",alpha,sep="")
    coefs.model<-coefs.model[alphas%in%alpha.i]
    alphas.ii<-alphas[alphas%in%alpha.i]
    taus<-taus[alphas%in%alpha.i]
    alphas<-alphas.ii
    }


  if ( is.null(alpha)&!is.null(tau)){
    tau.i<-paste("tau=",tau,sep="")
    coefs.model<-coefs.model[taus%in%tau.i]
    taus.ii<-taus[taus%in%tau.i]
    alphas<-alphas[taus%in%tau.i]
    taus<-taus.ii}


}

if (type.plot=="coefficients"){

if (object$criterion=="alpha"&!is.null(object$tau)) {
 ss<-lapply(split(coefs.model,list(taus,alphas)),function(x) {mm<-matrix(unlist(x),ncol=length(object$all.vars),nrow=object$num.boot,byrow=T);colnames(mm)<-object$all.vars;mm})
  if (!is.null(variable)) ss<-lapply(ss,function(x) {xi<-matrix(x[,colnames(x)%in%variable],ncol=sum(colnames(x)%in%variable),nrow=object$num.boot);colnames(xi)=object$all.vars[colnames(x)%in%variable] ;xi})
  par(mfcol=c(length(ss),ncol(ss[[1]])))

 for (i in 1:ncol(ss[[1]])){
   for (j in 1:length(ss)){
     hist(ss[[j]][,i],xlab=colnames(ss[[1]])[i],main=names(ss)[j],...)
   }
 }

}

if (object$criterion=="alpha"&is.null(object$tau)) {
  ss<-lapply(split(coefs.model,list(alphas)),function(x) {mm<-matrix(unlist(x),ncol=length(object$all.vars),nrow=object$num.boot,byrow=T);colnames(mm)<-object$all.vars;mm})
  if (!is.null(variable)) ss<-lapply(ss,function(x) {xi<-matrix(x[,colnames(x)%in%variable],ncol=sum(colnames(x)%in%variable),nrow=object$num.boot);colnames(xi)=object$all.vars[colnames(x)%in%variable] ;xi})

  par(mfcol=c(length(ss),ncol(ss[[1]])))

  for (i in 1:ncol(ss[[1]])){
    for (j in 1:length(ss)){
      hist(ss[[j]][,i],xlab=colnames(ss[[1]])[i],main=names(ss)[j],...)
    }
  }

}

if (object$criterion!="alpha"&!is.null(object$tau)){
  ss<-lapply(split(coefs.model,list(taus)),function(x) {mm<-matrix(unlist(x),ncol=length(object$all.vars),nrow=object$num.boot,byrow=T);colnames(mm)<-object$all.vars;mm})
  if (!is.null(variable)) ss<-lapply(ss,function(x) {xi<-matrix(x[,colnames(x)%in%variable],ncol=sum(colnames(x)%in%variable),nrow=object$num.boot);colnames(xi)=object$all.vars[colnames(x)%in%variable] ;xi})

  par(mfcol=c(length(ss),ncol(ss[[1]])))

  for (i in 1:ncol(ss[[1]])){
    for (j in 1:length(ss)){
      hist(ss[[j]][,i],xlab=colnames(ss[[1]])[i],main=names(ss)[j],...)
    }
  }

}


if (object$criterion!="alpha"&is.null(object$tau)) {

  ss<-lapply(coefs.model,function(x) {mm<-matrix(unlist(x),ncol=length(object$all.vars),nrow=object$num.boot,byrow=T);colnames(mm)<-object$all.vars;mm})
  if (!is.null(variable)) ss<-lapply(ss,function(x) {xi<-matrix(x[,colnames(x)%in%variable],ncol=sum(colnames(x)%in%variable),nrow=object$num.boot);colnames(xi)=object$all.vars[colnames(x)%in%variable] ;xi})

  par(mfcol=c(length(ss),ncol(ss[[1]])))

  for (i in 1:ncol(ss[[1]])){
      hist(ss[[1]][,i],xlab=colnames(ss[[1]])[i],main=names(ss)[1],...)
    }
  }

}


if (type.plot=="variables"){
  sum.obj<-summary(object)$var.rel.frequencies

  if (!is.null(variable)){
    sum.obji<-matrix(sum.obj[,colnames(sum.obj)%in%variable  ],ncol=sum(colnames(sum.obj)%in%variable),nrow=nrow(sum.obj))
    colnames(sum.obji)<-colnames(sum.obj)[colnames(sum.obj)%in%variable]
    rownames(sum.obji)<-rownames(sum.obj)
      sum.obj<-sum.obji
  }


  if (is.null(alpha)&is.null(tau)) cnm<-rownames(sum.obj) else {
    if (!is.null(alphas)&!is.null(taus)) {

      if (!is.null(tau)&!is.null(alpha)) cnm<-paste("tau=",tau,".alpha=",alpha,sep="")
      if (is.null(tau)&!is.null(alpha)) cnm<-paste("tau=",object$tau,".alpha=",alpha,sep="")
      if (!is.null(tau)&is.null(alpha)) cnm<-paste("tau=",tau,".alpha=",object$alpha,sep="")

      } else {
    if (is.null(alphas)) cnm<-paste("tau=",tau,sep="")
    if (is.null(tau)) cnm<-paste("alpha=",alpha,sep="")
    }
  }

    sum.obji<-matrix(sum.obj[ rownames(sum.obj)%in%cnm  ,  ],nrow=sum(rownames(sum.obj)%in%cnm),ncol=ncol(sum.obj))
    colnames(sum.obji)<-colnames(sum.obj)
    rownames(sum.obji)<-rownames(sum.obj)[rownames(sum.obj)%in%cnm]
    sum.obj<-sum.obji





    par(mfcol=c(nrow(sum.obj),1)   )
    for (i in 1:nrow(sum.obj)) barplot(sum.obj[i,],main=rownames(sum.obj)[i],names=colnames(sum.obj),... )


}


if (type.plot=="models"){
  if (!is.null(variable)) warning("Ploting relative frequencies of the final models but variable is not null. Ignoring the argument variable.")
  sum.obj<-summary(object)$model.rel.frequencies

  if (is.null(alpha)&is.null(tau)) cnm<-names(sum.obj) else {
    if (!is.null(alphas)&!is.null(taus)) {

      if (!is.null(tau)&!is.null(alpha)) cnm<-paste("tau=",tau,".alpha=",alpha,sep="")
      if (is.null(tau)&!is.null(alpha)) cnm<-paste("tau=",object$tau,".alpha=",alpha,sep="")
      if (!is.null(tau)&is.null(alpha)) cnm<-paste("tau=",tau,".alpha=",object$alpha,sep="")

    } else {
      if (is.null(alphas)) cnm<-paste("tau=",tau,sep="")
      if (is.null(tau)) cnm<-paste("alpha=",alpha,sep="")
    }
  }

  sum.obj<-sum.obj[names(sum.obj)%in%cnm]

  par(mfcol=c(length(sum.obj),1)   )
  for (i in 1:length(sum.obj)) barplot(sum.obj[[i]],main=names(sum.obj)[i],names=names(sum.obj[[i]]),... )



}

}




#' Pie Function
#'
#' Pie function for the resampled/bootstrapped version of ABE. Plots a pie chart of the model frequencies for specified values of \code{alpha} and \code{tau}.
#'
#' @param x an object of class \code{"abe"}, an object returned by a call to \code{\link{abe.resampling}}
#' @param alpha values of alpha for which the plot is to be made (can be a vector of length >1)
#' @param tau values of tau for which the plot is to be made (can be a vector of length >1)
#' @param labels plot labels, defaults to NA, i.e. no labels are ploted
#' @param ... Arguments to be passed to methods, such as graphical parameters (see \code{\link{pie}}, \code{\link{barplot}}, \code{\link{hist}}).
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @author Sladana Babic
#' @export
#' @seealso \code{\link{abe.resampling}}, \code{\link{summary.abe}}, \code{\link{plot.abe}}
#' @examples
#' set.seed(1)
#' n=100
#' x1<-runif(n)
#' x2<-runif(n)
#' x3<-runif(n)
#' y<--5+5*x1+5*x2+ rnorm(n,sd=5)
#' dd<-data.frame(y=y,x1=x1,x2=x2,x3=x3)
#' fit<-lm(y~x1+x2+x3,x=TRUE,y=TRUE,data=dd)
#'
#' fit.boot<-abe.resampling(fit,data=dd,include="x1",active="x2",
#' tau=c(0.05,0.1),exp.beta=FALSE,exact=TRUE,
#' criterion="alpha",alpha=c(0.2,0.05),type.test="Chisq",
#' num.boot=50,type.boot="bootstrap")
#'
#' pie.abe(fit.boot, alpha=0.2,tau=0.1)



pie.abe<-function(x,alpha=NULL,tau=NULL,labels=NA,...){
  object<-x

  ggff<-function(x){ tbx<-table(x); {for (jj in which((tbx>1)==T)) {x[x==names(tbx[jj])]<-paste(x[x==names(tbx[jj])],1:sum(x==names(tbx[jj])),sep="")  }} ;x  }


  object$all.vars<-ggff(object$all.vars)



  if (object$criterion!="alpha") alphas<-NULL else {

    if (!is.null(object$tau)) alphas<-paste("alpha=",rep(object$alpha,each=length(object$tau)*object$num.boot),sep="")
    if (is.null(object$tau)) alphas<-paste("alpha=",rep(object$alpha,each=1*object$num.boot),sep="")
  }


  if (is.null(object$tau)) taus<-NULL else {

    if (is.null(object$alpha)) {taus<- paste("tau=",rep(rep( object$tau,each=object$num.boot  ),1),sep="")} else {

      taus<- paste("tau=",rep(rep( object$tau,each=object$num.boot  ),length(object$alpha)),sep="")
    }
  }


  if (!is.null(alpha)) if (sum(paste("alpha=",alpha,sep="")%in%unique(alphas))!=length(alpha)) stop("This value of alpha was not considered when using abe.resampling.")
  if (!is.null(tau)) if (sum(paste("tau=",tau,sep="")%in%unique(taus))!=length(tau)) stop("This value of tau was not considered when using abe.resampling.")




  if (!is.null(object$tau)&!is.null(object$alpha)) boot.iter<- rep(1:object$num.boot,length(object$alpha)*length(object$tau))

  if (is.null(object$tau)&!is.null(object$alpha)) boot.iter<- rep(1:object$num.boot,length(object$alpha))

  if (!is.null(object$tau)&is.null(object$alpha)) boot.iter<- rep(1:object$num.boot,length(object$tau))


  vars.model<-lapply(object$models,function(x) if (is.numeric(x)) "empty model" else c(names(coef(x)),attributes(x$terms)$term.labels[my_grepl("strata",attributes(x$terms)$term.labels)]))
  vars.model.cf<-lapply(object$models,function(x) if (is.numeric(x)) "empty model" else names(coef(x)))


  ggff<-function(x){ tbx<-table(x); {for (jj in which((tbx>1)==T)) {x[x==names(tbx[jj])]<-paste(x[x==names(tbx[jj])],1:sum(x==names(tbx[jj])),sep="")  }} ;x  }


  vars.model<-lapply(vars.model,ggff  )
  vars.model.cf<-lapply(vars.model.cf,ggff   )

  coefs.model<-lapply(object$models,function(x) if (is.numeric(x)) NULL else  coef(x))


  object$all.vars<-ggff(object$all.vars)


  coefs.model<-lapply(1:length(vars.model.cf),function(i,x,y,z) {zz<-list();zz[[i]]<-rep(0,length(y));zz[[i]][y%in%x[[i]]==T]<-z[[i]];names(zz[[i]])<-y;zz[[i]]},  vars.model.cf,object$all.vars,coefs.model  )





  if (!is.null(alpha)|!is.null(tau)){

    if ( !is.null(alpha)&!is.null(tau)){
      tau.i<-paste("tau=",tau,sep="")
      alpha.i<-paste("alpha=",alpha,sep="")
      coefs.model<-coefs.model[alphas%in%alpha.i&taus%in%tau.i]
      alphas.ii<-alphas[alphas%in%alpha.i&taus%in%tau.i]
      taus<-taus[alphas%in%alpha.i&taus%in%tau.i]
      alphas<-alphas.ii
    }

    if ( !is.null(alpha)&is.null(tau)){
      alpha.i<-paste("alpha=",alpha,sep="")
      coefs.model<-coefs.model[alphas%in%alpha.i]
      alphas.ii<-alphas[alphas%in%alpha.i]
      taus<-taus[alphas%in%alpha.i]
      alphas<-alphas.ii
    }


    if ( is.null(alpha)&!is.null(tau)){
      tau.i<-paste("tau=",tau,sep="")
      coefs.model<-coefs.model[taus%in%tau.i]
      taus.ii<-taus[taus%in%tau.i]
      alphas<-alphas[taus%in%tau.i]
      taus<-taus.ii}


  }




    sum.obj<-summary(object)$model.rel.frequencies

    if (is.null(alpha)&is.null(tau)) cnm<-names(sum.obj) else {
      if (!is.null(alphas)&!is.null(taus)) {

        if (!is.null(tau)&!is.null(alpha)) cnm<-paste("tau=",tau,".alpha=",alpha,sep="")
        if (is.null(tau)&!is.null(alpha)) cnm<-paste("tau=",object$tau,".alpha=",alpha,sep="")
        if (!is.null(tau)&is.null(alpha)) cnm<-paste("tau=",tau,".alpha=",object$alpha,sep="")

      } else {
        if (is.null(alphas)) cnm<-paste("tau=",tau,sep="")
        if (is.null(tau)) cnm<-paste("alpha=",alpha,sep="")
      }
    }

    sum.obj<-sum.obj[names(sum.obj)%in%cnm]

    par(mfcol=c(length(sum.obj),1)   )
    for (i in 1:length(sum.obj)) pie(sum.obj[[i]],main=names(sum.obj)[i],labels=labels )





}





#' ABE for models which include only numeric covariates
#' @keywords internal
#' @examples
#' \dontshow{
#' set.seed(1)
#' n=100
#' x1<-runif(n)
#' x2<-runif(n)
#' x3<-runif(n)
#' y<--5+5*x1+5*x2+ rnorm(n,sd=5)
#' dd<-data.frame(y,x1,x2,x3)
#' fit<-lm(y~x1+x2+x3,x=TRUE,y=TRUE,data=dd)
#'
#' abe.fit<-abe.num(fit,data=dd,include="x1",active="x2",
#' tau=0.05,exp.beta=FALSE,exact=TRUE,criterion="alpha",alpha=0.2,
#' type.test="Chisq",verbose=FALSE)
#' summary(abe.fit)
#' }

abe.num<-function(fit,data,include=NULL,active=NULL,tau=0.05,exp.beta=TRUE,exact=FALSE,criterion="alpha",alpha=0.2,type.test="Chisq",verbose=T){


if (criterion[1]=="alpha") k<-qchisq(1-alpha,df=1)
if (criterion[1]=="AIC") k<-2

n<-nrow(fit$x)

if (criterion[1]=="BIC") k<-log(n)




 if (colnames(model.matrix(fit))[1]=="(Intercept)") xm<-as.matrix(fit$x)[,-1] else xm<-as.matrix(fit$x)


vcvx<-var(xm)
 rownames(vcvx)<-colnames(vcvx)<-attributes(fit$terms)$term.labels[!my_grepl("strata",attributes(fit$terms)$term.labels)]


cnms<-attributes(fit$terms)$term.labels


if (!is.null(include)) {
  include.l<-list()
  for (i in 1:length(include)) {
    if (sum(my_grepl(include[i],cnms))==0) stop("at least one include variable is not in the model")
    include.l[[i]]<-cnms[my_grep(include[i],cnms)]
    }
  include<-unlist(include.l)
}
if (!is.null(active))  {
  active.l<-list()
  for (i in 1:length(active)) {
    if (sum(my_grepl(active[i],cnms))==0) stop("at least one active variable is not in the model")
    active.l[[i]]<-cnms[my_grep(active[i],cnms)]
    }
  active<-unlist(active.l)
}


if (sum(include%in%active)!=0) stop("at least one include variable is also specified as active")


if ( sum(cnms%in%include)==length(cnms) ) stop("all variables are specified as pasive, cannot perform variable selection")
if ( sum(cnms%in%active)==length(cnms) ) {
		include=active=NULL
		warning("all variables are specified as active, treating all variables as active or pasive")
		}



if (is.null(include)&is.null(active)) {varnfix<- cnms;varpas<-cnms}
if (is.null(active)&!is.null(include)) {varnfix<- cnms[!cnms%in%include]; varpas<-cnms}
if (!is.null(active)&is.null(include))  {varnfix<-cnms;varpas<-cnms[!cnms%in%active]  }
if (!is.null(active)&!is.null(include)) {varnfix<-cnms[!cnms%in%include];varpas<-cnms[!cnms%in%active]}

varpas<-varpas[!my_grepl("strata",varpas)]

stop=F

while(stop==F){

vcvm<-vcov(fit)

if (colnames(model.matrix(fit))[1]=="(Intercept)") vcvm<-vcvm[-1,-1]
if (is.null(dim(vcvm))) {
		vcvm<-matrix(vcvm,ncol=1,nrow=1)
		if (colnames(model.matrix(fit))[1]=="(Intercept)") colnames(vcvm)<-rownames(vcvm)<-colnames(vcov(fit))[-1] else colnames(vcvm)<-rownames(vcvm)<-colnames(vcov(fit))
		}


if (verbose==TRUE) {
	cat("\n\nModel under investigation:\n")
	print(fit$call)
	}


 if (criterion!="alpha") bl<-drop1(fit,scope=as.formula(paste("~",paste(varnfix,collapse=" + ") )),k=k) else bl<-drop1(fit,scope=as.formula(paste("~",paste(varnfix,collapse=" + ") )),test=type.test)

if (verbose==TRUE)  if (criterion!="alpha")  cat("Criterion for non-passive variables: "  ,paste( paste(varnfix,round(bl$AIC[-1],4),sep=" : "),collapse=" , "),"\n"  ,sep=""  ) else cat("Criterion for non-passive variables: "  ,paste( paste(varnfix,round(bl[-1,pmatch("Pr",names(bl))],4),sep=" : "),collapse=" , "),"\n"  ,sep=""  )


if (criterion!="alpha") black.list.i<-varnfix[which(bl$AIC[-1]<bl$AIC[1])]  else   black.list.i<-varnfix[which(bl[-1,pmatch("Pr",names(bl))]>alpha)]


if (length(black.list.i)!=0){


  if (criterion=="alpha") {
     black.list.i<-varnfix[order(-bl[-1,pmatch("Pr",names(bl))])][1:length(black.list.i)]
    criterion.i<-bl[-1,pmatch("Pr",names(bl))][order(-bl[-1,pmatch("Pr",names(bl))])][1:length(black.list.i)]
  } else {
    black.list.i<-varnfix[order(bl$AIC[-1])][1:length(black.list.i)]
    criterion.i<-bl$AIC[-1][order(bl$AIC[-1])][1:length(black.list.i)]-bl$AIC[1]
  }




	if (verbose==TRUE) cat( "  ",paste("black list: ", paste(paste(black.list.i,round(criterion.i,4),sep=" : "),collapse =", ")),"\n" )

		flag=T
		i=1
		while (flag==T&i<=length(black.list.i)){



			if (exact==T){

				xf<-as.formula(paste("~.-",black.list.i[i]  ))

		 		fit.i<-update(fit,xf,evaluate=FALSE)
				fit.i<-eval.parent(fit.i)

				if (colnames(model.matrix(fit))[1]=="(Intercept)") {

					change.in.estimate<-abs(fit$coef[which(attributes(fit$terms)$term.labels[!my_grepl("strata",attributes(fit$terms)$term.labels)]%in%varpas[!varpas%in%black.list.i[i]])+1]-fit.i$coef[which(!attributes(fit.i$terms)$term.labels[!my_grepl("strata",attributes(fit.i$terms)$term.labels)]%in%active[!active%in%black.list.i[i]])+1])
					} else {
					change.in.estimate<-abs(fit$coef[which(attributes(fit$terms)$term.labels[!my_grepl("strata",attributes(fit$terms)$term.labels)]%in%varpas[!varpas%in%black.list.i[i]])]-fit.i$coef[which(!attributes(fit.i$terms)$term.labels[!my_grepl("strata",attributes(fit.i$terms)$term.labels)]%in%active[!active%in%black.list.i[i]])])
					}


				} else {
				if (colnames(model.matrix(fit))[1]=="(Intercept)") {
					change.in.estimate<- abs(fit$coef[which(attributes(fit$terms)$term.labels[!my_grepl("strata",attributes(fit$terms)$term.labels)]==black.list.i[i])+1]*vcvm[rownames(vcvm)==black.list.i[i],colnames(vcvm)%in%varpas[!varpas%in%black.list.i[i]]] /vcvm[rownames(vcvm)==black.list.i[i],colnames(vcvm)==black.list.i[i]])
					names(change.in.estimate)<-colnames(vcvm)[colnames(vcvm)%in%varpas[!varpas%in%black.list.i[i]]]
					} else {
					change.in.estimate<- abs(fit$coef[which(attributes(fit$terms)$term.labels[!my_grepl("strata",attributes(fit$terms)$term.labels)]==black.list.i[i])]*vcvm[rownames(vcvm)==black.list.i[i],colnames(vcvm)%in%varpas[!varpas%in%black.list.i[i]]] /vcvm[rownames(vcvm)==black.list.i[i],colnames(vcvm)==black.list.i[i]])
					names(change.in.estimate)<-colnames(vcvm)[colnames(vcvm)%in%varpas[!varpas%in%black.list.i[i]]]
					}
				}


			if (exp.beta==TRUE) if (length(varpas[!varpas%in%black.list.i[i]])>1) ch.in.est<-exp( change.in.estimate*sqrt(diag(vcvx[rownames(vcvx)%in%varpas[!varpas%in%black.list.i[i]],colnames(vcvx)%in%varpas[!varpas%in%black.list.i[i]]]))  ) else ch.in.est<-exp( change.in.estimate*sqrt((vcvx[rownames(vcvx)%in%varpas[!varpas%in%black.list.i[i]],colnames(vcvx)%in%varpas[!varpas%in%black.list.i[i]]])) )
			if (exp.beta==FALSE) if (length(varpas[!varpas%in%black.list.i[i]])>1) ch.in.est<-change.in.estimate*sqrt(diag(vcvx[rownames(vcvx)%in%varpas[!varpas%in%black.list.i[i]],colnames(vcvx)%in%varpas[!varpas%in%black.list.i[i]]]))/sd(fit$y) else ch.in.est<-change.in.estimate*sqrt((vcvx[rownames(vcvx)%in%varpas[!varpas%in%black.list.i[i]],colnames(vcvx)%in%varpas[!varpas%in%black.list.i[i]]]))/sd(fit$y)

			if (verbose==TRUE) cat("          ",paste("Investigating change in b or exp(b) due to omitting variable ",black.list.i[i]," ; ",paste(paste(names(ch.in.est),round(ch.in.est,4), sep=" : ") ,collapse=", ") )  ,"\n")

			if (exp.beta==TRUE) if( sum(ch.in.est >= 1+tau)==0) {
						flag=F
						varnfix<-varnfix[-which(varnfix==black.list.i[i])]
						if(!black.list.i[i]%in%include) if (sum(varpas==black.list.i[i])==0) varpas<-varpas else varpas<-varpas[-which(varpas==black.list.i[i])]
						}

			if (exp.beta==FALSE) if (sum(  ch.in.est  >=tau  )==0)  {
						flag=F
						varnfix<-varnfix[-which(varnfix==black.list.i[i])]
						if(!black.list.i[i]%in%include) if (sum(varpas==black.list.i[i])==0) varpas<-varpas else  varpas<-varpas[-which(varpas==black.list.i[i])]
						}

   			if (flag==F) {
   			  fit<-update(fit,as.formula(paste("~.-",black.list.i[i])),evaluate=FALSE)
   			  fit<-eval.parent(fit)
   			   }
			if (verbose==TRUE)  if (flag==T&i!=length(black.list.i))  cat("  ","updated black list:" , paste(paste(black.list.i[-(1:i)],round(criterion.i[-(1:i)],4),sep=" : "),collapse =", "),"\n" )

			i=i+1
			}


	} else {
	if (verbose==T) cat("black list: empty","\n")
	flag=T
	}

if (length(black.list.i)==0) stop=T
if (length(varnfix)==0) stop=T
if (length(include)>0) if (length(varpas)==length(include)& ( length(varnfix[!varnfix%in%varpas])==0 )) stop=T
if (flag==T) stop=T

}

if (verbose==T) {
	cat("\n\nFinal model:\n")
	print(fit$call)
	cat("\n\n")
	}

fit

}







#' ABE for model which include only numeric covariates, bootstrap version
#' @keywords internal
#' @examples
#' \dontshow{
#' set.seed(1)
#' n=100
#' x1<-runif(n)
#' x2<-runif(n)
#' x3<-runif(n)
#' y<--5+5*x1+5*x2+ rnorm(n,sd=5)
#' dd<-data.frame(y,x1,x2,x3)
#' fit<-lm(y~x1+x2+x3,x=TRUE,y=TRUE,data=dd)
#'
#' abe.fit<-abe.num.boot(fit,data=dd,include="x1",active="x2",
#' tau=0.05,exp.beta=FALSE,exact=TRUE,criterion="alpha",alpha=0.2,
#' type.test="Chisq",k=2)
#'
#' summary(abe.fit)
#' }

abe.num.boot<-function(fit,data,include=NULL,active=NULL,tau=0.05,exp.beta=TRUE,exact=FALSE,criterion="alpha",alpha=0.2,type.test="Chisq",k){




  if (colnames(model.matrix(fit))[1]=="(Intercept)") xm<-as.matrix(fit$x)[,-1] else xm<-as.matrix(fit$x)


  vcvx<-var(xm)
    rownames(vcvx)<-colnames(vcvx)<-attributes(fit$terms)$term.labels[!my_grepl("strata",attributes(fit$terms)$term.labels)]


  cnms<-attributes(fit$terms)$term.labels

  if (is.null(include)&is.null(active)) {varnfix<- cnms;varpas<-cnms}
  if (is.null(active)&!is.null(include)) {varnfix<- cnms[!cnms%in%include]; varpas<-cnms}
  if (!is.null(active)&is.null(include))  {varnfix<-cnms;varpas<-cnms[!cnms%in%active]  }
  if (!is.null(active)&!is.null(include)) {varnfix<-cnms[!cnms%in%include];varpas<-cnms[!cnms%in%active]}
  varpas<-varpas[!my_grepl("strata",varpas)]


  stop=F

  while(stop==F){

    vcvm<-vcov(fit)

    if (colnames(model.matrix(fit))[1]=="(Intercept)") vcvm<-vcvm[-1,-1]
    if (is.null(dim(vcvm))) {
      vcvm<-matrix(vcvm,ncol=1,nrow=1)
      if (colnames(model.matrix(fit))[1]=="(Intercept)") colnames(vcvm)<-rownames(vcvm)<-colnames(vcov(fit))[-1] else colnames(vcvm)<-rownames(vcvm)<-colnames(vcov(fit))
    }




      if (criterion!="alpha") bl<-drop1(fit,scope=as.formula(paste("~",paste(varnfix,collapse=" + ") )),k=k) else bl<-drop1(fit,scope=as.formula(paste("~",paste(varnfix,collapse=" + ") )),test=type.test)


      if (criterion!="alpha") black.list.i<-varnfix[which(bl$AIC[-1]<bl$AIC[1])]  else   black.list.i<-varnfix[which(bl[-1,pmatch("Pr",names(bl))]>alpha)]


    if (length(black.list.i)!=0){


      if (criterion=="alpha") {
         black.list.i<-varnfix[order(-bl[-1,pmatch("Pr",names(bl))])][1:length(black.list.i)]
        criterion.i<-bl[-1,pmatch("Pr",names(bl))][order(-bl[-1,pmatch("Pr",names(bl))])][1:length(black.list.i)]
      } else {
        black.list.i<-varnfix[order(bl$AIC[-1])][1:length(black.list.i)]
        criterion.i<-bl$AIC[-1][order(bl$AIC[-1])][1:length(black.list.i)]-bl$AIC[1]
      }




        flag=T
      i=1
      while (flag==T&i<=length(black.list.i)){



        if (exact==T){

          xf<-as.formula(paste("~.-",black.list.i[i]  ))

           fit.i<-update(fit,xf,evaluate=FALSE)
          fit.i<-eval.parent(fit.i)

          if (colnames(model.matrix(fit))[1]=="(Intercept)") {

            change.in.estimate<-abs(fit$coef[which(attributes(fit$terms)$term.labels[!my_grepl("strata",attributes(fit$terms)$term.labels)]%in%varpas[!varpas%in%black.list.i[i]])+1]-fit.i$coef[which(!attributes(fit.i$terms)$term.labels[!my_grepl("strata",attributes(fit.i$terms)$term.labels)]%in%active[!active%in%black.list.i[i]])+1])
          } else {
            change.in.estimate<-abs(fit$coef[which(attributes(fit$terms)$term.labels[!my_grepl("strata",attributes(fit$terms)$term.labels)]%in%varpas[!varpas%in%black.list.i[i]])]-fit.i$coef[which(!attributes(fit.i$terms)$term.labels[!my_grepl("strata",attributes(fit.i$terms)$term.labels)]%in%active[!active%in%black.list.i[i]])])
          }


        } else {
          if (colnames(model.matrix(fit))[1]=="(Intercept)") {
            change.in.estimate<- abs(fit$coef[which(attributes(fit$terms)$term.labels[!my_grepl("strata",attributes(fit$terms)$term.labels)]==black.list.i[i])+1]*vcvm[rownames(vcvm)==black.list.i[i],colnames(vcvm)%in%varpas[!varpas%in%black.list.i[i]]] /vcvm[rownames(vcvm)==black.list.i[i],colnames(vcvm)==black.list.i[i]])
            names(change.in.estimate)<-colnames(vcvm)[colnames(vcvm)%in%varpas[!varpas%in%black.list.i[i]]]
          } else {
            change.in.estimate<- abs(fit$coef[which(attributes(fit$terms)$term.labels[!my_grepl("strata",attributes(fit$terms)$term.labels)]==black.list.i[i])]*vcvm[rownames(vcvm)==black.list.i[i],colnames(vcvm)%in%varpas[!varpas%in%black.list.i[i]]] /vcvm[rownames(vcvm)==black.list.i[i],colnames(vcvm)==black.list.i[i]])
            names(change.in.estimate)<-colnames(vcvm)[colnames(vcvm)%in%varpas[!varpas%in%black.list.i[i]]]
          }
        }


        if (exp.beta==TRUE) if (length(varpas[!varpas%in%black.list.i[i]])>1) ch.in.est<-exp( change.in.estimate*sqrt(diag(vcvx[rownames(vcvx)%in%varpas[!varpas%in%black.list.i[i]],colnames(vcvx)%in%varpas[!varpas%in%black.list.i[i]]]))  ) else ch.in.est<-exp( change.in.estimate*sqrt((vcvx[rownames(vcvx)%in%varpas[!varpas%in%black.list.i[i]],colnames(vcvx)%in%varpas[!varpas%in%black.list.i[i]]])) )
        if (exp.beta==FALSE) if (length(varpas[!varpas%in%black.list.i[i]])>1) ch.in.est<-change.in.estimate*sqrt(diag(vcvx[rownames(vcvx)%in%varpas[!varpas%in%black.list.i[i]],colnames(vcvx)%in%varpas[!varpas%in%black.list.i[i]]]))/sd(fit$y) else ch.in.est<-change.in.estimate*sqrt((vcvx[rownames(vcvx)%in%varpas[!varpas%in%black.list.i[i]],colnames(vcvx)%in%varpas[!varpas%in%black.list.i[i]]]))/sd(fit$y)


        if (exp.beta==TRUE) if( sum(ch.in.est >= 1+tau)==0) {
          flag=F
          varnfix<-varnfix[-which(varnfix==black.list.i[i])]
          if(!black.list.i[i]%in%include) if (sum(varpas==black.list.i[i])==0) varpas<-varpas else varpas<-varpas[-which(varpas==black.list.i[i])]
        }

        if (exp.beta==FALSE) if (sum(  ch.in.est  >=tau  )==0)  {
          flag=F
          varnfix<-varnfix[-which(varnfix==black.list.i[i])]
          if(!black.list.i[i]%in%include) if (sum(varpas==black.list.i[i])==0) varpas<-varpas else  varpas<-varpas[-which(varpas==black.list.i[i])]
        }

        if (flag==F) {
          fit<-update(fit,as.formula(paste("~.-",black.list.i[i])),evaluate=FALSE)
          fit<-eval.parent(fit)

          }

        i=i+1
      }


    } else {
        flag=T
    }

    if (length(black.list.i)==0) stop=T
    if (length(varnfix)==0) stop=T
    if (length(include)>0) if (length(varpas)==length(include)& ( length(varnfix[!varnfix%in%varpas])==0 )) stop=T
    if (flag==T) stop=T

  }


  fit

}



#' ABE for model which includes categorical covariates, factor option
#' @keywords internal
#' @examples
#' \dontshow{
#' set.seed(1)
#' n=100
#' x1<-runif(n)
#' x2<-runif(n)
#' x3<-runif(n)
#' y<--5+5*x1+5*x2+ rnorm(n,sd=5)
#' dd<-data.frame(y,x1,x2,x3)
#' fit<-lm(y~x1+x2+x3,x=TRUE,y=TRUE,data=dd)
#'
#' abe.fit<-abe.fact1(fit,data=dd,include="x1",active="x2",
#' tau=0.05,exp.beta=FALSE,exact=TRUE,criterion="alpha",alpha=0.2,
#' type.test="Chisq",verbose=FALSE)
#' summary(abe.fit)
#' }



abe.fact1<-function(fit,data,include=NULL,active=NULL,tau=0.05,exp.beta=TRUE,exact=FALSE,criterion="alpha",alpha=0.2,type.test="Chisq",verbose=T){



if (exact==F) {warning("there are factors in the model, using approximate change-in-estimate with this type.factor is inappropriate; using exact change in estimate instead"); exact=T}



 if (criterion[1]=="AIC") k<-2

n<-nrow(fit$x)

if (criterion[1]=="BIC") k<-log(n)



 if (colnames(model.matrix(fit))[1]=="(Intercept)") xm<-as.matrix(fit$x)[,-1] else xm<-as.matrix(fit$x)


var.mod<-attributes(fit$terms)$term.labels

if (length(var.mod)==0) fit$assign<-0 else {

fit$assign<-list()
fit$assign[1]<-1
name.cf<-names(coef(fit))

  for (i in 1:length(var.mod)){
      fit$assign[[i+1]]<-length(grep(paste(unlist(strsplit(paste(unlist(strsplit(var.mod[i],split="\\(")),collapse="\\("),split="\\)")),collapse="\\)"),name.cf))
  }

fit$assign<-unlist(fit$assign)

fit$assign<-rep(0:length(var.mod),fit$assign)
}

cnms<-rep(attributes(fit$terms)$term.labels[!my_grepl("strata",var.mod)],table(fit$assign[-1]))


 vcvx<-var(xm)
rownames(vcvx)<-colnames(vcvx)<-cnms

if (sum(my_grepl("strata",var.mod))!=0) cnms<-c(cnms,var.mod[my_grepl("strata",var.mod)])


if (!is.null(include)) {
  include.l<-list()
  for (i in 1:length(include)) {
    if (sum(my_grepl(include[i],cnms))==0) stop("at least one include variable is not in the model")
    include.l[[i]]<-cnms[my_grep(include[i],cnms)]
  }
  include<-unlist(include.l)
}
if (!is.null(active))  {
  active.l<-list()
  for (i in 1:length(active)) {
    if (sum(my_grepl(active[i],cnms))==0) stop("at least one active variable is not in the model")
    active.l[[i]]<-cnms[my_grep(active[i],cnms)]
  }
  active<-unlist(active.l)
}



if (sum(include%in%active)!=0) stop("at least one include variable is also specified as active")



if ( sum(cnms%in%include)==length(cnms) ) stop("all variables are specified as pasive, cannot perform variable selection")
if ( sum(cnms%in%active)==length(cnms) ) {
		include=active=NULL
		warning("all variables are specified as active, treating all variables as active or pasive")
		}


if (is.null(include)&is.null(active)) {varnfix<- cnms;varpas<-cnms}
if (is.null(active)&!is.null(include)) {varnfix<- cnms[!cnms%in%include]; varpas<-cnms}
if (!is.null(active)&is.null(include))  {varnfix<-cnms;varpas<-cnms[!cnms%in%active]  }
if (!is.null(active)&!is.null(include)) {varnfix<-cnms[!cnms%in%include];varpas<-cnms[!cnms%in%active]}
varpas<-varpas[!my_grepl("strata",varpas)]

stop=F

while(stop==F){

vcvm<-vcov(fit)

if (colnames(model.matrix(fit))[1]=="(Intercept)") vcvm<-vcvm[-1,-1]


if (is.null(dim(vcvm))) {
		vcvm<-matrix(vcvm,ncol=1,nrow=1)
		colnames(vcvm)<-rownames(vcvm)<-rep(attributes(fit$terms)$term.labels,table(fit$assign[-1]))
		}


if (verbose==TRUE) {
	cat("\n\nModel under investigation:\n")
	print(fit$call)
	}


 if (criterion!="alpha") bl<-drop1(fit,scope=as.formula(paste("~",paste(unique(varnfix),collapse=" + ") )),k=k) else bl<-drop1(fit,scope=as.formula(paste("~",paste(unique(varnfix),collapse=" + ") )),test=type.test)
varnfixn<-unique(varnfix)
if (verbose==TRUE)  if (criterion!="alpha")  cat("Criterion for non-passive variables: "  ,paste( paste(varnfixn,round(bl$AIC[-1],4),sep=" : "),collapse=" , "),"\n"  ,sep=""  ) else cat("Criterion for non-passive variables: "  ,paste( paste(varnfixn,round(bl[-1,pmatch("Pr",names(bl))],4),sep=" : "),collapse=" , "),"\n"  ,sep=""  )


if (criterion!="alpha") black.list.i<-varnfixn[which(bl$AIC[-1]<bl$AIC[1])]  else   black.list.i<-varnfixn[which(bl[-1,pmatch("Pr",names(bl))]>alpha)]


if (length(black.list.i)!=0){


  if (criterion=="alpha") {
      black.list.i<-varnfixn[order(-bl[-1,pmatch("Pr",names(bl))])][1:length(black.list.i)]
    criterion.i<-bl[-1,pmatch("Pr",names(bl))][order(-bl[-1,pmatch("Pr",names(bl))])][1:length(black.list.i)]
  } else {
    black.list.i<-varnfixn[order(bl$AIC[-1])][1:length(black.list.i)]
    criterion.i<-bl$AIC[-1][order(bl$AIC[-1])][1:length(black.list.i)]-bl$AIC[1]
  }





	if (verbose==TRUE) cat( "  ",paste("black list: ", paste(paste(black.list.i,round(criterion.i,4),sep=" : "),collapse =", ")),"\n" )

		flag=T
		i=1
		while (flag==T&i<=length(black.list.i)){



				xf<-as.formula(paste("~.-",black.list.i[i]  ))

			 	fit.i<-update(fit,xf,evaluate=FALSE)
				fit.i<-eval.parent(fit.i)


				  var.mod.i<-attributes(fit.i$terms)$term.labels

				  if (length(var.mod.i)==0) fit.i$assign<-0 else {
				  fit.i$assign<-list()
				  fit.i$assign[1]<-1
				  name.cfi<-names(coef(fit.i))

				  for (ii in 1:length(var.mod.i)){
				    fit.i$assign[[ii+1]]<-length(grep(paste(unlist(strsplit(paste(unlist(strsplit(var.mod.i[ii],split="\\(")),collapse="\\("),split="\\)")),collapse="\\)"),name.cfi))
				  }

				  fit.i$assign<-unlist(fit.i$assign)

				  fit.i$assign<-rep(0:length(var.mod.i),fit.i$assign)
				}


				if (colnames(model.matrix(fit))[1]=="(Intercept)") {

					change.in.estimate<-abs(fit$coef[which(rep(attributes(fit$terms)$term.labels[!my_grepl("strata",var.mod)],table(fit$assign[-1]))%in%varpas[!varpas%in%black.list.i[i]])+1]-fit.i$coef[which(!rep(attributes(fit.i$terms)$term.labels[!my_grepl("strata",var.mod.i)],table(fit.i$assign[-1]))%in%active[!active%in%black.list.i[i]])+1])
					} else {
					change.in.estimate<-abs(fit$coef[which(rep(attributes(fit$terms)$term.labels[!my_grepl("strata",var.mod)],table(fit$assign[-1]))%in%varpas[!varpas%in%black.list.i[i]])]-fit.i$coef[which(!rep(attributes(fit.i$terms)$term.labels[!my_grepl("strata",var.mod.i)],table(fit.i$assign[-1]))%in%active[!active%in%black.list.i[i]])])
					}




			if (exp.beta==TRUE) if (length(varpas[!varpas%in%black.list.i[i]])>1) ch.in.est<-exp( change.in.estimate*sqrt(diag(vcvx[rownames(vcvx)%in%varpas[!varpas%in%black.list.i[i]],colnames(vcvx)%in%varpas[!varpas%in%black.list.i[i]]]))  ) else ch.in.est<-exp( change.in.estimate*sqrt((vcvx[rownames(vcvx)%in%varpas[!varpas%in%black.list.i[i]],colnames(vcvx)%in%varpas[!varpas%in%black.list.i[i]]])) )
			if (exp.beta==FALSE) if (length(varpas[!varpas%in%black.list.i[i]])>1) ch.in.est<-change.in.estimate*sqrt(diag(vcvx[rownames(vcvx)%in%varpas[!varpas%in%black.list.i[i]],colnames(vcvx)%in%varpas[!varpas%in%black.list.i[i]]]))/sd(fit$y) else ch.in.est<-change.in.estimate*sqrt((vcvx[rownames(vcvx)%in%varpas[!varpas%in%black.list.i[i]],colnames(vcvx)%in%varpas[!varpas%in%black.list.i[i]]]))/sd(fit$y)

			if (verbose==TRUE) cat("          ",paste("Investigating change in b or exp(b) due to omitting variable ",black.list.i[i]," ; ",paste(paste(names(ch.in.est),round(ch.in.est,4), sep=" : ") ,collapse=", ") )  ,"\n")

			if (exp.beta==TRUE) if( sum(ch.in.est >= 1+tau)==0) {
						flag=F
						varnfix<-varnfix[-which(varnfix==black.list.i[i])]
						if(!black.list.i[i]%in%include) if (sum(varpas==black.list.i[i])==0) varpas<-varpas else varpas<-varpas[-which(varpas==black.list.i[i])]
						}

			if (exp.beta==FALSE) if (sum(  ch.in.est  >=tau  )==0)  {
						flag=F
						varnfix<-varnfix[-which(varnfix==black.list.i[i])]
						if(!black.list.i[i]%in%include) if (sum(varpas==black.list.i[i])==0) varpas<-varpas else  varpas<-varpas[-which(varpas==black.list.i[i])]
						}

   			if (flag==F) {
   			 	  fit<-update(fit,as.formula(paste("~.-",black.list.i[i])),evaluate=FALSE)
   			  fit<-eval.parent(fit)


   			    var.mod<-attributes(fit$terms)$term.labels
   			    if (length(var.mod)==0) fit$assign<-0 else {
   			    fit$assign<-list()
   			    fit$assign[1]<-1
   			    name.cf<-names(coef(fit))

   			    for (i in 1:length(var.mod)){
   			      fit$assign[[i+1]]<-length(grep(paste(unlist(strsplit(paste(unlist(strsplit(var.mod[i],split="\\(")),collapse="\\("),split="\\)")),collapse="\\)"),name.cf))
   			    }

   			    fit$assign<-unlist(fit$assign)

   			    fit$assign<-rep(0:length(var.mod),fit$assign)
   			  }

   			}

			if (verbose==TRUE)  if (flag==T&i!=length(black.list.i))  cat("  ","updated black list:" , paste(paste(black.list.i[-(1:i)],round(criterion.i[-(1:i)],4),sep=" : "),collapse =", "),"\n" )

			i=i+1
			}


	} else {
	if (verbose==T) cat("black list: empty","\n")
	flag=T
	}

if (length(black.list.i)==0) stop=T
if (length(varnfix)==0) stop=T
if (length(include)>0) if (length(varpas)==length(include)& ( length(varnfix[!varnfix%in%varpas])==0 )) stop=T
if (flag==T) stop=T

}

if (verbose==T) {
	cat("\n\nFinal model:\n")
	print(fit$call)
	cat("\n\n")
	}

fit

}





#' ABE for model which includes categorical covariates, factor option, bootstrap version
#' @keywords internal
#' @examples
#' \dontshow{
#' set.seed(1)
#' n=100
#' x1<-runif(n)
#' x2<-runif(n)
#' x3<-runif(n)
#' y<--5+5*x1+5*x2+ rnorm(n,sd=5)
#' dd<-data.frame(y,x1,x2,x3)
#' fit<-lm(y~x1+x2+x3,x=TRUE,y=TRUE,data=dd)
#'
#' abe.fit<-abe.fact1.boot(fit,data=dd,include="x1",active="x2",
#' tau=0.05,exp.beta=FALSE,exact=TRUE,criterion="alpha",alpha=0.2,
#' type.test="Chisq",k=2)
#' summary(abe.fit)
#' }


abe.fact1.boot<-function(fit,data,include=NULL,active=NULL,tau=0.05,exp.beta=TRUE,exact=FALSE,criterion="alpha",alpha=0.2,type.test="Chisq",k){


  if (colnames(model.matrix(fit))[1]=="(Intercept)") xm<-as.matrix(fit$x)[,-1] else xm<-as.matrix(fit$x)


  var.mod<-attributes(fit$terms)$term.labels
  if (length(var.mod)==0) fit<-0 else {
    fit$assign<-list()
    fit$assign[1]<-1
    name.cf<-names(coef(fit))

    for (i in 1:length(var.mod)){
      fit$assign[[i+1]]<-length(grep(paste(unlist(strsplit(paste(unlist(strsplit(var.mod[i],split="\\(")),collapse="\\("),split="\\)")),collapse="\\)"),name.cf))
    }

    fit$assign<-unlist(fit$assign)

    fit$assign<-rep(0:length(var.mod),fit$assign)
  }
   cnms<-rep(attributes(fit$terms)$term.labels[!my_grepl("strata",var.mod)],table(fit$assign[-1]))


  vcvx<-var(xm)
  rownames(vcvx)<-colnames(vcvx)<-cnms

  if (sum(my_grepl("strata",var.mod))!=0) cnms<-c(cnms,var.mod[my_grepl("strata",var.mod)])



  if (!is.null(include)) {
    include.l<-list()
    for (i in 1:length(include)) {
      include.l[[i]]<-cnms[my_grep(include[i],cnms)]
    }
    include<-unlist(include.l)
  }
  if (!is.null(active))  {
    active.l<-list()
    for (i in 1:length(active)) {
      active.l[[i]]<-cnms[my_grep(active[i],cnms)]
    }
    active<-unlist(active.l)
  }




  if (is.null(include)&is.null(active)) {varnfix<- cnms;varpas<-cnms}
  if (is.null(active)&!is.null(include)) {varnfix<- cnms[!cnms%in%include]; varpas<-cnms}
  if (!is.null(active)&is.null(include))  {varnfix<-cnms;varpas<-cnms[!cnms%in%active]  }
  if (!is.null(active)&!is.null(include)) {varnfix<-cnms[!cnms%in%include];varpas<-cnms[!cnms%in%active]}
  varpas<-varpas[!my_grepl("strata",varpas)]


  stop=F

  while(stop==F){

    vcvm<-vcov(fit)

    if (colnames(model.matrix(fit))[1]=="(Intercept)") vcvm<-vcvm[-1,-1]



    if (is.null(dim(vcvm))) {
      vcvm<-matrix(vcvm,ncol=1,nrow=1)
      colnames(vcvm)<-rownames(vcvm)<-rep(attributes(fit$terms)$term.labels,table(fit$assign[-1]))
    }




     if (criterion!="alpha") bl<-drop1(fit,scope=as.formula(paste("~",paste(unique(varnfix),collapse=" + ") )),k=k) else bl<-drop1(fit,scope=as.formula(paste("~",paste(unique(varnfix),collapse=" + ") )),test=type.test)
    varnfixn<-unique(varnfix)


    if (criterion!="alpha") black.list.i<-varnfixn[which(bl$AIC[-1]<bl$AIC[1])]  else   black.list.i<-varnfixn[which(bl[-1,pmatch("Pr",names(bl))]>alpha)]


    if (length(black.list.i)!=0){

      if (criterion=="alpha") {
         black.list.i<-varnfixn[order(-bl[-1,pmatch("Pr",names(bl))])][1:length(black.list.i)]
        criterion.i<-bl[-1,pmatch("Pr",names(bl))][order(-bl[-1,pmatch("Pr",names(bl))])][1:length(black.list.i)]
      } else {
        black.list.i<-varnfixn[order(bl$AIC[-1])][1:length(black.list.i)]
        criterion.i<-bl$AIC[-1][order(bl$AIC[-1])][1:length(black.list.i)]-bl$AIC[1]
      }





         flag=T
      i=1
      while (flag==T&i<=length(black.list.i)){



          xf<-as.formula(paste("~.-",black.list.i[i]  ))

          fit.i<-update(fit,xf,evaluate=FALSE)
        fit.i<-eval.parent(fit.i)


        var.mod.i<-attributes(fit.i$terms)$term.labels

        if (length(var.mod.i)==0) fit.i$assign<-0 else {
          fit.i$assign<-list()
          fit.i$assign[1]<-1
          name.cfi<-names(coef(fit.i))

          for (ii in 1:length(var.mod.i)){
            fit.i$assign[[ii+1]]<-length(grep(paste(unlist(strsplit(paste(unlist(strsplit(var.mod.i[ii],split="\\(")),collapse="\\("),split="\\)")),collapse="\\)"),name.cfi))
          }

          fit.i$assign<-unlist(fit.i$assign)

          fit.i$assign<-rep(0:length(var.mod.i),fit.i$assign)
        }


        if (colnames(model.matrix(fit))[1]=="(Intercept)") {

          change.in.estimate<-abs(fit$coef[which(rep(attributes(fit$terms)$term.labels[!my_grepl("strata",var.mod)],table(fit$assign[-1]))%in%varpas[!varpas%in%black.list.i[i]])+1]-fit.i$coef[which(!rep(attributes(fit.i$terms)$term.labels[!my_grepl("strata",var.mod.i)],table(fit.i$assign[-1]))%in%active[!active%in%black.list.i[i]])+1])
        } else {
          change.in.estimate<-abs(fit$coef[which(rep(attributes(fit$terms)$term.labels[!my_grepl("strata",var.mod)],table(fit$assign[-1]))%in%varpas[!varpas%in%black.list.i[i]])]-fit.i$coef[which(!rep(attributes(fit.i$terms)$term.labels[!my_grepl("strata",var.mod.i)],table(fit.i$assign[-1]))%in%active[!active%in%black.list.i[i]])])
        }


        if (exp.beta==TRUE) if (length(varpas[!varpas%in%black.list.i[i]])>1) ch.in.est<-exp( change.in.estimate*sqrt(diag(vcvx[rownames(vcvx)%in%varpas[!varpas%in%black.list.i[i]],colnames(vcvx)%in%varpas[!varpas%in%black.list.i[i]]]))  ) else ch.in.est<-exp( change.in.estimate*sqrt((vcvx[rownames(vcvx)%in%varpas[!varpas%in%black.list.i[i]],colnames(vcvx)%in%varpas[!varpas%in%black.list.i[i]]])) )
        if (exp.beta==FALSE) if (length(varpas[!varpas%in%black.list.i[i]])>1) ch.in.est<-change.in.estimate*sqrt(diag(vcvx[rownames(vcvx)%in%varpas[!varpas%in%black.list.i[i]],colnames(vcvx)%in%varpas[!varpas%in%black.list.i[i]]]))/sd(fit$y) else ch.in.est<-change.in.estimate*sqrt((vcvx[rownames(vcvx)%in%varpas[!varpas%in%black.list.i[i]],colnames(vcvx)%in%varpas[!varpas%in%black.list.i[i]]]))/sd(fit$y)


        if (exp.beta==TRUE) if( sum(ch.in.est >= 1+tau)==0) {
          flag=F
          varnfix<-varnfix[-which(varnfix==black.list.i[i])]
          if(!black.list.i[i]%in%include) if (sum(varpas==black.list.i[i])==0) varpas<-varpas else varpas<-varpas[-which(varpas==black.list.i[i])]
        }

        if (exp.beta==FALSE) if (sum(  ch.in.est  >=tau  )==0)  {
          flag=F
          varnfix<-varnfix[-which(varnfix==black.list.i[i])]
          if(!black.list.i[i]%in%include) if (sum(varpas==black.list.i[i])==0) varpas<-varpas else  varpas<-varpas[-which(varpas==black.list.i[i])]
        }

        if (flag==F) {
            fit<-update(fit,as.formula(paste("~.-",black.list.i[i])),evaluate=FALSE)
          fit<-eval.parent(fit)


          var.mod<-attributes(fit$terms)$term.labels
          if (length(var.mod)==0) fit$assign<-0 else {
            fit$assign<-list()
            fit$assign[1]<-1
            name.cf<-names(coef(fit))

            for (i in 1:length(var.mod)){
              fit$assign[[i+1]]<-length(grep(paste(unlist(strsplit(paste(unlist(strsplit(var.mod[i],split="\\(")),collapse="\\("),split="\\)")),collapse="\\)"),name.cf))
            }

            fit$assign<-unlist(fit$assign)

            fit$assign<-rep(0:length(var.mod),fit$assign)
          }

        }


        i=i+1
      }


    } else {
       flag=T
    }

    if (length(black.list.i)==0) stop=T
    if (length(varnfix)==0) stop=T
    if (length(include)>0) if (length(varpas)==length(include)& ( length(varnfix[!varnfix%in%varpas])==0 )) stop=T
    if (flag==T) stop=T

  }



  fit

}







#' ABE for model which includes categorical covariates, individual option
#' @keywords internal
#' @examples
#' \dontshow{
#' set.seed(1)
#' n=100
#' x1<-runif(n)
#' x2<-runif(n)
#' x3<-runif(n)
#' y<--5+5*x1+5*x2+ rnorm(n,sd=5)
#' dd<-data.frame(y,x1,x2,x3)
#' fit<-lm(y~x1+x2+x3,x=TRUE,y=TRUE,data=dd)
#'
#' abe.fit<-abe.fact2(fit,data=dd,include="x1",active="x2",
#' tau=0.05,exp.beta=FALSE,exact=TRUE,criterion="alpha",alpha=0.2,
#' type.test="Chisq",verbose=FALSE)
#' summary(abe.fit)
#' }


abe.fact2<-function(fit,data,include=NULL,active=NULL,tau=0.05,exp.beta=TRUE,exact=FALSE,criterion="alpha",alpha=0.2,type.test="Chisq",verbose=T){


df<-as.data.frame(model.matrix(fit))

names(df)<-gsub("factor", replacement="", names(df), ignore.case = FALSE, perl = FALSE,
     fixed = FALSE, useBytes = FALSE)

names(df)<-gsub(")", replacement=".", names(df), ignore.case = FALSE, perl = FALSE,
     fixed = FALSE, useBytes = FALSE)

names(df)<-gsub("\\(", replacement="", names(df), ignore.case = FALSE, perl = FALSE,
     fixed = FALSE, useBytes = FALSE)

check.names<-names(df)
var.mod<-attributes(fit$terms)$term.labels
if (sum(my_grepl("strata",var.mod))!=0)  check.names<-c(check.names, var.mod[my_grepl("strata",var.mod)] )



if (!is.null(include)) {
  include.l<-list()
  for (i in 1:length(include)) {
    if (sum(my_grepl(include[i],check.names))==0) stop("at least one include variable is not in the model")
    include.l[[i]]<-check.names[my_grep(include[i],check.names)]
  }
  include<-unlist(include.l)
}
if (!is.null(active))  {
  active.l<-list()
  for (i in 1:length(active)) {
    if (sum(my_grepl(active[i],check.names))==0) stop("at least one active variable is not in the model")
    active.l[[i]]<-check.names[my_grep(active[i],check.names)]
  }
  active<-unlist(active.l)
}



if ( colnames(model.matrix(fit))[1]=="(Intercept)" )   updt.f<-as.formula(paste("~",paste(check.names[-1],collapse="+"))) else  updt.f<-as.formula(paste("~",paste(check.names,collapse="+"),"-1"))
if ( class(fit)[1]=="coxph" )    updt.f<-as.formula(paste("~",paste(check.names,collapse="+")))

df<-cbind(df,model.frame(fit),data)

fit<-my_update(fit, updt.f   ,data=df)

if (sum(include%in%active)!=0) stop("at least one include variable is also specified as active")

 if (criterion[1]=="AIC") k<-2

n<-nrow(fit$x)

if (criterion[1]=="BIC") k<-log(n)


 if (colnames(model.matrix(fit))[1]=="(Intercept)") xm<-as.matrix(fit$x)[,-1] else xm<-as.matrix(fit$x)


 vcvx<-var(xm)
rownames(vcvx)<-colnames(vcvx)<-attributes(fit$terms)$term.labels[!my_grepl("strata",attributes(fit$terms)$term.labels)]

cnms<-attributes(fit$terms)$term.labels

if ( sum(cnms%in%include)==length(cnms) ) stop("all variables are specified as pasive, cannot perform variable selection")
if ( sum(cnms%in%active)==length(cnms) ) {
		include=active=NULL
		warning("all variables are specified as active, treating all variables as active or pasive")
		}


if (is.null(include)&is.null(active)) {varnfix<- cnms;varpas<-cnms}
if (is.null(active)&!is.null(include)) {varnfix<- cnms[!cnms%in%include]; varpas<-cnms}
if (!is.null(active)&is.null(include))  {varnfix<-cnms;varpas<-cnms[!cnms%in%active]  }
if (!is.null(active)&!is.null(include)) {varnfix<-cnms[!cnms%in%include];varpas<-cnms[!cnms%in%active]}
varpas<-varpas[!my_grepl("strata",varpas)]


stop=F

while(stop==F){

vcvm<-vcov(fit)

if (colnames(model.matrix(fit))[1]=="(Intercept)") vcvm<-vcvm[-1,-1]
if (is.null(dim(vcvm))) {
		vcvm<-matrix(vcvm,ncol=1,nrow=1)
		if (colnames(model.matrix(fit))[1]=="(Intercept)") colnames(vcvm)<-rownames(vcvm)<-colnames(vcov(fit))[-1] else colnames(vcvm)<-rownames(vcvm)<-colnames(vcov(fit))
		}


if (verbose==TRUE) {
	cat("\n\nModel under investigation:\n")
 	print(formula(fit))
	}


 if (criterion!="alpha") bl<-drop1(fit,scope=as.formula(paste("~",paste(varnfix,collapse=" + ") )),k=k) else bl<-drop1(fit,scope=as.formula(paste("~",paste(varnfix,collapse=" + ") )),test=type.test)

if (verbose==TRUE)  if (criterion!="alpha")  cat("Criterion for non-passive variables: "  ,paste( paste(varnfix,round(bl$AIC[-1],4),sep=" : "),collapse=" , "),"\n"  ,sep=""  ) else cat("Criterion for non-passive variables: "  ,paste( paste(varnfix,round(bl[-1,pmatch("Pr",names(bl))],4),sep=" : "),collapse=" , "),"\n"  ,sep=""  )

if (criterion!="alpha") black.list.i<-varnfix[which(bl$AIC[-1]<bl$AIC[1])]  else   black.list.i<-varnfix[which(bl[-1,pmatch("Pr",names(bl))]>alpha)]


if (length(black.list.i)!=0){

  if (criterion=="alpha") {
     black.list.i<-varnfix[order(-bl[-1,pmatch("Pr",names(bl))])][1:length(black.list.i)]
    criterion.i<-bl[-1,pmatch("Pr",names(bl))][order(-bl[-1,pmatch("Pr",names(bl))])][1:length(black.list.i)]
  } else {
    black.list.i<-varnfix[order(bl$AIC[-1])][1:length(black.list.i)]
    criterion.i<-bl$AIC[-1][order(bl$AIC[-1])][1:length(black.list.i)]-bl$AIC[1]
  }



	if (verbose==TRUE) cat( "  ",paste("black list: ", paste(paste(black.list.i,round(criterion.i,4),sep=" : "),collapse =", ")),"\n" )

		flag=T
		i=1
		while (flag==T&i<=length(black.list.i)){



			if (exact==T){

				xf<-as.formula(paste("~.-",black.list.i[i]  ))

			 	fit.i<-my_update(fit,xf,data=df)

				if (colnames(model.matrix(fit))[1]=="(Intercept)") {

					change.in.estimate<-abs(fit$coef[which(attributes(fit$terms)$term.labels[!my_grepl("strata",attributes(fit$terms)$term.labels)]%in%varpas[!varpas%in%black.list.i[i]])+1]-fit.i$coef[which(!attributes(fit.i$terms)$term.labels[!my_grepl("strata",attributes(fit.i$terms)$term.labels)]%in%active[!active%in%black.list.i[i]])+1])
					} else {
					change.in.estimate<-abs(fit$coef[which(attributes(fit$terms)$term.labels[!my_grepl("strata",attributes(fit$terms)$term.labels)]%in%varpas[!varpas%in%black.list.i[i]])]-fit.i$coef[which(!attributes(fit.i$terms)$term.labels[!my_grepl("strata",attributes(fit.i$terms)$term.labels)]%in%active[!active%in%black.list.i[i]])])
					}


				} else {
				if (colnames(model.matrix(fit))[1]=="(Intercept)") {
					change.in.estimate<- abs(fit$coef[which(attributes(fit$terms)$term.labels[!my_grepl("strata",attributes(fit$terms)$term.labels)]==black.list.i[i])+1]*vcvm[rownames(vcvm)==black.list.i[i],colnames(vcvm)%in%varpas[!varpas%in%black.list.i[i]]] /vcvm[rownames(vcvm)==black.list.i[i],colnames(vcvm)==black.list.i[i]])
					names(change.in.estimate)<-colnames(vcvm)[colnames(vcvm)%in%varpas[!varpas%in%black.list.i[i]]]
					} else {
					change.in.estimate<- abs(fit$coef[which(attributes(fit$terms)$term.labels[!my_grepl("strata",attributes(fit$terms)$term.labels)]==black.list.i[i])]*vcvm[rownames(vcvm)==black.list.i[i],colnames(vcvm)%in%varpas[!varpas%in%black.list.i[i]]] /vcvm[rownames(vcvm)==black.list.i[i],colnames(vcvm)==black.list.i[i]])
					names(change.in.estimate)<-colnames(vcvm)[colnames(vcvm)%in%varpas[!varpas%in%black.list.i[i]]]
					}
				}


			if (exp.beta==TRUE) if (length(varpas[!varpas%in%black.list.i[i]])>1) ch.in.est<-exp( change.in.estimate*sqrt(diag(vcvx[rownames(vcvx)%in%varpas[!varpas%in%black.list.i[i]],colnames(vcvx)%in%varpas[!varpas%in%black.list.i[i]]]))  ) else ch.in.est<-exp( change.in.estimate*sqrt((vcvx[rownames(vcvx)%in%varpas[!varpas%in%black.list.i[i]],colnames(vcvx)%in%varpas[!varpas%in%black.list.i[i]]])) )
			if (exp.beta==FALSE) if (length(varpas[!varpas%in%black.list.i[i]])>1) ch.in.est<-change.in.estimate*sqrt(diag(vcvx[rownames(vcvx)%in%varpas[!varpas%in%black.list.i[i]],colnames(vcvx)%in%varpas[!varpas%in%black.list.i[i]]]))/sd(fit$y) else ch.in.est<-change.in.estimate*sqrt((vcvx[rownames(vcvx)%in%varpas[!varpas%in%black.list.i[i]],colnames(vcvx)%in%varpas[!varpas%in%black.list.i[i]]]))/sd(fit$y)

			if (verbose==TRUE) cat("          ",paste("Investigating change in b or exp(b) due to omitting variable ",black.list.i[i]," ; ",paste(paste(names(ch.in.est),round(ch.in.est,4), sep=" : ") ,collapse=", ") )  ,"\n")

			if (exp.beta==TRUE) if( sum(ch.in.est >= 1+tau)==0) {
						flag=F
						varnfix<-varnfix[-which(varnfix==black.list.i[i])]
						if(!black.list.i[i]%in%include) if (sum(varpas==black.list.i[i])==0) varpas<-varpas else varpas<-varpas[-which(varpas==black.list.i[i])]
						}

			if (exp.beta==FALSE) if (sum(  ch.in.est  >=tau  )==0)  {
						flag=F
						varnfix<-varnfix[-which(varnfix==black.list.i[i])]
						if(!black.list.i[i]%in%include) if (sum(varpas==black.list.i[i])==0) varpas<-varpas else  varpas<-varpas[-which(varpas==black.list.i[i])]
						}

   			if (flag==F) {
   			   fit<-my_update(fit,as.formula(paste("~.-",black.list.i[i])),data=df)
   			    }

			if (verbose==TRUE)  if (flag==T&i!=length(black.list.i))  cat("  ","updated black list:" , paste(paste(black.list.i[-(1:i)],round(criterion.i[-(1:i)],4),sep=" : "),collapse =", "),"\n" )

			i=i+1
			}


	} else {
	if (verbose==T) cat("black list: empty","\n")
	flag=T
	}

if (length(black.list.i)==0) stop=T
if (length(varnfix)==0) stop=T
if (length(include)>0) if (length(varpas)==length(include)& ( length(varnfix[!varnfix%in%varpas])==0 )) stop=T
if (flag==T) stop=T

}

if (verbose==T) {
	cat("\n\nFinal model:\n")
	print(formula(fit))
	cat("\n\n")
	}

fit<-my_update2(fit,data.n="df")
fit
}








#' ABE for model which includes categorical covariates, individual option, bootstrap version
#' @keywords internal
#' @examples
#' \dontshow{
#' set.seed(1)
#' n=100
#' x1<-runif(n)
#' x2<-runif(n)
#' x3<-runif(n)
#' y<--5+5*x1+5*x2+ rnorm(n,sd=5)
#' dd<-data.frame(y,x1,x2,x3)
#' fit<-lm(y~x1+x2+x3,x=TRUE,y=TRUE,data=dd)
#'
#' abe.fit<-abe.fact2.boot(fit,data=dd,include="x1",active="x2",
#' tau=0.05,exp.beta=FALSE,exact=TRUE,criterion="alpha",alpha=0.2,
#' type.test="Chisq",k=2)
#' summary(abe.fit)
#' }


abe.fact2.boot<-function(fit,data,include=NULL,active=NULL,tau=0.05,exp.beta=TRUE,exact=FALSE,criterion="alpha",alpha=0.2,type.test="Chisq",k){


  df<-as.data.frame(model.matrix(fit))

  names(df)<-gsub("factor", replacement="", names(df), ignore.case = FALSE, perl = FALSE,
                  fixed = FALSE, useBytes = FALSE)

  names(df)<-gsub(")", replacement=".", names(df), ignore.case = FALSE, perl = FALSE,
                  fixed = FALSE, useBytes = FALSE)

  names(df)<-gsub("\\(", replacement="", names(df), ignore.case = FALSE, perl = FALSE,
                  fixed = FALSE, useBytes = FALSE)

  check.names<-names(df)
  var.mod<-attributes(fit$terms)$term.labels
  if (sum(my_grepl("strata",var.mod))!=0)  check.names<-c(check.names, var.mod[my_grepl("strata",var.mod)] )



  if (!is.null(include)) {
    include.l<-list()
    for (i in 1:length(include)) {
       include.l[[i]]<-check.names[my_grep(include[i],check.names)]
    }
    include<-unlist(include.l)
  }
  if (!is.null(active))  {
    active.l<-list()
    for (i in 1:length(active)) {
       active.l[[i]]<-check.names[my_grep(active[i],check.names)]
    }
    active<-unlist(active.l)
  }

  if ( colnames(model.matrix(fit))[1]=="(Intercept)" )   updt.f<-as.formula(paste("~",paste(check.names[-1],collapse="+"))) else  updt.f<-as.formula(paste("~",paste(check.names,collapse="+"),"-1"))
  if ( class(fit)[1]=="coxph" )    updt.f<-as.formula(paste("~",paste(check.names,collapse="+")))

  df<-cbind(df,model.frame(fit),data)



  fit<-my_update(fit, updt.f   ,data=df)




  if (colnames(model.matrix(fit))[1]=="(Intercept)") xm<-as.matrix(fit$x)[,-1] else xm<-as.matrix(fit$x)


  vcvx<-var(xm)
  rownames(vcvx)<-colnames(vcvx)<-attributes(fit$terms)$term.labels[!my_grepl("strata",attributes(fit$terms)$term.labels)]


  cnms<-attributes(fit$terms)$term.labels


  if (is.null(include)&is.null(active)) {varnfix<- cnms;varpas<-cnms}
  if (is.null(active)&!is.null(include)) {varnfix<- cnms[!cnms%in%include]; varpas<-cnms}
  if (!is.null(active)&is.null(include))  {varnfix<-cnms;varpas<-cnms[!cnms%in%active]  }
  if (!is.null(active)&!is.null(include)) {varnfix<-cnms[!cnms%in%include];varpas<-cnms[!cnms%in%active]}
  varpas<-varpas[!my_grepl("strata",varpas)]


  stop=F

  while(stop==F){

    vcvm<-vcov(fit)

    if (colnames(model.matrix(fit))[1]=="(Intercept)") vcvm<-vcvm[-1,-1]
    if (is.null(dim(vcvm))) {
      vcvm<-matrix(vcvm,ncol=1,nrow=1)
      if (colnames(model.matrix(fit))[1]=="(Intercept)") colnames(vcvm)<-rownames(vcvm)<-colnames(vcov(fit))[-1] else colnames(vcvm)<-rownames(vcvm)<-colnames(vcov(fit))
    }




     if (criterion!="alpha") bl<-drop1(fit,scope=as.formula(paste("~",paste(varnfix,collapse=" + ") )),k=k) else bl<-drop1(fit,scope=as.formula(paste("~",paste(varnfix,collapse=" + ") )),test=type.test)


    if (criterion!="alpha") black.list.i<-varnfix[which(bl$AIC[-1]<bl$AIC[1])]  else   black.list.i<-varnfix[which(bl[-1,pmatch("Pr",names(bl))]>alpha)]


    if (length(black.list.i)!=0){

      if (criterion=="alpha") {
         black.list.i<-varnfix[order(-bl[-1,pmatch("Pr",names(bl))])][1:length(black.list.i)]
        criterion.i<-bl[-1,pmatch("Pr",names(bl))][order(-bl[-1,pmatch("Pr",names(bl))])][1:length(black.list.i)]
      } else {
        black.list.i<-varnfix[order(bl$AIC[-1])][1:length(black.list.i)]
        criterion.i<-bl$AIC[-1][order(bl$AIC[-1])][1:length(black.list.i)]-bl$AIC[1]
      }



       flag=T
      i=1
      while (flag==T&i<=length(black.list.i)){



        if (exact==T){

          xf<-as.formula(paste("~.-",black.list.i[i]  ))

            fit.i<-update(fit,xf,evaluate=FALSE)
          fit.i<-eval.parent(fit.i)

          if (colnames(model.matrix(fit))[1]=="(Intercept)") {

            change.in.estimate<-abs(fit$coef[which(attributes(fit$terms)$term.labels[!my_grepl("strata",attributes(fit$terms)$term.labels)]%in%varpas[!varpas%in%black.list.i[i]])+1]-fit.i$coef[which(!attributes(fit.i$terms)$term.labels[!my_grepl("strata",attributes(fit.i$terms)$term.labels)]%in%active[!active%in%black.list.i[i]])+1])
          } else {
            change.in.estimate<-abs(fit$coef[which(attributes(fit$terms)$term.labels[!my_grepl("strata",attributes(fit$terms)$term.labels)]%in%varpas[!varpas%in%black.list.i[i]])]-fit.i$coef[which(!attributes(fit.i$terms)$term.labels[!my_grepl("strata",attributes(fit.i$terms)$term.labels)]%in%active[!active%in%black.list.i[i]])])
          }


        } else {
          if (colnames(model.matrix(fit))[1]=="(Intercept)") {
            change.in.estimate<- abs(fit$coef[which(attributes(fit$terms)$term.labels[!my_grepl("strata",attributes(fit$terms)$term.labels)]==black.list.i[i])+1]*vcvm[rownames(vcvm)==black.list.i[i],colnames(vcvm)%in%varpas[!varpas%in%black.list.i[i]]] /vcvm[rownames(vcvm)==black.list.i[i],colnames(vcvm)==black.list.i[i]])
            names(change.in.estimate)<-colnames(vcvm)[colnames(vcvm)%in%varpas[!varpas%in%black.list.i[i]]]
          } else {
            change.in.estimate<- abs(fit$coef[which(attributes(fit$terms)$term.labels[!my_grepl("strata",attributes(fit$terms)$term.labels)]==black.list.i[i])]*vcvm[rownames(vcvm)==black.list.i[i],colnames(vcvm)%in%varpas[!varpas%in%black.list.i[i]]] /vcvm[rownames(vcvm)==black.list.i[i],colnames(vcvm)==black.list.i[i]])
            names(change.in.estimate)<-colnames(vcvm)[colnames(vcvm)%in%varpas[!varpas%in%black.list.i[i]]]
          }
        }


        if (exp.beta==TRUE) if (length(varpas[!varpas%in%black.list.i[i]])>1) ch.in.est<-exp( change.in.estimate*sqrt(diag(vcvx[rownames(vcvx)%in%varpas[!varpas%in%black.list.i[i]],colnames(vcvx)%in%varpas[!varpas%in%black.list.i[i]]]))  ) else ch.in.est<-exp( change.in.estimate*sqrt((vcvx[rownames(vcvx)%in%varpas[!varpas%in%black.list.i[i]],colnames(vcvx)%in%varpas[!varpas%in%black.list.i[i]]])) )
        if (exp.beta==FALSE) if (length(varpas[!varpas%in%black.list.i[i]])>1) ch.in.est<-change.in.estimate*sqrt(diag(vcvx[rownames(vcvx)%in%varpas[!varpas%in%black.list.i[i]],colnames(vcvx)%in%varpas[!varpas%in%black.list.i[i]]]))/sd(fit$y) else ch.in.est<-change.in.estimate*sqrt((vcvx[rownames(vcvx)%in%varpas[!varpas%in%black.list.i[i]],colnames(vcvx)%in%varpas[!varpas%in%black.list.i[i]]]))/sd(fit$y)


        if (exp.beta==TRUE) if( sum(ch.in.est >= 1+tau)==0) {
          flag=F
          varnfix<-varnfix[-which(varnfix==black.list.i[i])]
          if(!black.list.i[i]%in%include) if (sum(varpas==black.list.i[i])==0) varpas<-varpas else varpas<-varpas[-which(varpas==black.list.i[i])]
        }

        if (exp.beta==FALSE) if (sum(  ch.in.est  >=tau  )==0)  {
          flag=F
          varnfix<-varnfix[-which(varnfix==black.list.i[i])]
          if(!black.list.i[i]%in%include) if (sum(varpas==black.list.i[i])==0) varpas<-varpas else  varpas<-varpas[-which(varpas==black.list.i[i])]
        }

        if (flag==F) {
            fit<-update(fit,as.formula(paste("~.-",black.list.i[i])),evaluate=FALSE)
          fit<-eval.parent(fit)
        }


        i=i+1
      }


    } else {
        flag=T
    }

    if (length(black.list.i)==0) stop=T
    if (length(varnfix)==0) stop=T
    if (length(include)>0) if (length(varpas)==length(include)& ( length(varnfix[!varnfix%in%varpas])==0 )) stop=T
    if (flag==T) stop=T

  }


  fit


}



#' update function which searches for objects within the parent environment, gives a nicer output than my_update
#' @keywords internal
#' @examples
#' \dontshow{
#' set.seed(1)
#' n=100
#' x1<-runif(n)
#' x2<-runif(n)
#' x3<-runif(n)
#' y<--5+5*x1+5*x2+ rnorm(n,sd=5)
#' dd<-data.frame(y,x1,x2,x3)
#' fit<-lm(y~x1+x2+x3,x=TRUE,y=TRUE,data=dd)
#'
#' ddn<-dd[-1,]
#' my_update2(fit,data=ddn,data.n="ddn")
#' my_update2(fit,formula=as.formula(".~.-x1"),data=ddn,data.n="ddn")
#' }



my_update2 <- function(mod, formula = NULL, data = NULL,data.n=NULL) {
  call <- getCall(mod)
  if (is.null(call)) {
    stop("Model object does not support updating (no call)", call. = FALSE)
  }
  term <- terms(mod)
  if (is.null(term)) {
    stop("Model object does not support updating (no terms)", call. = FALSE)
  }

  if (!is.null(data)) call$data <- data
  if (!is.null(formula)) call$formula <- update.formula(call$formula, formula)
  env <- attr(term, ".Environment")

  fit<-eval(call, env, parent.frame())
  if (!is.null(data.n)) fit$call$data<-as.symbol(data.n)
  fit
}

#' update function which searches for objects within the parent environment
#' @keywords internal
#' @examples
#' \dontshow{
#' set.seed(1)
#' n=100
#' x1<-runif(n)
#' x2<-runif(n)
#' x3<-runif(n)
#' y<--5+5*x1+5*x2+ rnorm(n,sd=5)
#' dd<-data.frame(y,x1,x2,x3)
#' fit<-lm(y~x1+x2+x3,x=TRUE,y=TRUE,data=dd)
#'
#' ddn<-dd[-1,]
#' my_update(fit,data=ddn)
#' my_update(fit,formula=as.formula(".~.-x1"),data=ddn)
#' }


my_update <- function(mod, formula = NULL, data = NULL) {
  call <- getCall(mod)
  if (is.null(call)) {
    stop("Model object does not support updating (no call)", call. = FALSE)
  }
  term <- terms(mod)
  if (is.null(term)) {
    stop("Model object does not support updating (no terms)", call. = FALSE)
  }

  if (!is.null(data)) call$data <- data
  if (!is.null(formula)) call$formula <- update.formula(call$formula, formula)
  env <- attr(term, ".Environment")

  eval(call, env, parent.frame())

}


#' update function which searches for objects within the parent environment, bootstrap version, i.e. can only update the model based on a new dataset
#' @keywords internal
#' @examples
#' \dontshow{
#' set.seed(1)
#' n=100
#' x1<-runif(n)
#' x2<-runif(n)
#' x3<-runif(n)
#' y<--5+5*x1+5*x2+ rnorm(n,sd=5)
#' dd<-data.frame(y,x1,x2,x3)
#' fit<-lm(y~x1+x2+x3,x=TRUE,y=TRUE,data=dd)
#'
#' ddn<-dd[-1,]
#' my_update_boot(fit,data=ddn)
#' }

my_update_boot <- function(mod, data = NULL) {
   call <- getCall(mod)

   term <- terms(mod)


    call$data <- data
    call$formula <-  call$formula
   env <- attr(term, ".Environment")

   eval(call, env, parent.frame())
}


#' grepl function changed
#' @keywords internal
#' @examples
#' \dontshow{
#' my_grepl("x",c("xy","xz","ab"))
#' }
my_grepl<-function(...) grepl(fixed=TRUE,...)


#' grep function changed
#' @keywords internal
#' @examples
#' \dontshow{
#' my_grep("x",c("xy","xz","ab"))
#' }
my_grep<-function(...) grep(fixed=TRUE,...)

