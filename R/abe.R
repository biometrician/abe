

#' @import utils
if(getRversion() >= "2.15.1")  utils::globalVariables(c("value", "Variable", "VIF", "alpha.plot", "Model", "Frequency", "Var1", "Var2", "is"))



#' @title Augmented Backward Elimination
#'
#' @description  Function `abe` performs Augmented backward elimination where variable selection is based on the change-in-estimate and significance or information criteria.
#' It can also make a backward-selection based on significance or information criteria only by turning off the change-in-estimate criterion.
#'
#'
#' @param fit An object of class `"lm"`, `"glm"`, `"logistf"`, `"coxph"`, or `"survreg"` representing the fit.
#' Note, the functions should be fitted with argument `x=TRUE` and `y=TRUE` (or `model=TRUE` for `"logistf"` objects).
#' @param data data frame used when fitting the object `fit`.
#' @param include a vector containing the names of variables that will be included in the final model. These variables are used as only passive variables during modeling. These variables might be exposure variables of interest or known confounders.
#' They will never be dropped from the working model in the selection process,
#' but they will be used passively in evaluating change-in-estimate criteria of other variables.
#' Note, variables which are not specified as include or active in the model fit are assumed to be active and passive variables.
#' @param active a vector containing the names of active variables. These less important explanatory variables will only be used as active,
#' but not as passive variables when evaluating the change-in-estimate criterion.
#' @param tau  Value that specifies the threshold of the relative change-in-estimate criterion. Default is set to 0.05.
#' @param exact Logical, specifies if the method will use exact change-in-estimate or its approximation. Default is set to FALSE, which means that the method will use approximation proposed by Dunkler et al.
#' Note, setting to TRUE can severely slow down the algorithm, but setting to FALSE can in some cases lead to a poor approximation of the change-in-estimate criterion.
#' @param criterion String that specifies the strategy to select variables for the black list.
#' Currently supported options are significance level `'alpha'`, Akaike information criterion `'AIC'` and Bayesian information criterion `'BIC'`.
#' If you are using significance level, in that case you have to specify the value of 'alpha' (see parameter `alpha`) and type of the test statistic (see parameter `type.test`). Default is set to `"alpha"`.
#' @param alpha Value that specifies the level of significance as explained above. Default is set to 0.2.
#' @param type.test String that specifies which test should be performed in case the `criterion = "alpha"`.
#' Possible values are `"F"` and `"Chisq"` (default) for class `"lm"`, `"Rao"`, `"LRT"`, `"Chisq"` (default), `"F"` for class `"glm"` and `"Chisq"` for class `"coxph"`. See also [drop1()].
#' @param type.factor String that specifies how to treat factors, see details, possible values are `"factor"` and `"individual"`.
#' @param verbose Logical that specifies if the variable selection process should be printed. Note: this can severely slow down the algorithm. Default is set to TRUE.
#' @param ... Further arguments. Currently, this is primarily used to warn users about arguments that are no longer supported.
#' @details
#' Using the default settings ABE will perform augmented backward elimination based on significance.
#' The level of significance will be set to 0.2. All variables will be treated as "passive or active".
#' Approximated change-in-estimate will be used. Threshold of the relative change-in-estimate criterion will be 0.05.
#' Setting tau to a very large number (e.g. `Inf`) turns off the change-in-estimate criterion, and ABE will only perform backward elimination.
#' Specifying `"alpha" = 0` will include variables only because of the change-in-estimate criterion,
#' as then variables are not safe from exclusion because of their p-values.
#' Specifying `"alpha" = 1` will always include all variables.
#'
#' When using `type.factor="individual"` each dummy variable of a factor is treated as an individual explanatory variable, hence only this dummy variable can be removed from the model (warning: use sensible coding for the reference group).
#' Using `type.factor="factor"` will look at the significance of removing all dummy variables of the factor and can drop the entire variable from the model.
#'
#' In earlier versions, \code{abe} used to include an \code{exp.beta} argument. This is not supported anymore. Instead, the function now uses the exponential change in estimate for logistic, Cox, and parametric survival models only.
#'
#' @return An object of class `"lm"`, `"glm"`, `"coxph"`, or `"survreg"` representing the model chosen by abe method.
#' @references Daniela Dunkler, Max Plischke, Karen Lefondre, and Georg Heinze. Augmented backward elimination: a pragmatic and purposeful way to develop statistical models. PloS one, 9(11):e113677, 2014.
#'
#' @seealso [abe.resampling()], [lm()], [glm()] and [coxph()]
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
#' tau=0.05,exact=TRUE,criterion="alpha",alpha=0.2,
#' type.test="Chisq",verbose=TRUE)
#'
#' summary(abe.fit)
#'
#' # similar example, but turn off the change-in-estimate and perform
#' # only backward elimination
#'
#' abe.fit<-abe(fit,data=dd,include="x1",active="x2",
#' tau=Inf,exact=TRUE,criterion="alpha",alpha=0.2,
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
#' tau=0.05,exact=TRUE,criterion="alpha",alpha=0.2,
#' type.test="Chisq",verbose=TRUE,type.factor="factor")
#'
#' summary(abe.fit.fact)
#'
#' # treat each dummy of "x3" as a separate covariate:
#'
#' abe.fit.ind<-abe(fit,data=dd,include="x1",active="x2",
#' tau=0.05,exact=TRUE,criterion="alpha",alpha=0.2,
#' type.test="Chisq",verbose=TRUE,type.factor="individual")
#'
#' summary(abe.fit.ind)


abe<-function(fit,data=NULL,include=NULL,active=NULL,tau=0.05,exact=FALSE,criterion=c("alpha", "AIC", "BIC"),alpha=0.2,type.test=c("Chisq", "F", "Rao", "LRT"),type.factor=NULL,verbose=TRUE,...){
  if (is.null(data)) stop("Supply the data which were used when fitting the full model.")
assign(as.character(substitute(data)),data)

# match arguments
criterion <- match.arg(criterion)
type.test <- match.arg(type.test)

# check if user supplied the exp.beta argument and warn them if so
if("exp.beta" %in% names(list(...))) warning("Using exp.beta is not supported anymore. It is now automatically set to FALSE for linear models and TRUE for logistic and Cox models.")

# fix exp.beta depending on model type
exp.beta <- FALSE
if(class(fit)[1] == "glm" && fit$family$family=="binomial") exp.beta <- TRUE
if(class(fit)[1] == "coxph") exp.beta <- TRUE
if(inherits(fit, "logistf")) exp.beta <- TRUE
if(inherits(fit, "survreg")) exp.beta <- TRUE



# some checks and adjustments for logistf objects
if(inherits(fit, "logistf")){
  if(criterion != "alpha") stop("AIC and BIC selection are not supported for objects of class logistf")
  if(!("model" %in% names(fit))) stop("the model should be fitted with: model=TRUE")
  fit$x <- fit$model[, -1] # add x variables
  attr(fit$terms, "term.labels") <- if("(Intercept)" %in% fit$terms) fit$terms[-1] else fit$terms # add term.labels attribute
  attr(fit$terms, "dataClasses") <- sapply(fit$x, mode) # add dataClasses attribute
}

# check if model matrix is included in the fit
if (!"x"%in%names(fit)) stop("the model should be fitted with: x=T")
if (nrow(fit$x)!=nrow(data)) stop("Data object contains missing values. Remove all the missing values and refit the model.")

if (class(fit)[1]=="lm") if (!"y"%in%names(fit)) stop("the model should be fitted with: y=T")

if (!class(fit)[1]%in%c("lm","glm","coxph","brglmFit", "survreg", "logistf")) stop("this model is not supported")

if (sum(unlist(lapply(strsplit(colnames(fit$x),split=":"),function(x) length(x)!=1)))!=0) stop("interaction effects are not supported")

if (length(criterion)!=1) stop("you need to specify a single criterion")

if (colnames(fit$x)[1] == "(Intercept)") xm <- as.matrix(fit$x)[,-1] else xm <- as.matrix(fit$x)

if (!is.matrix(xm) ) stop("performing variable selection with a single variable in the model is meaningless")
if (sum(criterion%in%c("alpha","AIC","BIC"))==0) stop("valid criteria are alpha, AIC and BIC")
if (criterion=="alpha") if (alpha<0|alpha>1) stop("specify alpha between zero and one")
if (is.null(tau)) stop("Specify tau.")
if (tau<0) stop("Tau has to be >=0.")

nm.var<-ncol(xm)
if(class(fit)[1]=="lm"){
  n<-nrow(xm)
  epv<-n/nm.var
  if (epv<10) cat("Warning: Events per variable ratio is smaller than 10.")
}
if(class(fit)[1]=="glm"){
  if (fit$family$family=="binomial"){
  n<-min(table(fit$y))
  epv<-n/nm.var
  if (epv<10) cat("Warning: Events per variable ratio is smaller than 10.")
}
}
if(class(fit)[1]=="coxph"){
  n<-fit$nevent
  epv<-n/nm.var
  if (epv<10) cat("Warning: Events per variable ratio is smaller than 10.")
}

if(class(fit)[1]=="logistf"){
  n<-min(table(fit$y))
  epv<-n/nm.var
  if (epv<10) cat("Warning: Events per variable ratio is smaller than 10.")
}

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

if (length( my_grep("matrix",attributes(fit$terms)$dataClasses[-1]))!=0){ warning("The model contains a covariate which was fitted as matrix (splines etc.), treating it as a single covariate. \n Using type.factor=factor for factors. \n Using approximate change-in-estimate is inappropriate; using exact change in estimate instead.");bt<-abe.fact1(fit,data,include,active,tau,exp.beta,exact=TRUE,criterion,alpha,type.test,verbose)} else {

if (sum(attributes(fit$terms)$dataClasses[!my_grepl("strata",names(attributes(fit$terms)$dataClasses))]=="factor")>0)  {

if (is.null(type.factor)) {type.factor="factor"; warning("There are factors in the model but type.factor is not specified, using type.factor=factor")}

	if (type.factor=="factor") bt<-abe.fact1(fit,data,include,active,tau,exp.beta,exact,criterion,alpha,type.test,verbose) else bt<-abe.fact2(fit,data,include,active,tau,exp.beta,exact,criterion,alpha,type.test,verbose)


	} else  bt<-abe.num(fit,data,include,active,tau,exp.beta,exact,criterion,alpha,type.test,verbose)
}

eval(bt$call)
}







#' Resampled Augmented Backward Elimination
#'
#' Performs Augmented backward elimination on re-sampled datasets using different bootstrap and re-sampling techniques.
#'
#'
#' @param fit An object of class `"lm"`, `"glm"`, `"logistf"`, `"coxph"`, or `"survreg"` representing the fit.
#' Note, the functions should be fitted with argument `x=TRUE` and `y=TRUE` (or `model=TRUE` for `"logistf"` objects).
#' @param data data frame used when fitting the object `fit`.
#' @param include a vector containing the names of variables that will be included in the final model. These variables are used as passive variables during modeling. These variables might be exposure variables of interest or known confounders.
#' They will never be dropped from the working model in the selection process,
#' but they will be used passively in evaluating change-in-estimate criteria of other variables.
#' Note, variables which are not specified as include or active in the model fit are assumed to be active and passive variables.
#' @param active a vector containing the names of active variables. These less important explanatory variables will only be used as active,
#' but not as passive variables when evaluating the change-in-estimate criterion.
#' @param tau  Value that specifies the threshold of the relative change-in-estimate criterion. Default is set to 0.05.
#' @param exact Logical, specifies if the method will use exact change-in-estimate or approximated. Default is set to FALSE, which means that the method will use approximation proposed by Dunkler et al.
#' Note, setting to TRUE can severely slow down the algorithm, but setting to FALSE can in some cases lead to a poor approximation of the change-in-estimate criterion.
#' @param criterion String that specifies the strategy to select variables for the blacklist.
#' Currently supported options are significance level `'alpha'`, Akaike information criterion `'AIC'` and Bayesian information criterion `'BIC'`.
#' If you are using significance level, in that case you have to specify the value of 'alpha' (see parameter `alpha`). Default is set to `"alpha"`.
#' @param alpha Value that specifies the level of significance as explained above. Default is set to 0.2.
#' @param type.test String that specifies which test should be performed in case the `criterion = "alpha"`.
#' Possible values are `"F"` and `"Chisq"` (default) for class `"lm"`, `"Rao"`, `"LRT"`, `"Chisq"` (default), `"F"` for class `"glm"` and `"Chisq"` for class `"coxph"`. See also [drop1()].
#' @param type.factor String that specifies how to treat factors, see details, possible values are `"factor"` and `"individual"`.
#' @param num.resamples number of resamples.
#' @param type.resampling String that specifies the type of resampling. Possible values are `"Wallisch2021"`, `"bootstrap"`, `"mn.bootstrap"`, `"subsampling"`. Default is set to `"Wallisch2021"`. See details.
#' @param prop.sampling Sampling proportion. Only applicable for `type.boot="mn.bootstrap"` and `type.boot="subsampling"`, defaults to 0.5. See details.
#' @param save.out String that specifies if only the minimal output of the refitted models (`save.out="minimal"`) or the entire object (`save.out="complete"`) is to be saved. Defaults to `"minimal"`
#' @param parallel Logical, specifies if the calculations should be run in parallel `TRUE` or not `FALSE`. Defaults to `FALSE`. See details.
#' @param seed Numeric, a random seed to be used to form re-sampled datasets. Defaults to `NULL`. Can be used to assure complete reproducibility of the results, see Examples.
#' @param ... Further arguments. Currently, this is primarily used to warn users about arguments that are no longer supported.
#' @return an object of class `abe` for which `summary`, `plot` and `pie.abe` functions are available.
#' A list with the following elements:
#'
#' `coefficients` a matrix of coefficients of the final models obtained after performing ABE on re-sampled datasets; if using  `type.resampling="Wallisch2021"`, these models are obtained by using bootstrap.
#'
#' `coefficients.wallisch` if using `type.resampling="Wallisch2021"` the coefficients of the final models obtained after performing ABE using resampling with `prop.sampling` equal to 0.5; `NULL` when using any other option in `type.resampling`.
#'
#' `models` the final models obtained after performing ABE on re-sampled datasets, each object in the list is of the same class as `fit`; if using  `type.resampling="Wallisch2021"`, these models are obtained by using bootstrap. These are only returned if `save.out = "complete"`.
#'
#' `models.wallisch` similar as `models`; if using `type.resampling="Wallisch2021"` the coefficients and terms of the final models obtained after performing ABE using resampling with `prop.sampling` equal to 0.5; `NULL` when using any other option in `type.resampling`. These are only returned if `save.out = "complete"`.
#'
#' `model.parameters` a dataframe of alpha and tau values corresponding to the resampled models.
#'
#' `num.boot` number of resampled datasets
#'
#' `criterion` criterion used when constructing the black-list
#'
#' `all.vars` a list of variables used when estimating `fit`
#'
#' `fit.global` the initial model. In earlier versions of the package this parameter was called `fit.or`.
#'
#' `misc` the parameters of the call to `abe.resampling`
#'
#' `id` the rows of the data which were used when refitting the model; the list with elements `id1` (the rows used to refit the model; when `type.resampling="Wallisch2021"` these are based on bootstrap) and `id2` (`NULL` unless when `type.resampling="Wallisch2021"` in which case these are the rows used to refit the models based on subsampling)
#'
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @author Sladana Babic
#' @details `type.resampling` can be `bootstrap` (n observations drawn from the original data with replacement), `mn.bootstrap` (m out of n observations drawn from the original data with replacement), `subsampling` (m out of n observations drawn from the original data without replacement, where m is `prop.sampling*n` ) and `"Wallisch2021"`. When using `"Wallisch2021"` the resampling is done twice: first time using bootstrap (these results are contained in `models`) and the second time using resampling with `prop.sampling` equal to 0.5 (these results are contained in `models.wallisch`); see Walisch et al. (2021).
#' @details When using `parallel=TRUE` parallel backend must be registered before using `abe.resampling`. The parallel backends available will be system-specific; see [foreach()] for more details.
#' @details In earlier versions, \code{abe} used to include an \code{exp.beta} argument. This is not supported anymore. Instead, the function now uses the exponential change in estimate for logistic and Cox models only.
#' @references Daniela Dunkler, Max Plischke, Karen Lefondre, and Georg Heinze. Augmented backward elimination: a pragmatic and purposeful way to develop statistical models. PloS one, 9(11):e113677, 2014.
#' @references Riccardo De Bin, Silke Janitza, Willi Sauerbrei and Anne-Laure Boulesteix. Subsampling versus Bootstrapping in Resampling-Based Model Selection for Multivariable Regression. Biometrics 72, 272-280, 2016.
#' @references Wallisch C, Dunkler D, Rauch G, de Bin R, Heinze G. Selection of variables for multivariable models: Opportunities and limitations in quantifying model stability by resampling. Statistics in Medicine 40:369-381, 2021.
#' @seealso [abe()], [summary.abe()], [print.abe()], [plot.abe()], [pie.abe()]
#' @export
#' @import survival foreach
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
#' # use ABE on 50 re-samples considering different
#' # change-in-estimate thresholds and significance levels
#'
#' fit.resample<-abe.resampling(fit,data=dd,include="x1",active="x2",
#' tau=c(0.05,0.1),exact=TRUE,
#' criterion="alpha",alpha=c(0.2,0.05),type.test="Chisq",
#' num.resamples=50,type.resampling="Wallisch2021")
#'
#' summary(fit.resample)
#' print(fit.resample)
#'
#' # use ABE on 50 bootstrap re-samples considering different
#' # change-in-estimate thresholds and significance levels
#'
#' fit.resample<-abe.resampling(fit,data=dd,include="x1",active="x2",
#' tau=c(0.05,0.1),exact=TRUE,
#' criterion="alpha",alpha=c(0.2,0.05),type.test="Chisq",
#' num.resamples=50,type.resampling="bootstrap")
#'
#' summary(fit.resample)
#'
#' # use ABE on 50 subsamples randomly selecting 50% of subjects
#' # considering different change-in-estimate thresholds and
#' # significance levels
#'
#' fit.resample<-abe.resampling(fit,data=dd,include="x1",active="x2",
#' tau=c(0.05,0.1),exact=TRUE,
#' criterion="alpha",alpha=c(0.2,0.05),type.test="Chisq",
#' num.resamples=50,type.resampling="subsampling",prop.sampling=0.5)
#'
#' summary(fit.resample)
#'
#' #Assure reproducibility of the results
#'
#' fit.resample.1<-abe.resampling(fit,data=dd,include="x1",active="x2",
#' tau=c(0.05,0.1),exact=TRUE,
#' criterion="alpha",alpha=c(0.2,0.05),type.test="Chisq",
#' num.resamples=50,type.resampling="Wallisch2021")
#'
#' fit.resample.2<-abe.resampling(fit,data=dd,include="x1",active="x2",
#' tau=c(0.05,0.1),exact=TRUE,
#' criterion="alpha",alpha=c(0.2,0.05),type.test="Chisq",
#' num.resamples=50,type.resampling="Wallisch2021")
#'
#' #since different seed are used, fit.resample.1 and fit.resample.2 give different results
#'
#' fit.resample.1<-abe.resampling(fit,data=dd,include="x1",active="x2",
#' tau=c(0.05,0.1),exact=TRUE,
#' criterion="alpha",alpha=c(0.2,0.05),type.test="Chisq",
#' num.resamples=50,type.resampling="Wallisch2021",seed=87982)
#'
#' fit.resample.2<-abe.resampling(fit,data=dd,include="x1",active="x2",
#' tau=c(0.05,0.1),exact=TRUE,
#' criterion="alpha",alpha=c(0.2,0.05),type.test="Chisq",
#' num.resamples=50,type.resampling="Wallisch2021",seed=87982)
#'
#' #now fit.resample.1 and fit.resample.2 give exactly the same results
#'
#' #' Example to run parrallel computation on windows, using all but 2 cores
#'
#' #library(doParallel)
#' #N_CORES <- detectCores()
#' #cl <- makeCluster(N_CORES-2)
#' #registerDoParallel(cl)
#' #fit.resample<-abe.resampling(fit,data=dd,include="x1",active="x2",
#' #tau=c(0.05,0.1),exact=TRUE,
#' #criterion="alpha",alpha=c(0.2,0.05),type.test="Chisq",
#' #num.resamples=50,type.resampling="Wallisch2021")
#' #stopCluster(cl)


abe.resampling<-function(fit,data=NULL,include=NULL,active=NULL,tau=0.05,exact=FALSE,criterion=c("alpha", "AIC", "BIC"),alpha=0.2,type.test=c("Chisq", "F", "Rao", "LRT"),type.factor=NULL,num.resamples=100,type.resampling=c("Wallisch2021", "bootstrap", "mn.bootstrap", "subsampling"),prop.sampling=0.5,save.out=c("minimal", "complete"), parallel=FALSE,seed=NULL,...){
if (!is.null(seed)) set.seed(seed)


  # match arguments
  criterion <- match.arg(criterion)
  type.test <- match.arg(type.test)
  type.resampling <- match.arg(type.resampling)
  save.out <- match.arg(save.out)

  type.boot<-type.resampling
  num.boot<-num.resamples

  if (is.null(data)) stop("Supply the data which were used when fitting the full model.")

  # check if user supplied the exp.beta argument and warn them if so
  if("exp.beta" %in% names(list(...))) warning("Using exp.beta is not supported anymore. It is now automatically set to FALSE for linear models and TRUE for logistic and Cox models.")

  # fix exp.beta depending on model type
  exp.beta <- FALSE
  if(class(fit)[1] == "glm" && fit$family$family=="binomial") exp.beta <- TRUE
  if(class(fit)[1] == "coxph") exp.beta <- TRUE
  if(inherits(fit, "logistf")) exp.beta <- TRUE
  if(inherits(fit, "survreg")) exp.beta <- TRUE

  # some checks and adjustments for logistf objects
  if(inherits(fit, "logistf")){
    if(criterion != "alpha") stop("AIC and BIC selection are not supported for objects of class logistf")
    if(!("model" %in% names(fit))) stop("the model should be fitted with: model=TRUE")
    fit$x <- fit$model[, -1] # add x variables
    attr(fit$terms, "term.labels") <- if("(Intercept)" %in% fit$terms) fit$terms[-1] else fit$terms # add term.labels attribute
    attr(fit$terms, "dataClasses") <- sapply(fit$x, mode) # add dataClasses attribute
  }

  if (!"x"%in%names(fit)) stop("the model should be fitted with: x=T")
  if (nrow(fit$x)!=nrow(data)) stop("Data contains missing values. Remove all the missing values and refit the model.")

  if (class(fit)[1]=="lm") if (!"y"%in%names(fit)) stop("the model should be fitted with: y=T")

  if (!class(fit)[1]%in%c("lm","glm","coxph","brglmFit", "survreg", "logistf")) stop("this model is not supported")

  if (sum(unlist(lapply(strsplit(colnames(fit$x),split=":"),function(x) length(x)!=1)))!=0) stop("interaction effects are not supported")

  if (length(criterion)!=1) stop("you need to specify a single criterion")
  if (sum(criterion%in%c("alpha","AIC","BIC"))==0) stop("valid criteria are alpha, AIC and BIC")
  if (criterion=="alpha") if (sum(alpha<0)!=0|sum(alpha>1)!=0) stop("specify alphas between zero and one")
  if (is.null(tau)) stop("Specify tau.")
  if (sum(tau<0)!=0) stop("Taus has to be >=0.")

  if (length(type.boot)!=1) stop("you need to specify a single resampling method")

  if (criterion!="alpha") alpha=NULL

  if (colnames(fit$x)[1]=="(Intercept)") xm<-as.matrix(fit$x)[,-1] else xm<-as.matrix(fit$x)

  if (!is.matrix(xm) ) stop("performing variable selection with a single variable in the model is meaningless")


nm.var<-ncol(fit$x)
if (class(fit)[1]=="lm"){
  n<-nrow(model.matrix(fit))
  epv<-n/nm.var
  if (epv<10) cat("Warning: Events per variable ratio is smaller than 10.")
}
if (class(fit)[1]=="glm"){
  if (fit$family$family=="binomial"){
    n<-min(table(fit$y))
    epv<-n/nm.var
    if (epv<10) cat("Warning: Events per variable ratio is smaller than 10.")
  }
}
if (class(fit)[1]=="coxph"){
  n <-fit$nevent
  epv<-n/nm.var
  if (epv<10) cat("Warning: Events per variable ratio is smaller than 10.")

}


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


if (parallel==TRUE){

  if (type.boot!="Wallisch2021"){

    n<-nrow(fit$x)
    if (type.boot!="bootstrap") m<-round(n*prop.sampling,0) else m<-n
    if (criterion[1]=="BIC") k<-log(m)
    ids<-matrix(NA,ncol=m,nrow=num.boot)
    for (ii in 1:num.boot){
      if (type.boot=="bootstrap") idsi<-sample(1:n,n,replace=T)
      if (type.boot=="mn.bootstrap") idsi<-sample(1:n,m,replace=T)
      if (type.boot=="subsampling") idsi<-sample(1:n,m,replace=F)

      ids[ii,]<-idsi
    }

    type.boot.or<-type.boot
  #boot<-list()

  #i=0

  if (!is.null(alpha)&!is.null(tau)){

  boot<-foreach(a=alpha,.combine="c") %do% {#for (a in alpha){

    if (criterion[1]=="alpha") k<-qchisq(1-a,df=1)

  foreach(t=tau,.combine="c") %do% {#for (t in tau){
    #boot.i<-
    foreach(ii=1:num.boot,.packages=c("abe","survival")) %dopar% {#for (ii in 1:num.boot){
     #i=i+1


    data.boot<-data[ids[ii,],]

     fit.i<-my_update_boot(fit,data=data.boot)
    if (length( my_grep("matrix",attributes(fit$terms)$dataClasses[-1]))!=0){
      #boot[[i]]<-
        abe.fact1.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=a,type.test,k)
      } else {
    if (sum(attributes(fit$terms)$dataClasses[!my_grepl("strata",names(attributes(fit$terms)$dataClasses))]=="factor")>0)  {

    if (type.factor=="factor") {


                #boot[[i]]<-
                  abe.fact1.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=a,type.test,k) } else  {

                #boot[[i]]<-
                  abe.fact2.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=a,type.test,k)#,data=data.boot)


                      }


    } else  #boot[[i]]<-
        abe.num.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=a,type.test,k)

    }}
   }
  }}


  if (is.null(alpha)&!is.null(tau)){


     boot<-foreach(t=tau,.combine="c") %do% {# for (t in tau){
        #boot.i<-
       foreach(ii=1:num.boot,.packages=c("abe","survival") ) %dopar% {#for (ii in 1:num.boot){
          #i=i+1

          data.boot<-data[ids[ii,],]

            fit.i<-my_update_boot(fit,data=data.boot)

          if (length( my_grep("matrix",attributes(fit$terms)$dataClasses[-1]))!=0){
            #boot[[i]]<-
            abe.fact1.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=alpha,type.test,k)
            } else {

          if (sum(attributes(fit$terms)$dataClasses[!my_grepl("strata",names(attributes(fit$terms)$dataClasses))]=="factor")>0)  {

            if (type.factor=="factor") {


              #boot[[i]]<-
              abe.fact1.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=alpha,type.test,k) } else  {

               # boot[[i]]<-
                abe.fact2.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=alpha,type.test,k)#,data=data.boot)
                      }


          } else  #boot[[i]]<-
            abe.num.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=alpha,type.test,k)

        }}

    }}

  if (!is.null(alpha)&is.null(tau)){


    boot<-foreach(a=alpha,.combine="c") %do% {#for (a in alpha){

      if (criterion[1]=="alpha") k<-qchisq(1-a,df=1)

      #boot.i<-
      foreach(ii=1:num.boot ,.packages=c("abe","survival") ) %dopar% {#for (ii in 1:num.boot){
        #i=i+1

        data.boot<-data[ids[ii,],]

         fit.i<-my_update_boot(fit,data=data.boot)

        if (length( my_grep("matrix",attributes(fit$terms)$dataClasses[-1]))!=0){
          #boot[[i]]<-
          abe.fact1.boot(fit.i,data.boot,include,active,tau=tau,exp.beta,exact,criterion,alpha=a,type.test,k)
          } else {

        if (sum(attributes(fit$terms)$dataClasses[!my_grepl("strata",names(attributes(fit$terms)$dataClasses))]=="factor")>0)  {

          if (type.factor=="factor") {


            #boot[[i]]<-
            abe.fact1.boot(fit.i,data.boot,include,active,tau=tau,exp.beta,exact,criterion,alpha=a,type.test,k) } else  {

              #boot[[i]]<-
              abe.fact2.boot(fit.i,data.boot,include,active,tau=tau,exp.beta,exact,criterion,alpha=a,type.test,k)#,data=data.boot)

                       }


        } else  #boot[[i]]<-
          abe.num.boot(fit.i,data.boot,include,active,tau=tau,exp.beta,exact,criterion,alpha=a,type.test,k)

      }}

    }}

boot1<-boot
boot2<-NULL

  } else {

    type.boot.or<-type.boot
    type.boot<-"bootstrap"
    boot<-list()
    n<-nrow(fit$x)
    m<-n

    idsb<-matrix(NA,ncol=m,nrow=num.boot)
    for (ii in 1:num.boot){
      idsi<-sample(1:n,n,replace=T)
      idsb[ii,]<-idsi
    }

    if (criterion[1]=="BIC") k<-log(m)

    #i=0

    if (!is.null(alpha)&!is.null(tau)){

      boot<-foreach(a=alpha,.combine="c") %do% {#for (a in alpha){

        if (criterion[1]=="alpha") k<-qchisq(1-a,df=1)

        foreach(t=tau,.combine="c") %do% {#for (t in tau){
          foreach(ii=1:num.boot,.packages=c("abe","survival") ) %dopar% {#for (ii in 1:num.boot){
           # i=i+1
#source("R/abe.R") #used only to test is dopar works (since it looks for abe package on cran!)
#library(survival)

            data.boot<-data[idsb[ii,],]

            fit.i<-my_update_boot(fit,data=data.boot)
            if (length( my_grep("matrix",attributes(fit$terms)$dataClasses[-1]))!=0){
              #boot[[i]]<-
              abe.fact1.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=a,type.test,k)
              } else {
              if (sum(attributes(fit$terms)$dataClasses[!my_grepl("strata",names(attributes(fit$terms)$dataClasses))]=="factor")>0)  {

                if (type.factor=="factor") {


                  #boot[[i]]<-
                  abe.fact1.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=a,type.test,k) } else  {

                    #boot[[i]]<-
                    abe.fact2.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=a,type.test,k)#,data=data.boot)


                  }


              } else  #boot[[i]]<-
                abe.num.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=a,type.test,k)

            }}
        }
      }}


    if (is.null(alpha)&!is.null(tau)){


      boot<-foreach(t=tau,.combine="c") %do% {#for (t in tau){
        foreach(ii=1:num.boot,.packages=c("abe","survival") ) %dopar% {#for (ii in 1:num.boot){
          #i=i+1

          data.boot<-data[idsb[ii,],]

          fit.i<-my_update_boot(fit,data=data.boot)

          if (length( my_grep("matrix",attributes(fit$terms)$dataClasses[-1]))!=0){
            #boot[[i]]<-
            abe.fact1.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=alpha,type.test,k)
            } else {

            if (sum(attributes(fit$terms)$dataClasses[!my_grepl("strata",names(attributes(fit$terms)$dataClasses))]=="factor")>0)  {

              if (type.factor=="factor") {


                #boot[[i]]<-
                abe.fact1.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=alpha,type.test,k) } else  {

                  #boot[[i]]<-
                  abe.fact2.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=alpha,type.test,k)#,data=data.boot)
                }


            } else  #boot[[i]]<-
              abe.num.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=alpha,type.test,k)

          }}

      }}

    if (!is.null(alpha)&is.null(tau)){


      boot<-foreach(a=alpha,.combine="c") %do% {#for (a in alpha){

        if (criterion[1]=="alpha") k<-qchisq(1-a,df=1)

        foreach(ii=1:num.boot,.packages=c("abe","survival") ) %dopar% {#for (ii in 1:num.boot){
          #i=i+1

          data.boot<-data[idsb[ii,],]

          fit.i<-my_update_boot(fit,data=data.boot)

          if (length( my_grep("matrix",attributes(fit$terms)$dataClasses[-1]))!=0){
            #boot[[i]]<-
            abe.fact1.boot(fit.i,data.boot,include,active,tau=tau,exp.beta,exact,criterion,alpha=a,type.test,k)
            } else {

            if (sum(attributes(fit$terms)$dataClasses[!my_grepl("strata",names(attributes(fit$terms)$dataClasses))]=="factor")>0)  {

              if (type.factor=="factor") {


                #boot[[i]]<-
                abe.fact1.boot(fit.i,data.boot,include,active,tau=tau,exp.beta,exact,criterion,alpha=a,type.test,k) } else  {

                  #boot[[i]]<-
                  abe.fact2.boot(fit.i,data.boot,include,active,tau=tau,exp.beta,exact,criterion,alpha=a,type.test,k)#,data=data.boot)

                }


            } else  #boot[[i]]<-
              abe.num.boot(fit.i,data.boot,include,active,tau=tau,exp.beta,exact,criterion,alpha=a,type.test,k)

          }}

      }}

boot1<-boot
    type.boot<-"subsampling"
    n<-nrow(fit$x)
    m<-round(n*0.5,0)

    if (criterion[1]=="BIC") k<-log(m)

    idss<-matrix(NA,ncol=m,nrow=num.boot)
    for (ii in 1:num.boot){
      idsi<-sample(1:n,m,replace=F)
      idss[ii,]<-idsi
    }

    #boot<-list()

    #i=0

    if (!is.null(alpha)&!is.null(tau)){

     boot<-foreach(a=alpha,.combine="c") %do% {# for (a in alpha){

        if (criterion[1]=="alpha") k<-qchisq(1-a,df=1)

        foreach(t=tau,.combine="c") %do% {# for (t in tau){
          foreach(ii=1:num.boot ,.packages=c("abe","survival")) %dopar% {#for (ii in 1:num.boot){
           # i=i+1
#source("R/abe.R") #used only to test is dopar works (since it looks for abe package on cran!)
#library(survival)

            data.boot<-data[idss[ii,],]

            fit.i<-my_update_boot(fit,data=data.boot)
            if (length( my_grep("matrix",attributes(fit$terms)$dataClasses[-1]))!=0){
             # boot[[i]]<-
              abe.fact1.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=a,type.test,k)
              } else {
              if (sum(attributes(fit$terms)$dataClasses[!my_grepl("strata",names(attributes(fit$terms)$dataClasses))]=="factor")>0)  {

                if (type.factor=="factor") {


                  #boot[[i]]<-
                  abe.fact1.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=a,type.test,k) } else  {

                    #boot[[i]]<-
                    abe.fact2.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=a,type.test,k)#,data=data.boot)


                  }


              } else  #boot[[i]]<-
                abe.num.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=a,type.test,k)

            }}
        }
      }}


    if (is.null(alpha)&!is.null(tau)){


      boot<-foreach(t=tau,.combine="c") %do% {#for (t in tau){
        foreach(ii=1:num.boot,.packages=c("abe","survival") ) %dopar% {#for (ii in 1:num.boot){
          #i=i+1
          data.boot<-data[idss[ii,],]

          fit.i<-my_update_boot(fit,data=data.boot)

          if (length( my_grep("matrix",attributes(fit$terms)$dataClasses[-1]))!=0){
            #boot[[i]]<-
            abe.fact1.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=alpha,type.test,k)
            } else {

            if (sum(attributes(fit$terms)$dataClasses[!my_grepl("strata",names(attributes(fit$terms)$dataClasses))]=="factor")>0)  {

              if (type.factor=="factor") {


                #boot[[i]]<-
                abe.fact1.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=alpha,type.test,k) } else  {

                  #boot[[i]]<-
                  abe.fact2.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=alpha,type.test,k)#,data=data.boot)
                }


            } else  #boot[[i]]<-
              abe.num.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=alpha,type.test,k)

          }}

      }}

    if (!is.null(alpha)&is.null(tau)){


     boot<-foreach(a=alpha,.combine="c") %do% {# for (a in alpha){

        if (criterion[1]=="alpha") k<-qchisq(1-a,df=1)

        foreach(ii=1:num.boot ,.packages=c("abe","survival")) %dopar% {#for (ii in 1:num.boot){
          #i=i+1
           data.boot<-data[idss[ii,],]

          fit.i<-my_update_boot(fit,data=data.boot)

          if (length( my_grep("matrix",attributes(fit$terms)$dataClasses[-1]))!=0){
           # boot[[i]]<-
            abe.fact1.boot(fit.i,data.boot,include,active,tau=tau,exp.beta,exact,criterion,alpha=a,type.test,k)
            } else {

            if (sum(attributes(fit$terms)$dataClasses[!my_grepl("strata",names(attributes(fit$terms)$dataClasses))]=="factor")>0)  {

              if (type.factor=="factor") {


                #boot[[i]]<-
                abe.fact1.boot(fit.i,data.boot,include,active,tau=tau,exp.beta,exact,criterion,alpha=a,type.test,k) } else  {

                  #boot[[i]]<-
                  abe.fact2.boot(fit.i,data.boot,include,active,tau=tau,exp.beta,exact,criterion,alpha=a,type.test,k)#,data=data.boot)

                }


            } else  #boot[[i]]<-
              abe.num.boot(fit.i,data.boot,include,active,tau=tau,exp.beta,exact,criterion,alpha=a,type.test,k)

          }}

      }}

    boot2<-boot
      }


#  }

} else {


  if (type.boot!="Wallisch2021"){

    n<-nrow(fit$x)
    if (type.boot!="bootstrap") m<-round(n*prop.sampling,0) else m<-n
    if (criterion[1]=="BIC") k<-log(m)
    ids<-matrix(NA,ncol=m,nrow=num.boot)
    for (ii in 1:num.boot){
      if (type.boot=="bootstrap") idsi<-sample(1:n,n,replace=T)
      if (type.boot=="mn.bootstrap") idsi<-sample(1:n,m,replace=T)
      if (type.boot=="subsampling") idsi<-sample(1:n,m,replace=F)

      ids[ii,]<-idsi
    }

    type.boot.or<-type.boot
    boot<-list()

    i=0

    if (!is.null(alpha)&!is.null(tau)){

      for (a in alpha){

        if (criterion[1]=="alpha") k<-qchisq(1-a,df=1)

        for (t in tau){

          for (ii in 1:num.boot){
            i=i+1


            data.boot<-data[ids[ii,],]

            fit.i<-my_update_boot(fit,data=data.boot)
            if (length( my_grep("matrix",attributes(fit$terms)$dataClasses[-1]))!=0){
              boot[[i]]<-
              abe.fact1.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=a,type.test,k)
            } else {
              if (sum(attributes(fit$terms)$dataClasses[!my_grepl("strata",names(attributes(fit$terms)$dataClasses))]=="factor")>0)  {

                if (type.factor=="factor") {


                  boot[[i]]<-
                  abe.fact1.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=a,type.test,k) } else  {

                    boot[[i]]<-
                    abe.fact2.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=a,type.test,k)#,data=data.boot)


                  }


              } else  boot[[i]]<-
                abe.num.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=a,type.test,k)

            }}
        }
      }}


    if (is.null(alpha)&!is.null(tau)){


      for (t in tau){
        for (ii in 1:num.boot){
          i=i+1

          data.boot<-data[ids[ii,],]

          fit.i<-my_update_boot(fit,data=data.boot)

          if (length( my_grep("matrix",attributes(fit$terms)$dataClasses[-1]))!=0){
            boot[[i]]<-
            abe.fact1.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=alpha,type.test,k)
          } else {

            if (sum(attributes(fit$terms)$dataClasses[!my_grepl("strata",names(attributes(fit$terms)$dataClasses))]=="factor")>0)  {

              if (type.factor=="factor") {


                boot[[i]]<-
                abe.fact1.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=alpha,type.test,k) } else  {

                   boot[[i]]<-
                  abe.fact2.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=alpha,type.test,k)#,data=data.boot)
                }


            } else  boot[[i]]<-
              abe.num.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=alpha,type.test,k)

          }}

      }}

    if (!is.null(alpha)&is.null(tau)){


      for (a in alpha){

        if (criterion[1]=="alpha") k<-qchisq(1-a,df=1)

        for (ii in 1:num.boot){
          i=i+1

          data.boot<-data[ids[ii,],]

          fit.i<-my_update_boot(fit,data=data.boot)

          if (length( my_grep("matrix",attributes(fit$terms)$dataClasses[-1]))!=0){
            boot[[i]]<-
            abe.fact1.boot(fit.i,data.boot,include,active,tau=tau,exp.beta,exact,criterion,alpha=a,type.test,k)
          } else {

            if (sum(attributes(fit$terms)$dataClasses[!my_grepl("strata",names(attributes(fit$terms)$dataClasses))]=="factor")>0)  {

              if (type.factor=="factor") {


                boot[[i]]<-
                abe.fact1.boot(fit.i,data.boot,include,active,tau=tau,exp.beta,exact,criterion,alpha=a,type.test,k) } else  {

                  boot[[i]]<-
                  abe.fact2.boot(fit.i,data.boot,include,active,tau=tau,exp.beta,exact,criterion,alpha=a,type.test,k)#,data=data.boot)

                }


            } else  boot[[i]]<-
              abe.num.boot(fit.i,data.boot,include,active,tau=tau,exp.beta,exact,criterion,alpha=a,type.test,k)

          }}

      }}

    boot1<-boot
    boot2<-NULL

  } else {

    type.boot.or<-type.boot
    type.boot<-"bootstrap"
    boot<-list()
    n<-nrow(fit$x)
    m<-n

    idsb<-matrix(NA,ncol=m,nrow=num.boot)
    for (ii in 1:num.boot){
      idsi<-sample(1:n,n,replace=T)
      idsb[ii,]<-idsi
    }

    if (criterion[1]=="BIC") k<-log(m)

    i=0

    if (!is.null(alpha)&!is.null(tau)){

      for (a in alpha){

        if (criterion[1]=="alpha") k<-qchisq(1-a,df=1)

        for (t in tau){
          for (ii in 1:num.boot){
            i=i+1
            #source("R/abe.R") #used only to test is dopar works (since it looks for abe package on cran!)
            #library(survival)

            data.boot<-data[idsb[ii,],]

            fit.i<-my_update_boot(fit,data=data.boot)
            if (length( my_grep("matrix",attributes(fit$terms)$dataClasses[-1]))!=0){
              boot[[i]]<-
              abe.fact1.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=a,type.test,k)
            } else {
              if (sum(attributes(fit$terms)$dataClasses[!my_grepl("strata",names(attributes(fit$terms)$dataClasses))]=="factor")>0)  {

                if (type.factor=="factor") {


                  boot[[i]]<-
                  abe.fact1.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=a,type.test,k) } else  {

                    boot[[i]]<-
                    abe.fact2.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=a,type.test,k)#,data=data.boot)


                  }


              } else  boot[[i]]<-
                abe.num.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=a,type.test,k)

            }}
        }
      }}


    if (is.null(alpha)&!is.null(tau)){


      for (t in tau){
        for (ii in 1:num.boot){
          i=i+1

          data.boot<-data[idsb[ii,],]

          fit.i<-my_update_boot(fit,data=data.boot)

          if (length( my_grep("matrix",attributes(fit$terms)$dataClasses[-1]))!=0){
            boot[[i]]<-
            abe.fact1.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=alpha,type.test,k)
          } else {

            if (sum(attributes(fit$terms)$dataClasses[!my_grepl("strata",names(attributes(fit$terms)$dataClasses))]=="factor")>0)  {

              if (type.factor=="factor") {


                boot[[i]]<-
                abe.fact1.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=alpha,type.test,k) } else  {

                  boot[[i]]<-
                  abe.fact2.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=alpha,type.test,k)#,data=data.boot)
                }


            } else  boot[[i]]<-
              abe.num.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=alpha,type.test,k)

          }}

      }}

    if (!is.null(alpha)&is.null(tau)){


      for (a in alpha){

        if (criterion[1]=="alpha") k<-qchisq(1-a,df=1)

        for (ii in 1:num.boot){
          i=i+1

          data.boot<-data[idsb[ii,],]

          fit.i<-my_update_boot(fit,data=data.boot)

          if (length( my_grep("matrix",attributes(fit$terms)$dataClasses[-1]))!=0){
            boot[[i]]<-
            abe.fact1.boot(fit.i,data.boot,include,active,tau=tau,exp.beta,exact,criterion,alpha=a,type.test,k)
          } else {

            if (sum(attributes(fit$terms)$dataClasses[!my_grepl("strata",names(attributes(fit$terms)$dataClasses))]=="factor")>0)  {

              if (type.factor=="factor") {


                boot[[i]]<-
                abe.fact1.boot(fit.i,data.boot,include,active,tau=tau,exp.beta,exact,criterion,alpha=a,type.test,k) } else  {

                  boot[[i]]<-
                  abe.fact2.boot(fit.i,data.boot,include,active,tau=tau,exp.beta,exact,criterion,alpha=a,type.test,k)#,data=data.boot)

                }


            } else  boot[[i]]<-
              abe.num.boot(fit.i,data.boot,include,active,tau=tau,exp.beta,exact,criterion,alpha=a,type.test,k)

          }}

      }}

    boot1<-boot
    type.boot<-"subsampling"
    n<-nrow(fit$x)
    m<-round(n*0.5,0)

    if (criterion[1]=="BIC") k<-log(m)

    idss<-matrix(NA,ncol=m,nrow=num.boot)
    for (ii in 1:num.boot){
      idsi<-sample(1:n,m,replace=F)
      idss[ii,]<-idsi
    }

    boot<-list()

    i=0

    if (!is.null(alpha)&!is.null(tau)){

      for (a in alpha){

        if (criterion[1]=="alpha") k<-qchisq(1-a,df=1)

        for (t in tau){
          for (ii in 1:num.boot){
             i=i+1
            #source("R/abe.R") #used only to test is dopar works (since it looks for abe package on cran!)
            #library(survival)

            data.boot<-data[idss[ii,],]

            fit.i<-my_update_boot(fit,data=data.boot)
            if (length( my_grep("matrix",attributes(fit$terms)$dataClasses[-1]))!=0){
               boot[[i]]<-
              abe.fact1.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=a,type.test,k)
            } else {
              if (sum(attributes(fit$terms)$dataClasses[!my_grepl("strata",names(attributes(fit$terms)$dataClasses))]=="factor")>0)  {

                if (type.factor=="factor") {


                  boot[[i]]<-
                  abe.fact1.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=a,type.test,k) } else  {

                    boot[[i]]<-
                    abe.fact2.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=a,type.test,k)#,data=data.boot)


                  }


              } else  boot[[i]]<-
                abe.num.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=a,type.test,k)

            }}
        }
      }}


    if (is.null(alpha)&!is.null(tau)){


      for (t in tau){
        for (ii in 1:num.boot){
          i=i+1
          data.boot<-data[idss[ii,],]

          fit.i<-my_update_boot(fit,data=data.boot)

          if (length( my_grep("matrix",attributes(fit$terms)$dataClasses[-1]))!=0){
            boot[[i]]<-
            abe.fact1.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=alpha,type.test,k)
          } else {

            if (sum(attributes(fit$terms)$dataClasses[!my_grepl("strata",names(attributes(fit$terms)$dataClasses))]=="factor")>0)  {

              if (type.factor=="factor") {


                boot[[i]]<-
                abe.fact1.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=alpha,type.test,k) } else  {

                  boot[[i]]<-
                  abe.fact2.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=alpha,type.test,k)#,data=data.boot)
                }


            } else  boot[[i]]<-
              abe.num.boot(fit.i,data.boot,include,active,tau=t,exp.beta,exact,criterion,alpha=alpha,type.test,k)

          }}

      }}

    if (!is.null(alpha)&is.null(tau)){


       for (a in alpha){

        if (criterion[1]=="alpha") k<-qchisq(1-a,df=1)

        for (ii in 1:num.boot){
          i=i+1
          data.boot<-data[idss[ii,],]

          fit.i<-my_update_boot(fit,data=data.boot)

          if (length( my_grep("matrix",attributes(fit$terms)$dataClasses[-1]))!=0){
            boot[[i]]<-
            abe.fact1.boot(fit.i,data.boot,include,active,tau=tau,exp.beta,exact,criterion,alpha=a,type.test,k)
          } else {

            if (sum(attributes(fit$terms)$dataClasses[!my_grepl("strata",names(attributes(fit$terms)$dataClasses))]=="factor")>0)  {

              if (type.factor=="factor") {


                boot[[i]]<-
                abe.fact1.boot(fit.i,data.boot,include,active,tau=tau,exp.beta,exact,criterion,alpha=a,type.test,k) } else  {

                  boot[[i]]<-
                  abe.fact2.boot(fit.i,data.boot,include,active,tau=tau,exp.beta,exact,criterion,alpha=a,type.test,k)#,data=data.boot)

                }


            } else  boot[[i]]<-
              abe.num.boot(fit.i,data.boot,include,active,tau=tau,exp.beta,exact,criterion,alpha=a,type.test,k)

          }}

      }}

    boot2<-boot
  }



  } #end if parallel


  # add global fit
  fit.global <- fit


if (length( my_grep("matrix",attributes(fit$terms)$dataClasses[-1]))==0){
  if (sum(attributes(fit$terms)$dataClasses[!my_grepl("strata",names(attributes(fit$terms)$dataClasses))]=="factor")>0){
    if (type.factor=="individual") {
      df<-as.data.frame(fit$x)

      names(df)<-gsub("factor", replacement="", names(df), ignore.case = FALSE, perl = FALSE,
                      fixed = FALSE, useBytes = FALSE)

      names(df)<-gsub(")", replacement=".", names(df), ignore.case = FALSE, perl = FALSE,
                      fixed = FALSE, useBytes = FALSE)

      names(df)<-gsub("\\(", replacement="", names(df), ignore.case = FALSE, perl = FALSE,
                      fixed = FALSE, useBytes = FALSE)
      names(df)<-unlist(lapply(strsplit(names(df),split=""),function(x) {if(x[length(x)]==".") x<-x[-length(x)];paste(x,collapse="")}))
      names(df)<-unlist(lapply(strsplit(names(df),split="\\^"),paste,collapse=""))

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

        if ( colnames(fit$x)[1]=="(Intercept)" )   updt.f<-as.formula(paste("~",paste(check.names[-1],collapse="+"))) else  updt.f<-as.formula(paste("~",paste(check.names,collapse="+"),"-1"))

      df<-cbind(df,model.frame(fit),data)

      fit.global<-my_update(fit, updt.f   ,data=df)

      }
  }
}

  # add model parameters
  if(!is.null(alpha)){
    grid <- expand.grid(1:num.resamples, "tau" = tau, "alpha" = alpha)
    model.params <- grid[, c("alpha", "tau")]
  }
  if(is.null(alpha)){
    grid <- expand.grid(1:num.resamples, "tau" = tau)
    model.params <- grid[, c("tau"), drop = FALSE]
  }

  # add selected fit
  if(criterion != "alpha"){
    fit.selected <- lapply(unique(model.params$tau), function(tau.int){
      abe(fit, data = data, criterion = criterion, tau = tau.int, include = include,
          active = active, exact = exact, type.test = type.test, type.factor = type.factor, verbose = FALSE)
    })
  }
  if(criterion == "alpha"){
    fit.selected <- apply(unique(model.params), 1, function(x){
      abe(fit, data = data, criterion = criterion, alpha = x[1], tau = x[2], include = include,
          active = active, exact = exact, type.test = type.test, type.factor = type.factor, verbose = FALSE)
    })
  }
  names(fit.selected) <- 1:length(fit.selected)



  # create coefficient matrix
  coefficients <- matrix(0, nrow = length(boot1), ncol = length(fit.global$coefficients),
                         dimnames = list(1:length(boot1), names(fit.global$coefficients)))
  if(type.resampling=="Wallisch2021"){
    coefficients.wallisch <-  matrix(0, nrow = length(boot2), ncol = length(fit.global$coefficients),
                                     dimnames = list(1:length(boot2), names(fit.global$coefficients)))
  } else coefficients.wallisch <- NULL

  for (i in 1:nrow(coefficients)) {
    coefficients[i, names(boot1[[i]]$coefficients)] <- boot1[[i]]$coefficients
    if(type.resampling=="Wallisch2021") coefficients.wallisch[i, names(boot2[[i]]$coefficients)] <- boot2[[i]]$coefficients
  }


if(type.boot.or!="Wallisch2021") {id1<-ids;id2<-NULL} else {id1<-idsb;id2<-idss}

  misc<-list(tau=tau,criterion=criterion,alpha=alpha,type.boot=type.boot.or,prop.sampling=prop.sampling)
  id<-list(id.models=id1,id.wallisch=id2)

  if(save.out == "minimal"){
    res <- list(coefficients=coefficients,coefficients.wallisch=coefficients.wallisch, model.parameters = model.params,
                num.boot=num.boot,criterion=criterion,all.vars=names(coef(fit.global)),
                fit.global=fit.global, fit.selected = fit.selected,misc=misc,id=id,call=match.call())
  } else {
    res <- list(models=boot1,models.wallisch=boot2,coefficients=coefficients,
                coefficients.wallisch=coefficients.wallisch, model.parameters = model.params,
                num.boot=num.boot,criterion=criterion,all.vars=names(coef(fit.global)),
                fit.global=fit.global, fit.selected = fit.selected,misc=misc,id=id,call=match.call())
  }


  class(res)<-"abe"
  return(res)
}


#' Bootstrapped Augmented Backward Elimination
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function was deprecated, use `abe.resampling` instead.
#'
#' Performs Augmented backward elimination on re-sampled datasets using different bootstrap and re-sampling techniques.
#'
#'
#' @param fit An object of a class `"lm"`, `"glm"` or `"coxph"` representing the fit.
#' Note, the functions should be fitted with argument `x=TRUE` and `y=TRUE`.
#' @param data data frame used when fitting the object `fit`.
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
#' Currently supported options are significance level `'alpha'`, Akaike information criterion `'AIC'` and Bayesian information criterion `'BIC'`.
#' If you are using significance level, in that case you have to specify the value of 'alpha' (see parameter `alpha`). Default is set to `"alpha"`.
#' @param alpha Value that specifies the level of significance as explained above. Default is set to 0.2.
#' @param type.test String that specifies which test should be performed in case the `criterion = "alpha"`.
#' Possible values are `"F"` and `"Chisq"` (default) for class `"lm"`, `"Rao"`, `"LRT"`, `"Chisq"` (default), `"F"` for class `"glm"` and `"Chisq"` for class `"coxph"`. See also [drop1()].
#' @param type.factor String that specifies how to treat factors, see details, possible values are `"factor"` and `"individual"`.
#' @param num.boot number of bootstrap re-samples
#' @param type.boot String that specifies the type of bootstrap. Possible values are `"bootstrap"`, `"mn.bootstrap"`, `"subsampling"`,  see details
#' @param prop.sampling Sampling proportion. Only applicable for `type.boot="mn.bootstrap"` and `type.boot="subsampling"`, defaults to 0.5. See details.
#' @return an object of class `abe` for which `summary`, `plot` and `pie.abe` functions are available.
#' A list with the following elements:
#'
#' `models` the final models obtained after performing ABE on re-sampled datasets, each object in the list is of the same class as `fit`
#'
#' `alpha` the vector of significance levels used
#'
#' `tau` the vector of threshold values for the change-in-estimate
#'
#' `num.boot` number of re-sampled datasets
#'
#' `criterion` criterion used when constructing the black-list
#'
#' `all.vars` a list of variables used when estimating `fit`
#'
#' `fit.or` the initial model
#'
#' `misc` the parameters of the call to `abe.boot`
#'
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @author Sladana Babic
#' @details Used only for compatibility with the previous versions and will be removed at some point; see/use [abe.resampling()] instead.
#' @references Daniela Dunkler, Max Plischke, Karen Lefondre, and Georg Heinze. Augmented backward elimination: a pragmatic and purposeful way to develop statistical models. PloS one, 9(11):e113677, 2014.
#' @references Riccardo De Bin, Silke Janitza, Willi Sauerbrei and Anne-Laure Boulesteix. Subsampling versus Bootstrapping in Resampling-Based Model Selection for Multivariable Regression. Biometrics 72, 272-280, 2016.
#' @seealso [abe.resampling()]
#' @export
#' @import lifecycle
#' @keywords internal


abe.boot<-function(fit,data=NULL,include=NULL,active=NULL,tau=0.05,exp.beta=TRUE,exact=FALSE,criterion="alpha",alpha=0.2,type.test="Chisq",type.factor=NULL,num.boot=100,type.boot=c("bootstrap","mn.bootstrap","subsampling"),prop.sampling=0.5){

  # warn about deprecation
  lifecycle::deprecate_warn("5.1.1", "abe.boot()", "abe.resampling()")

  if (is.null(data)) stop("Supply the data which were used when fitting the full model.")

  if (!"x"%in%names(fit)) stop("the model should be fitted with: x=T")
  if (nrow(fit$x)!=nrow(data)) stop("Data contains missing values. Remove all the missing values and refit the model.")

  if (class(fit)[1]=="lm") if (!"y"%in%names(fit)) stop("the model should be fitted with: y=T")

  if (class(fit)[1]=="coxph"&exp.beta==F) stop("setting exp.beta=F for the cox model is not supported")

  if (!class(fit)[1]%in%c("lm","glm","coxph","brglmFit")) stop("this model is not supported")

  if (class(fit)[1]=="lm"&exp.beta==T&any(tau!=Inf)) warning("using change in estimate for exp(b) with linear model, try to use exp.beta=F")

  if (sum(unlist(lapply(strsplit(colnames(fit$x),split=":"),function(x) length(x)!=1)))!=0) stop("interaction effects are not supported")

  if (length(criterion)!=1) stop("you need to specify a single criterion")
  if (sum(criterion%in%c("alpha","AIC","BIC"))==0) stop("valid criteria are alpha, AIC and BIC")
  if (criterion=="alpha") if (sum(alpha<0)!=0|sum(alpha>1)!=0) stop("specify alphas between zero and one")
  if (is.null(tau)) stop("Specify tau.")
  if (sum(tau<0)!=0) stop("Taus has to be >=0.")

  if (length(type.boot)!=1) stop("you need to specify a single resampling method")
  nm.var<-ncol(fit$x)
  if (class(fit)[1]=="lm"){
    n<-nrow(model.matrix(fit))
    epv<-n/nm.var
    if (epv<10) cat("Warning: Events per variable ratio is smaller than 10.")
  }
  if (class(fit)[1]=="glm"){
    if (fit$family$family=="binomial"){
      n<-min(table(fit$y))
      epv<-n/nm.var
      if (epv<10) cat("Warning: Events per variable ratio is smaller than 10.")
    }
  }
  if (class(fit)[1]=="coxph"){
    n <-fit$nevent
    epv<-n/nm.var
    if (epv<10) cat("Warning: Events per variable ratio is smaller than 10.")

  }
  if (criterion!="alpha") alpha=NULL

  if (colnames(fit$x)[1]=="(Intercept)") xm<-as.matrix(fit$x)[,-1] else xm<-as.matrix(fit$x)

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
        df<-as.data.frame(fit$x)

        names(df)<-gsub("factor", replacement="", names(df), ignore.case = FALSE, perl = FALSE,
                        fixed = FALSE, useBytes = FALSE)

        names(df)<-gsub(")", replacement=".", names(df), ignore.case = FALSE, perl = FALSE,
                        fixed = FALSE, useBytes = FALSE)

        names(df)<-gsub("\\(", replacement="", names(df), ignore.case = FALSE, perl = FALSE,
                        fixed = FALSE, useBytes = FALSE)
        names(df)<-unlist(lapply(strsplit(names(df),split=""),function(x) {if(x[length(x)]==".") x<-x[-length(x)];paste(x,collapse="")}))
        names(df)<-unlist(lapply(strsplit(names(df),split="\\^"),paste,collapse=""))

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

        if ( colnames(fit$x)[1]=="(Intercept)" )   updt.f<-as.formula(paste("~",paste(check.names[-1],collapse="+"))) else  updt.f<-as.formula(paste("~",paste(check.names,collapse="+"),"-1"))

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




#' Summary Function
#'
#' makes a summary of a resampled version of ABE
#'
#' @param object an object of class `"abe"`, an object returned by a call to [abe.resampling()]
#' @param conf.level the confidence level, defaults to 0.95, see `details`
#' @param pval significance level to be used to determine a significant deviation from the expected pairwise inclusion frequency under independence.
#' @param alpha the alpha value for which the output is to be printed. If `NULL`, the output is printed for all alpha values.
#' @param tau the tau value for which the output is to be printed. If `NULL`, the output is printed for all tau values.
#' @param models.n controls the number of models printed for `model.rel.frequencies`. See details.
#' @param ... additional arguments affecting the summary produced.

#' @return a list with the following elements:
#'
#' `var.rel.frequencies`: inclusion relative frequencies for all variables from the initial model; if using `type.resampling="Wallisch2021"` in a call to [abe.resampling()] these results are based on subsampling with sampling proportion equal to 0.5, otherwise by using the method as specified by `type.sampling`
#'
#' `model.rel.frequencies`: relative frequencies of the final models; if using `type.resampling="Wallisch2021"` in a call to [abe.resampling()] these results are based on subsampling with sampling proportion equal to 0.5, otherwise by using the method as specified by `type.sampling`
#'
#' `var.coefs`: coefficient estimates and standard errors from the global and the selected model and medians, means, percentiles and standard deviations for the resampled estimates for each variable from the initial model; if using `type.resampling="Wallisch2021"` in a call to [abe.resampling()] these results are based on bootstrap, otherwise by using the method as specified by `type.sampling`
#'
#' `pair.rel.frequencies`: pairwise selection frequencies (in percent) for all pairs of variables. The significance of the deviation from the expected pairwise inclusion under independence is tested using a chi-squared test. If using `type.resampling="Wallisch2021"` in a call to [abe.resampling()] these results are based on subsampling with sampling proportion equal to 0.5, otherwise by using the method as specified by `type.sampling`
#'
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @author Sladana Babic
#' @author Gregor Steiner
#' @details Parameter `conf.level` defines the lower and upper quantile of the bootstrapped/resampled distribution such that equal proportion of values are smaller and larger than the lower and the upper quantile, respectively.
#' @details The `models.n` parameter controls the number of models printed in `model.rel.frequencies`. One option is to directly specify the number of models to return (i.e. an integer larger than 1). Alternatively, if `models.n` is set to a number less than (or equal to) 1, the number of models returned is such that the cumulative frequency attains that value. By default (`models.n = NULL`), the top 20 models or all models up to a cumulative frequency of 0.8, whichever is shorter, are returned. The selected model is marked with an asterisk. If it is not among the printed models, it is added as the last model.
#' @seealso [abe.resampling()], [print.abe()], [plot.abe()], [pie.abe()]
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
#' fit.resample<-abe.resampling(fit,data=dd,include="x1",active="x2",
#' tau=c(0.05,0.1),exact=TRUE,
#' criterion="alpha",alpha=c(0.2,0.05),type.test="Chisq",
#' num.resamples=50,type.resampling="Wallisch2021")
#'
#' summary(fit.resample)


summary.abe <- function(object, conf.level = 0.95, pval = 0.01, alpha = NULL, tau = NULL, models.n = NULL,...){

  # get different criterion combinations
  if(object$criterion != "alpha") names <- paste0(object$criterion, ", tau = ", unique(object$model.parameters[, "tau"]))
  if(object$criterion == "alpha") {
    uniq.comb <- unique(object$model.parameters)
    names <- paste0("alpha = ", uniq.comb[, "alpha"], ", ", "tau = ", uniq.comb[, "tau"])
  }

  ind <- 1:nrow(object$model.parameters) # index vector
  ind.split <- split(ind, ceiling(seq_along(ind) / object$num.boot)) # split index for different alpha/tau combinations (useful later)
  bool <- rep(TRUE, length(ind.split))

  # if alpha or tau is specified only compute for those values
  if(!is.null(alpha)) if(!all(alpha %in% object$misc$alpha)) stop("This value of alpha was not considered when using abe.resampling.")
  if(!is.null(tau)) if(!all(tau %in% object$misc$tau)) stop("This value of tau was not considered when using abe.resampling.")
  if(!is.null(alpha) | !is.null(tau)){
    if(object$criterion != "alpha"){
      if(is.null(tau)) tau <- unique(object$model.parameters[, "tau"])
      names.des <- paste0(object$criterion, ", tau = ", tau)
    }
    if(object$criterion == "alpha"){
      if(is.null(tau)) tau <- unique(object$model.parameters[, "tau"])
      if(is.null(alpha)) alpha <- unique(object$model.parameters[, "alpha"])
      names.des <- paste0("alpha = ", alpha, ", tau = ", tau)
    }

    bool <- names %in% names.des
    ind.split <- ind.split[bool]
    names <- names[bool]
  }


  # get relevant coef matrix depending on resampling type
  coef_matrix <- object$coefficients # if Wallisch, these are bootstrapped coefficients
  coef_matrix_sub <- object$coefficients # Here the same object is assigned, this object should be identical to coef.matrix if type.boot != "Wallisch2021". Else we need two different coefficient matrices
  if(object$misc$type.boot == "Wallisch2021") coef_matrix_sub <- object$coefficients.wallisch # if Wallisch, these are subsampled coefficients

  # Relative selection frequency for each variable
  var.rel.frequencies <- mapply(function(ind.models){
    # extract relevant models
    coef_matrix_int <- coef_matrix_sub[ind.models, ]

    # compute relative selection frequencies
    res <- colMeans(coef_matrix_int != 0)
    return(res)

  }, ind.split) |> t()

  rownames(var.rel.frequencies) <- names


  # Model selection frequencies
  model.frequency <- Map(function(ind.int, fit.sel){

    # get frequencies
    model.preds <- apply(coef_matrix_sub[ind.int, ], 1, function(x){
      # get predictors
      preds <- names(x[x != 0])
      preds <- preds[preds != "(Intercept)"] # exclude Intercept

      # paste the together and return string
      res <- paste(preds, collapse = " ")
      return(res)
    })

    # create output table
    res <- data.frame("Predictors" = unique(model.preds))

    # get counts
    res$Count <- sapply(res$Predictors, function(x) sum(x == model.preds))
    res$Percent <- res$Count / length(model.preds) * 100

    # sort by count
    res <- res[order(res$Count, decreasing = TRUE), ]
    rownames(res) <- 1:nrow(res)
    res$"Cumulative Percent" <- cumsum(res$Percent)

    # use default values if models.n is not specified
    if(is.null(models.n)){
      ind.cum.freq <- which((res$"Cumulative Percent" / 100) >= 0.8)[1]
      ind <- min(ind.cum.freq, 20) # return the shorter (=> min) of the two
    }

    # if models.n <= 1 use it as a cumulative frequency up until which to return
    if(!is.null(models.n) && models.n <= 1){
      ind.cum.freq <- which((res$"Cumulative Percent" / 100) >= models.n)[1]
      ind <- ind.cum.freq
    }

    # if models.n > 1 use it as the absolute number to return
    if(!is.null(models.n) && models.n > 1){
      ind <- min(models.n, nrow(res)) # if models.n is larger than the number of models return that number instead
    }

    # highlight the selected model
    terms.selected <- paste(names(fit.sel$coefficients)[names(fit.sel$coefficients) != "(Intercept)"], collapse = " ")
    ind.selected <- which(terms.selected == res$Predictors) # get index of selected model
    res[ind.selected, "Predictors"] <- paste0(res[ind.selected, "Predictors"], " *") # add asterisk for selected model

    # if selected model is not among the printed models add it as last model
    if(length(ind.selected) == 0){
      res <- res[1:ind, ]
      res <- rbind(res, c("Predictors" = paste0(terms.selected, " *"),
                          "Count" = 0, "Percent" = 0,
                          "Cumulative Percent" = res[ind, "Cumulative Percent"]))
    } else{
      if(ind.selected <= ind) res <- res[1:ind, ]
      if(ind.selected > ind) res <- res[c(1:ind, ind.selected), ]
    }

    # only return the specified number of models
    return(res)

  }, ind.split, object$fit.selected[bool]) |> setNames(names)


  # variable coefficient table
  var.coefs <- Map(function(ind.models, r.names, fit.sel){
    # extract relevant models
    coef_matrix_int <- coef_matrix[ind.models, ]

    # compute quantities of interest column wise
    res <- mapply(function(x, i){
      quantiles <- quantile(x, probs = c((1-conf.level)/2, conf.level+(1-conf.level)/2))
      vif <- mean(x != 0)

      return(c("Median" = median(x),
               "Lower quantile" = unname(quantiles[1]),
               "Upper quantile" = unname(quantiles[2]),
               "Mean" = mean(x),
               "SE" = sd(x),
               "RMSD Ratio" = unname(sqrt( sum((x - object$fit.global$coefficients[i])^2 / length(x)) / diag(vcov(object$fit.global))[i])),
               "Rel. cond. Bias" = unname(sum(x) / (object$fit.global$coefficients[i] * vif * length(x)) - 1)))

    }, as.data.frame(coef_matrix_int), 1:ncol(coef_matrix_int))

    # ad VIF
    res <- rbind("VIF" = var.rel.frequencies[r.names, ], res)

    # add resampling to names
    if(object$misc$type.boot %in% c("bootstrap", "mn.bootstrap")) resamp <- "Bootstrap"
    if(object$misc$type.boot == "subsampling") resamp <- "Subsampling"
    if(object$misc$type.boot == "Wallisch2021") resamp <- c("Subsampling", rep("Bootstrap", length(rownames(res)) - 1))
    rownames(res) <- paste(resamp, rownames(res))

    # add global and selected estimates
    coefs.selected <- rep(0, ncol(res)) |> setNames(colnames(res))
    se.selected <- rep(NA, ncol(res)) |> setNames(colnames(res))

    coefs.selected[names(fit.sel$coefficients)] <- fit.sel$coefficients
    se.selected[names(fit.sel$coefficients)] <- sqrt(diag(vcov(fit.sel)))

    mat <- rbind(object$fit.global$coefficients, sqrt(diag(vcov(object$fit.global))),
                 coefs.selected, se.selected)
    rownames(mat) <- c("Global estimate", "Global SE",
                       "Selected estimate", "Selected SE")

    res <- rbind(mat, res)

    # order by vif
    res <- res[, order(res[paste(resamp[1], "VIF"), ], decreasing = TRUE)]

    return(t(res))

  }, ind.split, names, object$fit.selected[bool]) |> setNames(names)



  # pairwise selection frequencies
  pair.rel.freqs <- Map(function(ind.models, r.names){

    pred <- colnames(var.rel.frequencies)
    pred_order <- order(var.rel.frequencies[r.names, ], decreasing = T)

    resampling_number <- object$num.boot

    coef_matrix_int <- coef_matrix_sub[ind.models, pred[pred_order]]
    coef_matrix_01 <- coef_matrix_int != 0

    resampling_VIF <- var.rel.frequencies[r.names, pred[pred_order]] * 100

    resampling_pairfreq <- matrix(100, ncol = length(pred), nrow = length(pred),
                                  dimnames = list(pred[pred_order],  pred[pred_order]))

    expect_pairfreq <- matrix(NA, ncol = length(pred), nrow = length(pred))

    # diagonal elements are individual selection frequencies
    diag(resampling_pairfreq) <- resampling_VIF

    for (i in 1:length(pred_order)){
      for (j in 1:length(pred_order)){
        # upper triangular matrix should be the pairwise selection frequencies
        if(j > i){
          # compute in how many of the models both variables appear together
          resampling_pairfreq[i, j] <- mean(apply(coef_matrix_01[, c(i, j)], 1, sum) == 2) * 100
          # compute expected selection frequency under independence (product of individual frequencies)
          expect_pairfreq[i, j] <- (resampling_VIF[i] / 100) * (resampling_VIF[j] / 100) * 100

          # lower triangular matrix should be the result of the chi-square test
          cont.table <- table(coef_matrix_01[, i], coef_matrix_01[, j]) # get contingency table

          if(ncol(cont.table) > 1 & nrow(cont.table) > 1){ # check if table has more than one row and column (i.e. dont test for variables with 100% selection frequency)
            sig <- suppressWarnings(chisq.test(cont.table)$p.value < pval) # perform chi-squared test
            text <- ifelse(!sig, "",
                           ifelse(as.numeric(resampling_pairfreq[i, j]) - as.numeric(expect_pairfreq[i, j]) > 0, "+", "-"))
          } else text <- ""


          resampling_pairfreq[j, i] <- text
        }
      }
    }

    return(resampling_pairfreq)

  }, ind.split, names) |> setNames(names)



  # return everything
  return(list(var.rel.frequencies = var.rel.frequencies,
              model.rel.frequencies = model.frequency,
              var.coefs = var.coefs,
              pair.rel.frequencies = pair.rel.freqs))

}


#' Print Function
#'
#' Prints a summary table of a bootstrapped/resampled version of ABE.
#' The table displays the relative inclusion frequencies of the covariates from the initial model,
#' the coefficient estimates and standard errors from the initial model (model with all covariates), the selected model,
#' resampled median and percentiles for the estimates of the regression coefficients for each variable from the initial model,
#' root mean squared difference ratio (RMSD) and relative bias conditional on selection (RBCS), see `details`.
#'
#' @param x an object of class `"abe"`, an object returned by a call to [abe.resampling()]
#' @param type the type of the output. `type = "coefficients"` prints summary statistics for each coefficient, `type = "coefficients reporting"` prints a reduced version of the coefficient statistics, and `type = "models"` reports model selection frequencies.
#' @param conf.level the confidence level, defaults to 0.95, see `details`
#' @param alpha the alpha value for which the output is to be printed, defaults to `NULL`
#' @param tau the tau value for which the output is to be printed, defaults to `NULL`
#' @param models.n controls the number of models printed if `type = "models"`. See details.
#' @param digits integer, indicating the number of digits to display in the table. Defaults to 2
#' @param ... additional arguments affecting the summary produced.
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @author Sladana Babic
#' @author Gregor Steiner
#' @details When using `type.resampling="Wallisch2021"` in a call to [abe.resampling()], the results for the relative inclusion frequencies of the covariates from the initial model are based on subsampling with sampling propotion equal to 0.5 and the other results are based on bootstrap as suggested by Wallisch et al. (2021); otherwise all the results are obtained by using the method as specified in `type.resampling`.
#' Parameter `conf.level` defines the lower and upper quantile of the bootstrapped/resampled distribution such that equal proportion of values are smaller and larger than the lower and the upper quantile, respectively.
#'
#' If `type = "models"`, the `models.n` parameter controls the number of models printed. One option is to directly specify the number of models to return (i.e. an integer larger than 1). Alternatively, if `models.n` is set to a number less than (or equal to) 1, the number of models returned is such that the cumulative frequency attains that value. By default (`models.n = NULL`), the top 20 models or all models up to a cumulative frequency of 0.8, whichever is shorter, are returned. The selected model is marked with an asterisk. If it is not among the printed models, it is added as the last model.
#' @references Wallisch C, Dunkler D, Rauch G, de Bin R, Heinze G. Selection of variables for multivariable models: Opportunities and limitations in quantifying model stability by resampling. Statistics in Medicine 40:369-381, 2021.
#' @seealso [abe.resampling()], [summary.abe()], [plot.abe()], [pie.abe()]
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
#' fit.resample<-abe.resampling(fit,data=dd,include="x1",active="x2",
#' tau=c(0.05,0.1),exact=TRUE,
#' criterion="alpha",alpha=c(0.2,0.05),type.test="Chisq",
#' num.resamples=50,type.resampling="Wallisch2021")
#'
#' print(fit.resample,conf.level=0.95,alpha=0.2,tau=0.05)


print.abe <- function(x, type = c("coefficients", "coefficients reporting", "models"), models.n = NULL, conf.level = 0.95, alpha = NULL, tau = NULL, digits = 3,...){

  # match arguments
  type <- match.arg(type)

  object <- x

  # coefficient table
  if(type == "coefficients"){

    sum.obj <- summary(object, conf.level = conf.level, alpha = alpha, tau = tau)

    res <- sum.obj$var.coefs
    # round
    res <- lapply(res, round, digits = digits)

    names(res) <- names(sum.obj$var.coefs)
    return(res)

  }

  if(type == "coefficients reporting"){

    sum.obj <- summary(object, conf.level = conf.level, alpha = alpha, tau = tau)

    # only keep essential columns
    res <- lapply(sum.obj$var.coefs, function(x) x[, c(1:5, 10)])
    # round
    res <- lapply(res, round, digits = digits)

    names(res) <- names(sum.obj$var.coefs)
    return(res)
  }

  # model selection frequencies
  if(type == "models"){

    res <- summary(object, conf.level = conf.level, alpha = alpha, tau = tau, models.n = models.n)$model.rel.frequencies
    return(res)

  }

}


#' Plot Function
#'
#' Plot function for the resampled/bootstrapped version of ABE.
#'
#' @param x an object of class `"abe"`, an object returned by a call to [abe.resampling()]
#' @param type.plot string which specifies the type of the plot. See details.
#' @param alpha values of alpha for which the plot is to be made (can be a vector of length >1)
#' @param tau values of tau for which the plot is to be made (can be a vector of length >1)
#' @param variable variables for which the plot is to be made (can be a vector of length >1)
#' @param type.stability string which specifies the type of stability plot. See details.
#' @param pval significance level to be used to determine a significant deviation from the expected pairwise inclusion frequency under independence (default 0.01). Only relevant if `type.plot="pairwise"`.
#' @param ... Arguments to be passed to methods, such as graphical parameters.
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @author Sladana Babic
#' @author Gregor Steiner
#' @details When using `type.plot="coefficients"` the function plots a histogram of the estimated regression coefficients for the specified variables, alpha(s) and tau(s) obtained from different re-sampled datasets.
#' When the variable is not included in the final model, its regression coefficient is set to zero. When using `type.resampling="Wallisch2021"` the plot is based on bootstrap, otherwise as specified in `type.resampling`.
#'
#' When using `type.plot="variables"` the function plots a barplot of the relative inclusion frequencies of the specified variables, for the specified values of alpha and tau. When using `type.resampling="Wallisch2021"` the plot is based on subsampling with sampling proportion equal to 0.5, otherwise as specified in `type.resampling`.
#'
#' When using `type.plot="models"` the function plots a barplot of the relative frequencies of the final models for specified alpha(s) and tau(s). When using `type.resampling="Wallisch2021"` the plot is based on subsampling with sampling proportion equal to 0.5, otherwise as specified in `type.resampling`.
#'
#' When using `type.plot="stability"` the function plots variable inclusion frequencies for each value of alpha. `type.stability` specifies if inclusion frequencies should be plotted as a function of alpha (default) or tau.
#'
#' When using `type.plot="pairwise"` the function plots a heatmap of differences between observed pairwise inclusion frequencies and the expected pairwise inclusion frequencies under independence. A high value indicates overselection, i.e. the pair of variables is selected together more often than expected under independence. Selection frequencies (in %) are displayed on top of the heatmap. See `summary.abe` for more details.
#' @import stats ggplot2 reshape2 tidytext
#' @export
#' @seealso [abe.resampling()], [summary.abe()], [pie.abe()]
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
#' fit.resample<-abe.resampling(fit,data=dd,include="x1",active="x2",
#' tau=c(0.05,0.1),exact=TRUE,
#' criterion="alpha",alpha=c(0.2,0.05),type.test="Chisq",
#' num.resamples=50,type.resampling="Wallisch2021")
#'
#' plot(fit.resample,type.plot="coefficients",
#' alpha=0.2,tau=0.1,variable=c("x1","x3"),
#' col="light blue")
#'
#' plot(fit.resample,type.plot="variables",
#' alpha=0.2,tau=0.1,variable=c("x1","x2","x3"),
#' col="light blue",horiz=TRUE,las=1)
#'
#' par(mar=c(4,6,4,2))
#' plot(fit.resample,type.plot="models",
#' alpha=0.2,tau=0.1,col="light blue",horiz=TRUE,las=1)
#'
#' fit.resample<-abe.resampling(fit,data=dd,include="x1",active="x2",
#' tau=c(0.05,0.1),exact=TRUE,
#' criterion="alpha",alpha=c(0.2,0.05),type.test="Chisq",
#' num.resamples=50,type.resampling="bootstrap")
#'
#' plot(fit.resample,type.plot="coefficients",
#' alpha=0.2,tau=0.1,variable=c("x1","x3"),
#' col="light blue")
#'
#' fit.resample<-abe.resampling(fit,data=dd,include="x1",active="x2",
#' tau=c(0.05,0.1),exact=TRUE,
#' criterion="alpha",alpha=c(0.2,0.05),type.test="Chisq",
#' num.resamples=50,type.resampling="subsampling")
#'
#' plot(fit.resample,type.plot="variables",
#' alpha=0.2,tau=0.1,variable=c("x1","x2","x3"),
#' col="light blue",horiz=TRUE,las=1)
#'
#' par(mar=c(4,6,4,2))
#' plot(fit.resample,type.plot="models",
#' alpha=0.2,tau=0.1,col="light blue",horiz=TRUE,las=1)

plot.abe<-function(x,type.plot=c("coefficients", "variables", "models", "stability", "pairwise"),alpha=NULL,tau=NULL,variable=NULL, type.stability = c("alpha", "tau"), pval = 0.01, ...){
  object<-x

  # match arguments
  type.plot <- match.arg(type.plot)
  type.stability <- match.arg(type.stability)

  if (type.plot=="coefficients"){

    # get coefficient matrix
    # if type.boot == "Wallisch2021", we want to use the bootstrapped coefficients here, so no need to use the coefficients wallisch
    coef_matrix <- object$coefficients

    # filter specified alpha and tau values if not NULL
    model.parameters <- object$model.parameters
    if(!is.null(alpha)){
      if(!all(alpha %in% object$misc$alpha)) stop("This value of alpha was not considered when using abe.resampling.")
      coef_matrix <- coef_matrix[model.parameters$alpha %in% alpha, ]
      model.parameters <- model.parameters[model.parameters$alpha %in% alpha, ]
    }
    if(!is.null(tau)){
      if(!all(tau %in% object$misc$tau)) stop("This value of tau was not considered when using abe.resampling.")
      coef_matrix <- coef_matrix[model.parameters$tau %in% tau, ]
      model.parameters <- model.parameters[model.parameters$tau %in% tau, ]
    }

    # filter for specified variables if not NULL
    if(!is.null(variable)){
      if(!(all(variable %in% colnames(coef_matrix)))) stop("At least one specified variable was not included in the initial model.")
      coef_matrix <- coef_matrix[, variable]
    }

    # get different criterion combinations
    if(object$criterion != "alpha"){
      coef.df <- data.frame(coef_matrix, "Model" = paste0(object$criterion, ", tau = ", model.parameters[, "tau"]))
    }
    if(object$criterion == "alpha") {
      coef.df <- data.frame(coef_matrix, "Model" = paste0("alpha = ", model.parameters[, "alpha"], ", tau = ", model.parameters[, "tau"]))
    }

    # reshape
    d.plot <- reshape2::melt(coef.df, id.vars = "Model")
    d.plot <- d.plot[!grepl("Intercept", d.plot$variable), ]

    # plot
    p <- ggplot(d.plot) +
      geom_histogram(aes(value), bins = 30) +
      geom_vline( xintercept = 0, color = "blue") +
      facet_wrap(~ Model + variable, scales = "free") +
      theme_bw() +
      xlab("Regression coefficient values") +
      ylab("Number of resamples") +
      theme(strip.text.x = element_text(size = 11))

  }

  if (type.plot=="variables"){

    sum.obj <- data.frame(summary(object, alpha = alpha, tau = tau)$var.rel.frequencies)

    if(!is.null(variable)){
      if(!(all(variable %in% colnames(sum.obj)))) stop("At least one specified variable was not included in the initial model.")
      sum.obj <- sum.obj[, variable]
    }

    d.plot <- reshape2::melt(data.frame(rownames(sum.obj), sum.obj), id.vars = 1)
    colnames(d.plot) <- c("Model", "Variable", "VIF")

    if(object$criterion == "AIC") d.plot$alpha.plot <- 0.157
    if(object$criterion == "BIC") d.plot$alpha.plot <- 1-pchisq(log(nrow(object$fit.global$x)), df=1)
    if(object$criterion == "alpha"){
      d.plot$alpha.plot <- as.numeric(sapply(strsplit(as.character(d.plot$Model), "alpha = |[,]"), "[[", 2))
    }

    d.plot$Model <- factor(d.plot$Model)
    d.plot$Variable <- tidytext::reorder_within(d.plot$Variable, d.plot$VIF, d.plot$Model)

    p <- ggplot(d.plot) +
      geom_col(aes(y = reorder(Variable, +VIF, max), x = VIF)) +
      geom_vline(data = d.plot, aes(xintercept = alpha.plot), col = 4) +
      facet_wrap( ~ Model, scales = "free") +
      tidytext::scale_y_reordered() +
      labs(y = NULL, x = "VIF") +
      theme_bw()

  }


  if (type.plot=="models"){
    if (!is.null(variable)) warning("Ploting relative frequencies of the final models but variable is not null. Ignoring the argument variable.")

    sum.obj<-summary(object, alpha = alpha, tau = tau)$model.rel.frequencies

    d.plot <- do.call(rbind, Map(function(x, y){
      d <- data.frame("Parameters" = y, "Model" = x$Predictors, "Frequency" = as.numeric(x$Percent) / 100)
      d
    }, sum.obj, names(sum.obj)))

    d.plot$Parameters <- factor(d.plot$Parameters)
    d.plot$Model <- tidytext::reorder_within(d.plot$Model, d.plot$Frequency, d.plot$Parameters)

    p <- ggplot(d.plot) +
      geom_col(aes(y = reorder(Model, +Frequency, max), x = Frequency)) +
      facet_wrap( ~ Parameters, scales = "free") +
      labs(y = NULL, x = "Frequency") +
      theme_bw() +
      tidytext::scale_y_reordered()

  }



  if(type.plot == "stability"){

    if(object$criterion == "alpha"){
      alphas <- sort(object$misc$alpha)
      if(!is.null(alpha)) alphas <- alpha
    }
    if(object$criterion == "AIC") alphas <- c("0.157")
    if(object$criterion == "BIC") alphas <- c(1-pchisq(log(nrow(object$fit.global$x)), df=1))
    taus <- sort(object$misc$tau)
    if(!is.null(tau)) taus <- tau

    if(length(alphas) > 1 & length(taus) == 1){
      if(type.stability == "tau") warning("type.stability = 'tau' requires more than 1 tau value, type.stability = 'alpha' is used instead.")
      type.stability <- "alpha"
    }
    if(length(taus) > 1 & length(alphas) == 1) type.stability <- "tau"

    sum.obj <- summary(object, alpha = alpha, tau = tau)
    var_rel_freqABE <- data.frame(sum.obj$var.rel.frequencies)[, -1]

    grid <- expand.grid( "tau" = taus, "alpha" = alphas)
    var_rel_freqABE <- cbind(var_rel_freqABE, grid)

    data_longABE <- reshape2::melt(var_rel_freqABE, id.vars = c("alpha", "tau"))

    if(type.stability == "alpha"){

      if(length(alphas) <= 1) stop("Stability plots require more than one alpha value.")

      p <- ggplot(data_longABE) +
        geom_line(aes(x = alpha, y = value, col = variable), linewidth = 0.75) +
        facet_wrap(~ paste0("Tau = ", tau)) +
        theme_bw() +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        labs(x = expression(alpha), y = "Inclusion frequencies", col = "") +
        ylim(0, 1)
    }

    if(type.stability == "tau"){

      if(length(taus) <= 1) stop("Stability plots require more than one tau value.")

      if(object$criterion != "alpha") data_longABE$alpha <- paste0(object$criterion, " (Alpha = ", data_longABE$alpha, ")")
      if(object$criterion == "alpha") data_longABE$alpha <- paste0("Alpha = ", data_longABE$alpha)

      p <- ggplot(data_longABE) +
        geom_line(aes(x = tau, y = value, col = variable), linewidth = 0.75) +
        facet_wrap(~ alpha) +
        scale_x_reverse() +
        theme_bw() +
        labs(x = expression(tau), y = "Inclusion Frequencies", col = "") +
        ylim(0, 1)
    }

  }


  if(type.plot == "pairwise"){

    sumobj <- summary(object, alpha = alpha, tau = tau, pval = pval)$pair.rel.frequencies

    d.plot <- do.call(rbind, Map(function(resampling_pairfreq, model){

      resampling_VIF <- as.numeric(diag(resampling_pairfreq))

      # matrix of actual selection frequencies
      m <- suppressWarnings(matrix(as.numeric(resampling_pairfreq),
                                   ncol = ncol(resampling_pairfreq),
                                   dimnames = dimnames(resampling_pairfreq)))
      diag(m) <- NA
      m[!is.na(m)] <- m[!is.na(m)]
      m[is.na(m)] <- 0
      m <- m + t(m)
      diag(m) <- resampling_VIF

      # compute expected selection frequencies under independence
      m.expect <- matrix(NA, ncol = ncol(m), nrow = nrow(m))
      for (i in 1:nrow(m)) {
        for (j in 1:ncol(m)) {
          m.expect[i, j] <- 100 * (m[i, i] / 100) * (m[j, j] / 100)
        }
      }

      # overselection is difference of actual and expected frequency
      m.diff <- m - m.expect

      # set diagonal elements to 0
      diag(m.diff) <- 0

      d.plot <- reshape2::melt(m.diff)
      d.plot$model <- model
      d.plot$text <- melt(t(resampling_pairfreq))$value

      return(d.plot)

    }, sumobj, names(sumobj)))

    d.plot$order <- 1:nrow(d.plot)
    d.plot$model <- factor(d.plot$model)
    d.plot$Var1 <- tidytext::reorder_within(d.plot$Var1, d.plot$order, d.plot$model)
    d.plot$Var2 <- tidytext::reorder_within(d.plot$Var2, d.plot$order, d.plot$model)


    p <- ggplot(d.plot, aes(x = Var1, y = ordered(factor(Var2), levels = rev(levels(factor(Var2)))))) +
      geom_tile(aes(fill = value)) +
      geom_tile(d.plot[d.plot$Var1 == d.plot$Var2, ],
                mapping = aes(x = Var1, y = ordered(factor(Var2), levels = rev(levels(factor(Var2))))),
                fill = NA, color = "black") + # frame diagonal elements
      geom_text(aes(label = text)) +
      facet_wrap(~ model, scales = "free") +
      scale_fill_gradient2(low = "#FFC20A", mid = "white", high = "#0C7BDC") +
      labs(x = "", y = "", fill = "Overselection") +
      theme_bw() +
      tidytext::scale_x_reordered() +
      tidytext::scale_y_reordered()


  }


  return(p)

}




#' Pie Function
#'
#' Pie function for the resampled/bootstrapped version of ABE. Plots a pie chart of the model frequencies for specified values of `alpha` and `tau`.
#'
#' @param x an object of class `"abe"`, an object returned by a call to [abe.resampling()]
#' @param alpha values of alpha for which the plot is to be made (can be a vector of length >1)
#' @param tau values of tau for which the plot is to be made (can be a vector of length >1)
#' @param labels plot labels, defaults to NA, i.e. no labels are ploted
#' @param ... Arguments to be passed to methods, such as graphical parameters (see [pie()], [barplot()], [hist()]).
#' @details When using `type.resampling="Wallisch2021"` the plot is based on subsampling with sampling proportion equal to 0.5, otherwise as specified in `type.resampling`.
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @author Sladana Babic
#' @import graphics
#' @export
#' @seealso [abe.resampling()], [summary.abe()], [plot.abe()]
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
#' fit.resample<-abe.resampling(fit,data=dd,include="x1",active="x2",
#' tau=c(0.05,0.1),exact=TRUE,
#' criterion="alpha",alpha=c(0.2,0.05),type.test="Chisq",
#' num.resamples=50,type.resampling="Wallisch2021")
#'
#' pie.abe(fit.resample, alpha=0.2,tau=0.1)
#'
#' fit.resample<-abe.resampling(fit,data=dd,include="x1",active="x2",
#' tau=c(0.05,0.1),exact=TRUE,
#' criterion="alpha",alpha=c(0.2,0.05),type.test="Chisq",
#' num.resamples=50,type.resampling="subsampling")
#'
#' pie.abe(fit.resample, alpha=0.2,tau=0.1)



pie.abe<-function(x,alpha=NULL,tau=NULL,labels=NA,...){

  # get model selection frequencies
  sum.obj <- summary(x, alpha = alpha, tau = tau)$model.rel.frequencies

  # plot
  par(mfcol = c(length(sum.obj), 1), mar = c(1, 1, 4, 1))
  invisible(
    Map(function(x.int, i){
      # drop zero entries
      x.int <- x.int[x.int$Percent > 0, ]

      pie(as.numeric(x.int$Percent), main = names(sum.obj)[i], labels=labels,...)
    }, sum.obj, 1:length(sum.obj))
  )

}





#' ABE for models which include only numeric covariates
#' @keywords internal
#' @examples
#' \dontrun{
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
#' tau=0.05,exact=TRUE,criterion="alpha",alpha=0.2,
#' type.test="Chisq",verbose=FALSE)
#' summary(abe.fit)
#' }

abe.num<-function(fit,data,include=NULL,active=NULL,tau=0.05,exp.beta=TRUE,exact=FALSE,criterion="alpha",alpha=0.2,type.test="Chisq",verbose=TRUE){


if (criterion[1]=="alpha") k<-qchisq(1-alpha,df=1)
if (criterion[1]=="AIC") k<-2

n<-nrow(fit$x)

if (criterion[1]=="BIC") k<-log(n)




 if (colnames(fit$x)[1]=="(Intercept)") xm<-as.matrix(fit$x)[,-1] else xm<-as.matrix(fit$x)


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


if ( sum(cnms%in%include)==length(cnms) ) stop("all variables are specified as passive, cannot perform variable selection")
if ( sum(cnms%in%active)==length(cnms) ) {
		include=active=NULL
		warning("all variables are specified as active, treating all variables as active or passive")
		}



if (is.null(include)&is.null(active)) {varnfix<- cnms;varpas<-cnms}
if (is.null(active)&!is.null(include)) {varnfix<- cnms[!cnms%in%include]; varpas<-cnms}
if (!is.null(active)&is.null(include))  {varnfix<-cnms;varpas<-cnms[!cnms%in%active]  }
if (!is.null(active)&!is.null(include)) {varnfix<-cnms[!cnms%in%include];varpas<-cnms[!cnms%in%active]}

varpas<-varpas[!my_grepl("strata",varpas)]

stop=F

while(stop==F){

# some necessary adjustments for logistf objects
if(inherits(fit, "logistf")){
  fit$x <- fit$model[, -1] # add x variables
  attr(fit$terms, "term.labels") <- if("(Intercept)" %in% fit$terms) fit$terms[-1] else fit$terms # add term.labels attribute
  attr(fit$terms, "dataClasses") <- sapply(fit$x, mode) # add dataClasses attribute
}

vcvm<-vcov(fit)

if (colnames(fit$x)[1]=="(Intercept)") vcvm<-vcvm[-1,-1]
if (is.null(dim(vcvm))) {
		vcvm<-matrix(vcvm,ncol=1,nrow=1)
		if (colnames(fit$x)[1]=="(Intercept)") colnames(vcvm)<-rownames(vcvm)<-colnames(vcov(fit))[-1] else colnames(vcvm)<-rownames(vcvm)<-colnames(vcov(fit))
		}


if (verbose==TRUE) {
	cat("\n\nModel under investigation:\n")
	print(fit$call)
	}

if(inherits(fit, "logistf")){
  scope <- varnfix # drop1.logistf() requires the scope as a vector of variables and not as a formula
} else {
  scope <- as.formula(paste("~",paste(varnfix,collapse=" + ") ))
}
if (criterion!="alpha") bl<-drop1(fit,scope=scope,k=k) else bl<-drop1(fit,scope=scope,test=type.test)

if (verbose==TRUE)  if (criterion!="alpha")  cat("Criterion for non-passive variables: "  ,paste( paste(varnfix,round(bl$AIC[-1],4),sep=" : "),collapse=" , "),"\n"  ,sep=""  ) else cat("Criterion for non-passive variables: "  ,paste( paste(varnfix,round(bl[-1,pmatch("Pr",names(bl))],4),sep=" : "),collapse=" , "),"\n"  ,sep=""  )


if (criterion!="alpha"){
  black.list.i<-varnfix[which(bl$AIC[-1]<bl$AIC[1])]
}
else {
  if(inherits(fit, "logistf")){
    ind <- pmatch("P",colnames(bl))
    black.list.i<-varnfix[which(bl[, ind]>alpha)]
  }  else {
    ind <- pmatch("Pr", names(bl))
    black.list.i<-varnfix[which(bl[-1, ind]>alpha)]
  }
}


if (length(black.list.i)!=0){

  if (criterion=="alpha") {
    if(inherits(fit, "logistf")){
       ind <- pmatch("P",colnames(bl))
       black.list.i <- varnfix[order(-bl[, ind])][1:length(black.list.i)]
       criterion.i <- bl[, ind][order(-bl[, ind])][1:length(black.list.i)]
    } else {
      ind <- pmatch("Pr", names(bl))
      black.list.i <- varnfix[order(-bl[-1, ind])][1:length(black.list.i)]
      criterion.i <- bl[-1, ind][order(-bl[-1, ind])][1:length(black.list.i)]
    }
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

				if (colnames(fit$x)[1]=="(Intercept)") {

					change.in.estimate<-abs(fit$coef[which(attributes(fit$terms)$term.labels[!my_grepl("strata",attributes(fit$terms)$term.labels)]%in%varpas[!varpas%in%black.list.i[i]])+1]-fit.i$coef[which(!attributes(fit.i$terms)$term.labels[!my_grepl("strata",attributes(fit.i$terms)$term.labels)]%in%active[!active%in%black.list.i[i]])+1])
					} else {
					change.in.estimate<-abs(fit$coef[which(attributes(fit$terms)$term.labels[!my_grepl("strata",attributes(fit$terms)$term.labels)]%in%varpas[!varpas%in%black.list.i[i]])]-fit.i$coef[which(!attributes(fit.i$terms)$term.labels[!my_grepl("strata",attributes(fit.i$terms)$term.labels)]%in%active[!active%in%black.list.i[i]])])
					}


				} else {
				if (colnames(fit$x)[1]=="(Intercept)") {
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
#' \dontrun{
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

  # some necessary adjustments for logistf objects
  if(inherits(fit, "logistf")){
    fit$x <- fit$model[, -1] # add x variables
    attr(fit$terms, "term.labels") <- if("(Intercept)" %in% fit$terms) fit$terms[-1] else fit$terms # add term.labels attribute
    attr(fit$terms, "dataClasses") <- sapply(fit$x, mode) # add dataClasses attribute
  }

  if (colnames(fit$x)[1]=="(Intercept)") xm<-as.matrix(fit$x)[,-1] else xm<-as.matrix(fit$x)


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

    # some necessary adjustments for logistf objects
    if(inherits(fit, "logistf")){
      if(!("model" %in% names(fit))) stop("the model should be fitted with: model=TRUE")
      fit$x <- fit$model[, -1] # add x variables
      attr(fit$terms, "term.labels") <- if("(Intercept)" %in% fit$terms) fit$terms[-1] else fit$terms # add term.labels attribute
      attr(fit$terms, "dataClasses") <- sapply(fit$x, mode) # add dataClasses attribute
    }


    vcvm<-vcov(fit)

    if (colnames(fit$x)[1]=="(Intercept)") vcvm<-vcvm[-1,-1]
    if (is.null(dim(vcvm))) {
      vcvm<-matrix(vcvm,ncol=1,nrow=1)
      if (colnames(fit$x)[1]=="(Intercept)") colnames(vcvm)<-rownames(vcvm)<-colnames(vcov(fit))[-1] else colnames(vcvm)<-rownames(vcvm)<-colnames(vcov(fit))
    }



      if(class(fit)[1] == "logistf"){
        scope <- varnfix # drop1.logistf() requires the scope as a vector of variables and not as a formula
      } else {
        scope <- as.formula(paste("~",paste(varnfix,collapse=" + ") ))
      }
      if (criterion!="alpha") bl<-drop1(fit,scope=scope,k=k) else bl<-drop1(fit,scope=scope,test=type.test)

    if (criterion!="alpha"){
      black.list.i<-varnfix[which(bl$AIC[-1]<bl$AIC[1])]
    }   else {
      if(inherits(fit, "logistf")){
        ind <- pmatch("P",colnames(bl))
        black.list.i<-varnfix[which(bl[, ind]>alpha)]
      }  else {
        ind <- pmatch("Pr", names(bl))
        black.list.i<-varnfix[which(bl[-1, ind]>alpha)]
      }
    }

    if (length(black.list.i)!=0){

      if (criterion=="alpha") {
        if(inherits(fit, "logistf")){
          ind <- pmatch("P",colnames(bl))
          black.list.i <- varnfix[order(-bl[, ind])][1:length(black.list.i)]
          criterion.i <- bl[, ind][order(-bl[, ind])][1:length(black.list.i)]
        } else {
          ind <- pmatch("Pr", names(bl))
          black.list.i <- varnfix[order(-bl[-1, ind])][1:length(black.list.i)]
          criterion.i <- bl[-1, ind][order(-bl[-1, ind])][1:length(black.list.i)]
        }
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

          if (colnames(fit$x)[1]=="(Intercept)") {

            change.in.estimate<-abs(fit$coef[which(attributes(fit$terms)$term.labels[!my_grepl("strata",attributes(fit$terms)$term.labels)]%in%varpas[!varpas%in%black.list.i[i]])+1]-fit.i$coef[which(!attributes(fit.i$terms)$term.labels[!my_grepl("strata",attributes(fit.i$terms)$term.labels)]%in%active[!active%in%black.list.i[i]])+1])
          } else {
            change.in.estimate<-abs(fit$coef[which(attributes(fit$terms)$term.labels[!my_grepl("strata",attributes(fit$terms)$term.labels)]%in%varpas[!varpas%in%black.list.i[i]])]-fit.i$coef[which(!attributes(fit.i$terms)$term.labels[!my_grepl("strata",attributes(fit.i$terms)$term.labels)]%in%active[!active%in%black.list.i[i]])])
          }


        } else {
          if (colnames(fit$x)[1]=="(Intercept)") {
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
#' \dontrun{
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



abe.fact1<-function(fit,data,include=NULL,active=NULL,tau=0.05,exp.beta=TRUE,exact=FALSE,criterion="alpha",alpha=0.2,type.test="Chisq",verbose=TRUE){



if (exact==F) {warning("there are factors in the model, using approximate change-in-estimate with this type.factor is inappropriate; using exact change in estimate instead"); exact=T}

# some necessary adjustments for logistf objects
if(inherits(fit, "logistf")){
  fit$x <- fit$model[, -1] # add x variables
  attr(fit$terms, "term.labels") <- if("(Intercept)" %in% fit$terms) fit$terms[-1] else fit$terms # add term.labels attribute
  attr(fit$terms, "dataClasses") <- sapply(fit$x, mode) # add dataClasses attribute
}


if (criterion[1]=="AIC") k<-2

n<-nrow(fit$x)

if (criterion[1]=="BIC") k<-log(n)



 if (colnames(fit$x)[1]=="(Intercept)") xm<-as.matrix(fit$x)[,-1] else xm<-as.matrix(fit$x)


var.mod<-attributes(fit$terms)$term.labels

if (length(var.mod)==0) fit$assign<-0 else {

fit$assign<-list()
fit$assign[1]<-1
name.cf<-names(coef(fit))
name.cfi<-name.cf
for(i in 1:length(name.cf)){
  name.cfi[i]<-paste(unlist(strsplit(paste(unlist(strsplit(name.cf[i],split="\\(")),collapse="\\("),split="\\)")),collapse = "\\)")

}



  for (i in 1:length(var.mod)){

    var.modi<-paste(unlist(strsplit(paste(unlist(strsplit(var.mod[i],split="\\(")),collapse="\\("),split="\\)")),collapse = "\\)")

    rep.var.names<-sum(grepl(var.modi,var.mod))-1
if (sum(rep.var.names)<0) rep.var.names<-0
    if (grepl("strata",var.mod[i])) fit$assign[[i+1]]<-0 else fit$assign[[i+1]]<-sum(my_grepl(var.modi,name.cfi))-rep.var.names

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

# some necessary adjustments for logistf objects
if(inherits(fit, "logistf")){
  if(!("model" %in% names(fit))) stop("the model should be fitted with: model=TRUE")
  fit$x <- fit$model[, -1] # add x variables
  attr(fit$terms, "term.labels") <- if("(Intercept)" %in% fit$terms) fit$terms[-1] else fit$terms # add term.labels attribute
  attr(fit$terms, "dataClasses") <- sapply(fit$x, mode) # add dataClasses attribute
}

vcvm<-vcov(fit)

if (colnames(fit$x)[1]=="(Intercept)") vcvm<-vcvm[-1,-1]


if (is.null(dim(vcvm))) {
		vcvm<-matrix(vcvm,ncol=1,nrow=1)
		colnames(vcvm)<-rownames(vcvm)<-rep(attributes(fit$terms)$term.labels,table(fit$assign[-1]))
		}


if (verbose==TRUE) {
	cat("\n\nModel under investigation:\n")
	print(fit$call)
	}

if(inherits(fit, "logistf")){
  scope <- varnfix # drop1.logistf() requires the scope as a vector of variables and not as a formula
} else {
  scope <- as.formula(paste("~",paste(varnfix,collapse=" + ") ))
}

if (criterion!="alpha") bl<-drop1(fit,scope=scope,k=k) else bl<-drop1(fit,scope=scope,test=type.test)
varnfixn<-unique(varnfix)
if (verbose==TRUE)  if (criterion!="alpha")  cat("Criterion for non-passive variables: "  ,paste( paste(varnfixn,round(bl$AIC[-1],4),sep=" : "),collapse=" , "),"\n"  ,sep=""  ) else cat("Criterion for non-passive variables: "  ,paste( paste(varnfixn,round(bl[-1,pmatch("Pr",names(bl))],4),sep=" : "),collapse=" , "),"\n"  ,sep=""  )


if (criterion!="alpha"){
  black.list.i<-varnfixn[which(bl$AIC[-1]<bl$AIC[1])]
}   else {
  if(inherits(fit, "logistf")){
    ind <- pmatch("P",colnames(bl))
    black.list.i<-varnfixn[which(bl[, ind]>alpha)]
  }  else {
    ind <- pmatch("Pr", names(bl))
    black.list.i<-varnfixn[which(bl[-1, ind]>alpha)]
  }
}


if (length(black.list.i)!=0){

  if (criterion=="alpha") {
    if(inherits(fit, "logistf")){
      ind <- pmatch("P",colnames(bl))
      black.list.i <- varnfixn[order(-bl[, ind])][1:length(black.list.i)]
      criterion.i <- bl[, ind][order(-bl[, ind])][1:length(black.list.i)]
    } else {
      ind <- pmatch("Pr", names(bl))
      black.list.i <- varnfixn[order(-bl[-1, ind])][1:length(black.list.i)]
      criterion.i <- bl[-1, ind][order(-bl[-1, ind])][1:length(black.list.i)]
    }
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
				  name.cfii<-name.cfi
				  for(iii in 1:length(name.cfi)){
				    name.cfii[iii]<-paste(unlist(strsplit(paste(unlist(strsplit(name.cfi[iii],split="\\(")),collapse="\\("),split="\\)")),collapse = "\\)")

				  }



				  for (iii in 1:length(var.mod.i)){


				  	    var.modii<-paste(unlist(strsplit(paste(unlist(strsplit(var.mod.i[iii],split="\\(")),collapse="\\("),split="\\)")),collapse = "\\)")

				    rep.var.namesi<-sum(grepl(var.modii,var.mod.i))-1
				    if (sum(rep.var.namesi)<0) rep.var.namesi<-0
				    if (grepl("strata",var.mod[iii])) fit.i$assign[[iii+1]]<-0 else fit.i$assign[[iii+1]]<-sum(my_grepl(var.modii,name.cfii))-rep.var.namesi
				   }

				  fit.i$assign<-unlist(fit.i$assign)

				  fit.i$assign<-rep(0:length(var.mod.i),fit.i$assign)
		}




				if (colnames(fit$x)[1]=="(Intercept)") {

					change.in.estimate<-abs(fit$coef[which(rep(attributes(fit$terms)$term.labels[!my_grepl("strata",var.mod)],table(fit$assign[-1]))%in%varpas[!varpas%in%black.list.i[i]])+1]-fit.i$coef[which(!rep(attributes(fit.i$terms)$term.labels[!my_grepl("strata",var.mod.i)],table(fit.i$assign[-1]))%in%active[!active%in%black.list.i[i]])+1])
					} else {
					change.in.estimate<-abs(fit$coef[which(rep(attributes(fit$terms)$term.labels[!my_grepl("strata",var.mod)],table(fit$assign[-1]))%in%varpas[!varpas%in%black.list.i[i]])]-
					                          fit.i$coef[which(!rep(attributes(fit.i$terms)$term.labels[!my_grepl("strata",var.mod.i)],table(fit.i$assign[-1]))%in%active[!active%in%black.list.i[i]])])
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
   			      name.cfi<-name.cf
   			      for(ii in 1:length(name.cf)){
   			        name.cfi[ii]<-paste(unlist(strsplit(paste(unlist(strsplit(name.cf[ii],split="\\(")),collapse="\\("),split="\\)")),collapse = "\\)")

   			      }



   			      for (ii in 1:length(var.mod)){

   			        var.modi<-paste(unlist(strsplit(paste(unlist(strsplit(var.mod[ii],split="\\(")),collapse="\\("),split="\\)")),collapse = "\\)")

   			        rep.var.names<-sum(grepl(var.modi,var.mod))-1
   			        if (sum(rep.var.names)<0) rep.var.names<-0
   			        if (grepl("strata",var.mod[ii])) fit$assign[[ii+1]]<-0 else fit$assign[[ii+1]]<-sum(my_grepl(var.modi,name.cfi))-rep.var.names
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
#' \dontrun{
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

  # some necessary adjustments for logistf objects
  if(inherits(fit, "logistf")){
    fit$x <- fit$model[, -1] # add x variables
    attr(fit$terms, "term.labels") <- if("(Intercept)" %in% fit$terms) fit$terms[-1] else fit$terms # add term.labels attribute
    attr(fit$terms, "dataClasses") <- sapply(fit$x, mode) # add dataClasses attribute
  }

  if (colnames(fit$x)[1]=="(Intercept)") xm<-as.matrix(fit$x)[,-1] else xm<-as.matrix(fit$x)


  var.mod<-attributes(fit$terms)$term.labels
  if (length(var.mod)==0) fit<-0 else {
    fit$assign<-list()
    fit$assign[1]<-1
    name.cf<-names(coef(fit))
    name.cfi<-name.cf
    for(i in 1:length(name.cf)){
      name.cfi[i]<-paste(unlist(strsplit(paste(unlist(strsplit(name.cf[i],split="\\(")),collapse="\\("),split="\\)")),collapse = "\\)")

    }



    for (i in 1:length(var.mod)){

      var.modi<-paste(unlist(strsplit(paste(unlist(strsplit(var.mod[i],split="\\(")),collapse="\\("),split="\\)")),collapse = "\\)")

      rep.var.names<-sum(grepl(var.modi,var.mod))-1
      if (sum(rep.var.names)<0) rep.var.names<-0
      if (grepl("strata",var.mod[i])) fit$assign[[i+1]]<-0 else fit$assign[[i+1]]<-sum(my_grepl(var.modi,name.cfi))-rep.var.names

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

    # some necessary adjustments for logistf objects
    if(inherits(fit, "logistf")){
      if(!("model" %in% names(fit))) stop("the model should be fitted with: model=TRUE")
      fit$x <- fit$model[, -1] # add x variables
      attr(fit$terms, "term.labels") <- if("(Intercept)" %in% fit$terms) fit$terms[-1] else fit$terms # add term.labels attribute
      attr(fit$terms, "dataClasses") <- sapply(fit$x, mode) # add dataClasses attribute
    }


    vcvm<-vcov(fit)

    if (colnames(fit$x)[1]=="(Intercept)") vcvm<-vcvm[-1,-1]



    if (is.null(dim(vcvm))) {
      vcvm<-matrix(vcvm,ncol=1,nrow=1)
      colnames(vcvm)<-rownames(vcvm)<-rep(attributes(fit$terms)$term.labels,table(fit$assign[-1]))
    }


    if(inherits(fit, "logistf")){
      scope <- varnfix # drop1.logistf() requires the scope as a vector of variables and not as a formula
    } else {
      scope <- as.formula(paste("~",paste(varnfix,collapse=" + ") ))
    }

    if (criterion!="alpha") bl<-drop1(fit,scope=scope,k=k) else bl<-drop1(fit,scope=scope,test=type.test)
    varnfixn<-unique(varnfix)


    if (criterion!="alpha"){
      black.list.i<-varnfixn[which(bl$AIC[-1]<bl$AIC[1])]
    }
    else {
      if(inherits(fit, "logistf")){
        ind <- pmatch("P",colnames(bl))
        black.list.i<-varnfixn[which(bl[, ind]>alpha)]
      }  else {
        ind <- pmatch("Pr", names(bl))
        black.list.i<-varnfixn[which(bl[-1, ind]>alpha)]
      }
    }

    if (length(black.list.i)!=0){

      if (criterion=="alpha") {
        if(inherits(fit, "logistf")){
          ind <- pmatch("P",colnames(bl))
          black.list.i <- varnfixn[order(-bl[, ind])][1:length(black.list.i)]
          criterion.i <- bl[, ind][order(-bl[, ind])][1:length(black.list.i)]
        } else {
          ind <- pmatch("Pr", names(bl))
          black.list.i <- varnfixn[order(-bl[-1, ind])][1:length(black.list.i)]
          criterion.i <- bl[-1, ind][order(-bl[-1, ind])][1:length(black.list.i)]
        }
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
          name.cfii<-name.cfi
          for(iii in 1:length(name.cfi)){
            name.cfii[iii]<-paste(unlist(strsplit(paste(unlist(strsplit(name.cfi[iii],split="\\(")),collapse="\\("),split="\\)")),collapse = "\\)")

          }



          for (iii in 1:length(var.mod.i)){


            var.modii<-paste(unlist(strsplit(paste(unlist(strsplit(var.mod.i[iii],split="\\(")),collapse="\\("),split="\\)")),collapse = "\\)")

            rep.var.namesi<-sum(grepl(var.modii,var.mod.i))-1
            if (sum(rep.var.namesi)<0) rep.var.namesi<-0
            if (grepl("strata",var.mod[iii])) fit.i$assign[[iii+1]]<-0 else fit.i$assign[[iii+1]]<-sum(my_grepl(var.modii,name.cfii))-rep.var.namesi
          }

          fit.i$assign<-unlist(fit.i$assign)

          fit.i$assign<-rep(0:length(var.mod.i),fit.i$assign)
        }


        if (colnames(fit$x)[1]=="(Intercept)") {

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
            name.cfi<-name.cf
            for(ii in 1:length(name.cf)){
              name.cfi[ii]<-paste(unlist(strsplit(paste(unlist(strsplit(name.cf[ii],split="\\(")),collapse="\\("),split="\\)")),collapse = "\\)")

            }



            for (ii in 1:length(var.mod)){

              var.modi<-paste(unlist(strsplit(paste(unlist(strsplit(var.mod[ii],split="\\(")),collapse="\\("),split="\\)")),collapse = "\\)")

              rep.var.names<-sum(grepl(var.modi,var.mod))-1
              if (sum(rep.var.names)<0) rep.var.names<-0
              if (grepl("strata",var.mod[ii])) fit$assign[[ii+1]]<-0 else fit$assign[[ii+1]]<-sum(my_grepl(var.modi,name.cfi))-rep.var.names
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
#' \dontrun{
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


abe.fact2<-function(fit,data,include=NULL,active=NULL,tau=0.05,exp.beta=TRUE,exact=FALSE,criterion="alpha",alpha=0.2,type.test="Chisq",verbose=TRUE){


# some necessary adjustments for logistf objects
if(inherits(fit, "logistf")){
  fit$x <- fit$model[, -1] # add x variables
  attr(fit$terms, "term.labels") <- if("(Intercept)" %in% fit$terms) fit$terms[-1] else fit$terms # add term.labels attribute
  attr(fit$terms, "dataClasses") <- sapply(fit$x, mode) # add dataClasses attribute
}

df<-as.data.frame(fit$x)

names(df)<-gsub("factor", replacement="", names(df), ignore.case = FALSE, perl = FALSE,
     fixed = FALSE, useBytes = FALSE)

names(df)<-gsub(")", replacement=".", names(df), ignore.case = FALSE, perl = FALSE,
     fixed = FALSE, useBytes = FALSE)

names(df)<-gsub("\\(", replacement="", names(df), ignore.case = FALSE, perl = FALSE,
     fixed = FALSE, useBytes = FALSE)


names(df)<-unlist(lapply(strsplit(names(df),split=""),function(x) {if(x[length(x)]==".") x<-x[-length(x)];paste(x,collapse="")}))
names(df)<-unlist(lapply(strsplit(names(df),split="\\^"),paste,collapse=""))
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



if ( colnames(fit$x)[1]=="(Intercept)" )   updt.f<-as.formula(paste("~",paste(check.names[-1],collapse="+"))) else  updt.f<-as.formula(paste("~",paste(check.names,collapse="+"),"-1"))
if ( class(fit)[1]=="coxph" )    updt.f<-as.formula(paste("~",paste(check.names,collapse="+")))

df<-cbind(df,model.frame(fit),data)

fit<-my_update(fit, updt.f   ,data=df)

if (sum(include%in%active)!=0) stop("at least one include variable is also specified as active")

 if (criterion[1]=="AIC") k<-2

n<-nrow(fit$x)

if (criterion[1]=="BIC") k<-log(n)


 if (colnames(fit$x)[1]=="(Intercept)") xm<-as.matrix(fit$x)[,-1] else xm<-as.matrix(fit$x)


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

# some necessary adjustments for logistf objects
if(inherits(fit, "logistf")){
  if(!("model" %in% names(fit))) stop("the model should be fitted with: model=TRUE")
  fit$x <- fit$model[, -1] # add x variables
  attr(fit$terms, "term.labels") <- if("(Intercept)" %in% fit$terms) fit$terms[-1] else fit$terms # add term.labels attribute
  attr(fit$terms, "dataClasses") <- sapply(fit$x, mode) # add dataClasses attribute
}


vcvm<-vcov(fit)

if (colnames(fit$x)[1]=="(Intercept)") vcvm<-vcvm[-1,-1]
if (is.null(dim(vcvm))) {
		vcvm<-matrix(vcvm,ncol=1,nrow=1)
		if (colnames(fit$x)[1]=="(Intercept)") colnames(vcvm)<-rownames(vcvm)<-colnames(vcov(fit))[-1] else colnames(vcvm)<-rownames(vcvm)<-colnames(vcov(fit))
		}


if (verbose==TRUE) {
	cat("\n\nModel under investigation:\n")
 	print(formula(fit))
	}


if(inherits(fit, "logistf")){
  scope <- varnfix # drop1.logistf() requires the scope as a vector of variables and not as a formula
} else {
  scope <- as.formula(paste("~",paste(varnfix,collapse=" + ") ))
}

if (criterion!="alpha") bl<-drop1(fit,scope=scope,k=k) else bl<-drop1(fit,scope=scope,test=type.test)

if (verbose==TRUE)  if (criterion!="alpha")  cat("Criterion for non-passive variables: "  ,paste( paste(varnfix,round(bl$AIC[-1],4),sep=" : "),collapse=" , "),"\n"  ,sep=""  ) else cat("Criterion for non-passive variables: "  ,paste( paste(varnfix,round(bl[-1,pmatch("Pr",names(bl))],4),sep=" : "),collapse=" , "),"\n"  ,sep=""  )

if (criterion!="alpha"){
  black.list.i<-varnfix[which(bl$AIC[-1]<bl$AIC[1])]
}
else {
  if(inherits(fit, "logistf")){
    ind <- pmatch("P",colnames(bl))
    black.list.i<-varnfix[which(bl[, ind]>alpha)]
  }  else {
    ind <- pmatch("Pr", names(bl))
    black.list.i<-varnfix[which(bl[-1, ind]>alpha)]
  }
}


if (length(black.list.i)!=0){

  if (criterion=="alpha") {
    if(inherits(fit, "logistf")){
      ind <- pmatch("P",colnames(bl))
      black.list.i <- varnfix[order(-bl[, ind])][1:length(black.list.i)]
      criterion.i <- bl[, ind][order(-bl[, ind])][1:length(black.list.i)]
    } else {
      ind <- pmatch("Pr", names(bl))
      black.list.i <- varnfix[order(-bl[-1, ind])][1:length(black.list.i)]
      criterion.i <- bl[-1, ind][order(-bl[-1, ind])][1:length(black.list.i)]
    }
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

				if (colnames(fit$x)[1]=="(Intercept)") {

					change.in.estimate<-abs(fit$coef[which(attributes(fit$terms)$term.labels[!my_grepl("strata",attributes(fit$terms)$term.labels)]%in%varpas[!varpas%in%black.list.i[i]])+1]-fit.i$coef[which(!attributes(fit.i$terms)$term.labels[!my_grepl("strata",attributes(fit.i$terms)$term.labels)]%in%active[!active%in%black.list.i[i]])+1])
					} else {
					change.in.estimate<-abs(fit$coef[which(attributes(fit$terms)$term.labels[!my_grepl("strata",attributes(fit$terms)$term.labels)]%in%varpas[!varpas%in%black.list.i[i]])]-fit.i$coef[which(!attributes(fit.i$terms)$term.labels[!my_grepl("strata",attributes(fit.i$terms)$term.labels)]%in%active[!active%in%black.list.i[i]])])
					}


				} else {
				if (colnames(fit$x)[1]=="(Intercept)") {
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
#' \dontrun{
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

  # some necessary adjustments for logistf objects
  if(inherits(fit, "logistf")){
    fit$x <- fit$model[, -1] # add x variables
    attr(fit$terms, "term.labels") <- if("(Intercept)" %in% fit$terms) fit$terms[-1] else fit$terms # add term.labels attribute
    attr(fit$terms, "dataClasses") <- sapply(fit$x, mode) # add dataClasses attribute
  }

  df<-as.data.frame(fit$x)

  names(df)<-gsub("factor", replacement="", names(df), ignore.case = FALSE, perl = FALSE,
                  fixed = FALSE, useBytes = FALSE)

  names(df)<-gsub(")", replacement=".", names(df), ignore.case = FALSE, perl = FALSE,
                  fixed = FALSE, useBytes = FALSE)

  names(df)<-gsub("\\(", replacement="", names(df), ignore.case = FALSE, perl = FALSE,
                  fixed = FALSE, useBytes = FALSE)

  names(df)<-unlist(lapply(strsplit(names(df),split=""),function(x) {if(x[length(x)]==".") x<-x[-length(x)];paste(x,collapse="")}))
  names(df)<-unlist(lapply(strsplit(names(df),split="\\^"),paste,collapse=""))

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

  if ( colnames(fit$x)[1]=="(Intercept)" )   updt.f<-as.formula(paste("~",paste(check.names[-1],collapse="+"))) else  updt.f<-as.formula(paste("~",paste(check.names,collapse="+"),"-1"))
  if ( class(fit)[1]=="coxph" )    updt.f<-as.formula(paste("~",paste(check.names,collapse="+")))

  df<-cbind(df,model.frame(fit),data)



  fit<-my_update(fit, updt.f   ,data=df)
  # some necessary adjustments for logistf objects
  if(inherits(fit, "logistf")){
    if(!("model" %in% names(fit))) stop("the model should be fitted with: model=TRUE")
    fit$x <- fit$model[, -1] # add x variables
    attr(fit$terms, "term.labels") <- if("(Intercept)" %in% fit$terms) fit$terms[-1] else fit$terms # add term.labels attribute
    attr(fit$terms, "dataClasses") <- sapply(fit$x, mode) # add dataClasses attribute
  }


  if (colnames(fit$x)[1]=="(Intercept)") xm<-as.matrix(fit$x)[,-1] else xm<-as.matrix(fit$x)


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

    # some necessary adjustments for logistf objects
    if(inherits(fit, "logistf")){
      if(!("model" %in% names(fit))) stop("the model should be fitted with: model=TRUE")
      fit$x <- fit$model[, -1] # add x variables
      attr(fit$terms, "term.labels") <- if("(Intercept)" %in% fit$terms) fit$terms[-1] else fit$terms # add term.labels attribute
      attr(fit$terms, "dataClasses") <- sapply(fit$x, mode) # add dataClasses attribute
    }

    vcvm<-vcov(fit)

    if (colnames(fit$x)[1]=="(Intercept)") vcvm<-vcvm[-1,-1]
    if (is.null(dim(vcvm))) {
      vcvm<-matrix(vcvm,ncol=1,nrow=1)
      if (colnames(fit$x)[1]=="(Intercept)") colnames(vcvm)<-rownames(vcvm)<-colnames(vcov(fit))[-1] else colnames(vcvm)<-rownames(vcvm)<-colnames(vcov(fit))
    }

    if(inherits(fit, "logistf")){
      scope <- varnfix # drop1.logistf() requires the scope as a vector of variables and not as a formula
    } else {
      scope <- as.formula(paste("~",paste(varnfix,collapse=" + ") ))
    }

    if (criterion!="alpha") bl<-drop1(fit,scope=scope,k=k) else bl<-drop1(fit,scope=scope,test=type.test)
    varnfixn<-unique(varnfix)


    if (criterion!="alpha"){
      black.list.i<-varnfix[which(bl$AIC[-1]<bl$AIC[1])]
    }
    else {
      if(inherits(fit, "logistf")){
        ind <- pmatch("P",colnames(bl))
        black.list.i<-varnfix[which(bl[, ind]>alpha)]
      }  else {
        ind <- pmatch("Pr", names(bl))
        black.list.i<-varnfix[which(bl[-1, ind]>alpha)]
      }
    }


    if (length(black.list.i)!=0){

      if (criterion=="alpha") {
        if(inherits(fit, "logistf")){
          ind <- pmatch("P",colnames(bl))
          black.list.i <- varnfix[order(-bl[, ind])][1:length(black.list.i)]
          criterion.i <- bl[, ind][order(-bl[, ind])][1:length(black.list.i)]
        } else {
          ind <- pmatch("Pr", names(bl))
          black.list.i <- varnfix[order(-bl[-1, ind])][1:length(black.list.i)]
          criterion.i <- bl[-1, ind][order(-bl[-1, ind])][1:length(black.list.i)]
        }
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

          if (colnames(fit$x)[1]=="(Intercept)") {

            change.in.estimate<-abs(fit$coef[which(attributes(fit$terms)$term.labels[!my_grepl("strata",attributes(fit$terms)$term.labels)]%in%varpas[!varpas%in%black.list.i[i]])+1]-fit.i$coef[which(!attributes(fit.i$terms)$term.labels[!my_grepl("strata",attributes(fit.i$terms)$term.labels)]%in%active[!active%in%black.list.i[i]])+1])
          } else {
            change.in.estimate<-abs(fit$coef[which(attributes(fit$terms)$term.labels[!my_grepl("strata",attributes(fit$terms)$term.labels)]%in%varpas[!varpas%in%black.list.i[i]])]-fit.i$coef[which(!attributes(fit.i$terms)$term.labels[!my_grepl("strata",attributes(fit.i$terms)$term.labels)]%in%active[!active%in%black.list.i[i]])])
          }


        } else {
          if (colnames(fit$x)[1]=="(Intercept)") {
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
#' \dontrun{
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

  #added to fix the issue with shrink, it would be probably be better to solve the issue with model.frame not working once we call upd2
  if (class(fit)=="coxph"){
  if (class(try(weights(fit),silent = TRUE))=="try-error") fit$weights<-rep(1L, fit$n) 
  }
  #end added
	
  fit
}

#' update function which searches for objects within the parent environment
#' @keywords internal
#' @examples
#' \dontrun{
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
#' \dontrun{
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
#' \dontrun{
#' my_grepl("x",c("xy","xz","ab"))
#' }
my_grepl<-function(...) grepl(fixed=TRUE,...)


#' grep function changed
#' @keywords internal
#' @examples
#' \dontrun{
#' my_grep("x",c("xy","xz","ab"))
#' }
my_grep<-function(...) grep(fixed=TRUE,...)

