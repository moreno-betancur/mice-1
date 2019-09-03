#' # --------------------MICE.IMPUTE.NORM.MNAR----------------------------------
#' 
#' #'For NARFCS procedure only: Not-at-Random Imputation by Bayesian linear regression
#' #'
#' #'Imputes univariate missing data using Bayesian linear regression analysis in 
#' #'exactly the same way as in \code{mice.impute.norm} except the drawn imputations are 
#' #'then modified by adding a quantity as specified for the 
#' #'unidentified part of the model through argument passed to \code{blots} - see Details.
#' #'This results in "not-at-random" imputations.
#' #'
#' #'@aliases mice.impute.mnar.norm mnar.norm
#' #'@param y Vector to be imputed
#' #'@param ry Vector of missing data pattern (\code{FALSE}=missing,
#' #'\code{TRUE}=observed).
#' #'@param x Numeric design matrix with \code{length(y)} rows with predictors for 
#' #'\code{y}. Matrix \code{x} may have no missing values.
#' #'@param wy Logical vector of length \code{length(y)}. A \code{TRUE} value 
#' #'indicates locations in \code{y} for which imputations are created.
#' #'@param user Expression specifying the unidentifiable part of the model (MNAR delt-adjustment)
#' #'@param ... Other named arguments.
#' #'@return Vector with imputed data, same type as \code{y}, and of length 
#' #'\code{sum(wy)}
#' #'@author M. Moreno-Betancur, F. Leacy, D. Tompsett, I. White, 2017
#' #'@seealso \code{link{mice.impute.norm}}
#' #'@references 
#' #'
#' #'Moreno-Betancur M, Leacy FP, Tompsett D, White I. "mice: The NARFCS procedure for sensitivity analyses" 
#' #'(available at: \url{https://rawgit.com/moreno-betancur/NARFCS/master/README.html})
#' #'
#' #'Leacy FP. Multiple imputation under missing not at random assumptions via fully conditional 
#' #'specification 2016 (PhD thesis).
#' #'
#' #'Tompsett D, Leacy FP, Moreno-Betancur M, White I. On the use of the not at random fully conditional
#' #'specification procedure (NARFCS) in practice. Statistics in Medicine, 2018. 37(15):2338-2353.
#' #'@family univariate imputation functions
#' #'@export


mice.impute.mnar.norm <- function(y, ry, x, wy = NULL, user, ...) {
  
  xold<-x
  
  ## Identifiable part: exactly the same as mice.impute.norm
  if (is.null(wy)) wy <- !ry
  x <- cbind(1, as.matrix(x))
  parm <- .norm.draw(y, ry, x, ...)
  
  ## Unidentifiable part  #e.g. user<-"3*intercept+2*bmi"
  mnar<-strsplit(unlist(strsplit(user,"+",fixed=T)),"*",fixed=T)
  mnar.parm<-as.numeric(unlist(lapply(mnar, function(x)x[1])))        #e.g. c("3","2")
  mnar.vars<-unlist(lapply(mnar, function(x)x[2]))        #e.g. c("intercept","bmi")

  
  xmnar<- cbind(1,as.matrix(xold[,mnar.vars[!mnar.vars=="intercept"]]))

return(x[wy, ] %*% parm$beta + rnorm(sum(wy)) * parm$sigma+ as.matrix(xmnar[wy, ])%*%mnar.parm)
}

