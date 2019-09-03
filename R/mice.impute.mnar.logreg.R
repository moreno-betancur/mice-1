
mice.impute.mnar.logreg<-
  function (y, ry, x, wy = NULL, user,...) 
  {
    
    if (is.null(wy)) wy <- !ry
    
    wyold <- wy
    xold<-x
    
    ## Identifiable part: exactly the same as mice.impute.logreg
    
    # augment data in order to evade perfect prediction
    aug <- mice:::augment(y, ry, x, wy)
    x <- aug$x
    y <- aug$y
    ry <- aug$ry
    wy <- aug$wy
    w <- aug$w
    
    # fit model
    x <- cbind(1, as.matrix(x))
    expr <- expression(glm.fit(x = x[ry, , drop = FALSE], 
                               y = y[ry], 
                               family = quasibinomial(link = logit), 
                               weights = w[ry]))
    fit <- eval(expr)
    fit.sum <- summary.glm(fit)
    beta <- coef(fit)
    rv <- t(chol(mice:::sym(fit.sum$cov.unscaled)))
    beta.star <- beta + rv %*% rnorm(ncol(rv))
    
    
    ## Unidentifiable part  #e.g. user<-"3*intercept+2*bmi"
    
    mnar<-strsplit(unlist(strsplit(user,"+",fixed=T)),"*",fixed=T)
    mnar.parm<-as.numeric(unlist(lapply(mnar, function(x)x[1])))        #e.g. c("3","2")
    mnar.vars<-unlist(lapply(mnar, function(x)x[2]))        #e.g. c("intercept","bmi")
    xmnar<- cbind(1,as.matrix(xold[,mnar.vars[!mnar.vars=="intercept"]]))
    
    
    # draw imputations
    p <- 1/(1 + exp(-(x[wy, , drop = FALSE] %*% beta.star+as.matrix(xmnar[wyold,])%*%mnar.parm)))
    vec <- (runif(nrow(p)) <= p)
    vec[vec] <- 1
    if (is.factor(y)) {
      vec <- factor(vec, c(0, 1), levels(y))
    }
    
    return(vec)
  }
    
    

   