coxmeta.fixed <- 

function(
ipd.formula,
meta.formula,
ipd.data,
meta.data,
sigma2,
study.group.interaction,
init.beta0,
init.beta,
error=1/1000)

{

    q = length(init.beta0)
    p = length(init.beta)
    beta0 = init.beta0
    beta = init.beta

    meta.formula <- as.formula(paste(c(meta.formula,
    		 as.character(ipd.formula)[3]),collapse="+"))

    #KNOWN COVARIANCE-VARIANCE FOR SURVIVAL OUTCOMES
    s <- meta.data[,all.vars(meta.formula)[[1]]]
    Sigma <- surv.covariance.variance(s,sigma2,study.group.interaction)

####NEWTON-RAPHSON MAXIMIZATION

    converge <- error*2
    loglik <- 0
    last <- FALSE

    monitor <- list(coef=c(),loglik=c())

    while(converge>error|!last){
    
    if(converge<error) last = TRUE

    #PATIENT LEVEL SCORE/INFO/LIKELIHOOD
    NR.cox <- coxph.with.offset(ipd.formula,ipd.data,beta=beta)        

    #STUDY LEVEL SCORE/INFO/LIKELIHOOD
    NR.surv <- surv.with.offset(meta.formula,meta.data,beta=c(beta0,beta),Sigma)

    #COMBINED
    NR <- combined.score.info(
       	  p=p,q=q,
	  U.cox=NR.cox$score,
	  I.cox=NR.cox$info,
	  loglik.cox=NR.cox$loglik,
	  U.surv=NR.surv$score,
	  I.surv=NR.surv$info,
	  loglik.surv=NR.surv$loglik)

    #UPDATE
    info <- matrix(NR$info,p+q,p+q)
    update <- as.vector(solve(info)%*%NR$score)
 

    beta0 <- beta0+update[1:q]
    beta <- beta+update[(q+1):(p+q)]
    converge = abs(NR$loglik-loglik)
    loglik = NR$loglik
  
    monitor$coef <- cbind(monitor$coef,c(beta0,beta))
    monitor$loglik <- c(monitor$loglik,loglik)
    }

return(list(
	coef=c(beta0,beta),
	var=solve(matrix(info,p+q,p+q)),
	loglik=loglik,
	monitor=monitor))
}

####DEPENDENT FUNCTIONS

combined.score.info <- 

function(p,q,U.cox,I.cox,loglik.cox,U.surv,I.surv,loglik.surv){

   U <- U.surv
   U[-(1:q)] <- U[-(1:q)]+U.cox	#REMOVE THE PARTS NOT SHARED

   m <- matrix(0,p+q,p+q)
   r <- row(m)
   c <- col(m)
   shared.index <- as.vector(r>q&c>q)

   I <- I.surv
   I[shared.index] <- I[shared.index]+I.cox

   loglik <- loglik.cox+loglik.surv

   return(list(score=U,info=I,loglik=loglik))
}

coxph.with.offset <- function(formula,data,beta,offset){

					#FORMULA CONSTRUCTION
   names <- all.vars(terms(formula))
   time <- data[,names[1]]
   event <- data[,names[2]]

   S <- Surv(time,event)
   covariates <- paste(attr(terms(formula),"term.labels"),collapse="+")

   f <- as.formula(paste(c("S~",covariates),c=""))

   if(!missing(offset)) f <- as.formula(paste(c(f,"offset(offset)"),
							collapse="+"))

  
					#MODEL ESTIMATION AT BETA/OFFSET
   fit <- coxph(f,data=data,init=beta,iter=0)
   fit.detail <- coxph.detail(fit)

					#SCORE/INFO/LOGLIK

   p <- length(beta)

   U <- fit.detail$score
   if(is.matrix(U)) U <- colSums(U) else U <- sum(U)

   I <- matrix(as.vector(fit.detail$imat),nrow=p^2)
   I <- apply(I,1,sum)   


   return(list(score=U,info=I,loglik=fit$loglik[2]))
}


surv.with.offset <- function(formula,data,beta,Sigma,offset){

   X <- model.matrix(terms(formula),data)

   S <- data[,all.vars(terms(formula))[1]]
   
   R <- log(-log(S))-X%*%beta
   
   if(!missing(offset)) R <- R-offset

   U <- as.vector(t(X)%*%solve(Sigma)%*%R)
   I <- as.vector(t(X)%*%solve(Sigma)%*%X)

   loglik <- -1/2*t(R)%*%solve(Sigma)%*%R

   return(list(score=U,info=I,loglik=loglik))
}

