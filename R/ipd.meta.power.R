ipd.meta.power <- function(
               y0,
               y1,
               var0,
               var1,
               x0,
               x1,
               s20,
               s21,
               n0,
               n1,
               interaction,
               alpha=.05,
               logOR=TRUE,
               level=95)
{

#REPEATED CALLS TO IPD.META.POWER.BASE.FUNCTION

    types <- c("point","lower","upper")

    estimates <- lapply(types,function(x){
      ipd.meta.power.base.function(type=x,
	y0=y0,y1=y1,var0=var0,var1=var1,x0=x0,x1=x1,s20=s20,s21=s21,
	n0=n0,n1=n1,interaction=interaction,alpha=alpha,logOR=logOR,level=level)
 })

    SE <- estimates[[1]]$estimated.se
    SE.l <- estimates[[2]]$estimated.se
    SE.u <- estimates[[3]]$estimated.se

    power <- estimates[[1]]$estimated.power
    power.lower <- estimates[[3]]$estimated.power
    power.upper <- estimates[[2]]$estimated.power

return(
      list(
         estimated.power=power,
         power.lower=power.lower,
         power.upper=power.upper,
         estimated.se=SE,
         se.lower=SE.l,
         se.upper=SE.u,
         sigma=estimates[[1]]$sigma,
         sigma0=estimates[[1]]$sigma0,
         sigma1=estimates[[1]]$sigma1,
         level=level)
         )
}


ipd.meta.power.base.function <- function(
                             y0,
                             y1,
                             var0,
                             var1,
                             x0,
                             x1,
                             s20,
                             s21,
                             n0,
                             n1,
                             interaction,
                             alpha=.05,
                             logOR=TRUE,
                             level=95,
                             type=c("point","lower","upper"))
{

######DEPENDENT FUNCTIONS

power <- function(beta,se,alpha=.05){
      z <- qnorm(1-alpha/2)
      return(pnorm(-z+beta/se)+pnorm(-z-beta/se))
}

get.weight <- function(logOR=TRUE,z){

 if(logOR){

 theta.bar <- sum(exp(y0)/(1+exp(y0))*n0+exp(y1)/(1+exp(y1))*n1)/sum(n0+n1)
 w.se <- (1-2*theta.bar)*sqrt(theta.bar*(1-theta.bar)/sum(n0+n1))
 weight <- theta.bar*(1-theta.bar)+z*w.se
 return(1/weight)

  }
 else{

 return(1)
 }
}

construct.beta.info.matrix <- function(n0,n1,x0,x1,s20,s21,d,e,f,g){

#RETURNS THE X'SIGMA^(-1)X INFORMATION FOR BETA OF IPD STUDY
#GIVEN COMPONENTS OF WEIGHT MATRIX SIGMA^(-1) COVARIATE SUMMARIES
#d, e, f, g are vectors that are the length of the number of studies K
#n0 control study sample sizes
#n1 treatment study sample sizes
#x0 control study means for covariate
#x1 treatment study means for covariate
#s20 control sample variances for covariate
#s21 treatment sample variances

c.factors <- function(n0,n1,x0,x1,e,f,g){

#control = 0
#treatment = 1
#x are study means

c.ctrl <- e*n0*x0+f*n1*x1
c.trt <- f*n0*x0+g*n1*x1

return(list(c0=c.ctrl,c1=c.trt))
}


######C FACTORS FOR EACH COLUMN VECTOR

c0 <- c.factors(n0,n1,rep(1,length(n0)),rep(1,length(n1)),e,f,g)
c1 <- c.factors(n0,n1,rep(0,length(n0)),rep(1,length(n1)),e,f,g)
c2 <- c.factors(n0,n1,x0,x1,e,f,g)
c3 <- c.factors(n0,n1,rep(0,length(n0)),x1,e,f,g)

#FIRST ROW

X00 <- sum(n0*(d+c0$c0)+n1*(d+c0$c1))
X01 <- sum(n1*(d+c0$c1))
X02 <- sum(n0*x0*(d+c0$c0)+n1*x1*(d+c0$c1))
X03 <- sum(n1*x1*(d+c0$c1))

#SECOND ROW

X11 <- sum(n1*(d+c1$c1))
X12 <- sum(n1*x1*(d+c1$c1)+n0*c1$c0*x0)
X13 <- sum(n1*x1*(d+c1$c1))

#THIRD ROW

X22 <-
sum(d*(n0*(s20+x0^2)+n1*(s21+x1^2)-(s20+s21))+n0*c2$c0*x0+n1*c2$c1*x1)
X23 <- sum(d*(n1*(s21+x1^2)-(s21))+n1*c2$c1*x1)

#LAST COMPONENT
X33 <- sum(d*(n1*(s21+x1^2)-(s21))+n1*c3$c1*x1)

Sigma <- 

matrix(
c(X00,X01,X02,X03,
X01,X11,X12,X13,
X02,X12,X22,X23,
X03,X13,X23,X33),4,4
)
return(Sigma)
}

me.inverse <- function(a,b,c,n,m){

#Determine the study inverse of the marginal variance for a mixed effects model
#mixed effects for intercept and treatment group

#a+b is diagonal for control
#b is control diagonal (covariance)
#c is the off diagonal elements for treatment
#m treated subjects; n control subjects

block.matrix.inverse <- function(n,a,b){

#Returns the 

a.inverse <- 1/a
b.inverse <- -b/(a*(a+n*b))

return(c(a.inverse,off.diag=b.inverse,diag=a.inverse+b.inverse))
}

####

#E matrix (bottom left)

A.inverse.components <- block.matrix.inverse(n,a,b)

a.breve <- A.inverse.components[1]
b.breve <- A.inverse.components[2]

gamma <- b^2*n*(a.breve+n*b.breve)

E.inverse <- block.matrix.inverse(m,a,c-gamma)

#FE.inverse (off diagonals)

a.tilde <- E.inverse[1]
c.tilde <- E.inverse[2]

FE.inverse <- -b*(a.breve+n*b.breve)*(a.tilde+m*c.tilde)

#A.inverse (upper left)

gamma.breve <- b*m*(a.breve+n*b.breve)*(-FE.inverse)
diag=a.breve+b.breve+gamma.breve
off=b.breve+gamma.breve

d=a.breve
e=b.breve+gamma.breve
f=FE.inverse
g=E.inverse[2]

return(
c(d,e,f,g)
)
}

#########

 if(!logOR){
  #USE VARIANCES OF THE RESPONSE TO GET POOLED RESIDUAL VARIANCE EST
  #WHEN DIFFERENCE IN MEANS THIS IS THE SAMPLE VARIANCE BY GROUP
 sigma <- sum(var0*(n0-1)+var1*(n1-1))/sum(n0+n1-2)
  }
else{
 sigma <- 1
 }

fit0 <- rma(y0,var0/n0,mods=c(),method="REML")
fit1 <- rma(y1,var1/n1,mods=c(),method="REML")

 #OBTAINING BETWEEN-STUDY VARIANCE COMPONENTS

if(type=="point"){
 sigma0 <- fit0$tau2
 sigma1 <- fit1$tau2-sigma0
 sigma1 <- max(c(0,sigma1))
 weight <- get.weight(logOR=logOR,z=0)
 }
else if(type=="lower"){

sigma0 <- confint(fit0,level=level)$ci.lb[1]
sigma1 <- fit1$tau2-sigma0
se.sigma1 <- fit1$se.tau2[1]+fit0$se.tau2[1]
z.norm <- qnorm(1-level/200)
sigma1 <- max(c(0,sigma1+z.norm*se.sigma1))
 #Reciprocal is smaller for upper bound on binomial variance
weight <- get.weight(logOR=logOR,z=z.norm)
 }
else{

sigma0 <- confint(fit0,level=level)$ci.ub[1]
sigma1 <- fit1$tau2-sigma0

se.sigma1 <- fit1$se.tau2[1]+fit0$se.tau2[1]
z.norm <- qnorm(1-level/200)
sigma1 <- max(c(0,sigma1-z.norm*se.sigma1))

#Reciprocal is larger for lower bound on binomial variance
weight <- get.weight(logOR=logOR,z=-z.norm)
 }

inverse.components <- mapply(me.inverse,n=n0,m=n1,
MoreArgs=list(a=sigma,b=sigma0,c=(sigma0+sigma1)))

d <- inverse.components[1,]
e <- inverse.components[2,]
f <- inverse.components[3,]
g <- inverse.components[4,]

Info <- construct.beta.info.matrix(n0,n1,x0,x1,s20,s21,d,e,f,g)

SE <- sqrt(solve(Info)[4,4])*sqrt(weight)
power <- power(beta=interaction,se=SE,alpha=alpha)

return(
       list(
       estimated.power=power,
       estimated.se=SE,
       sigma=sigma,
       sigma0=sigma0,
       sigma1=sigma1,
       level=level)
       )
}

