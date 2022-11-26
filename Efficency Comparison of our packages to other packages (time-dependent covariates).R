library(fastcmprsk)
library(cmprskHD)
library(cmprsk)
library(crrp)
library(glmnet)
library(survival)
library(parallel)
library(rbenchmark)
library(MASS)
library(tictoc)
source("cmprskData_gen.R")
source("Pseudo_Entries_for_glmnet_timedependent.R")

## Test the running times for each package based on the same number of iterations, lambdas and same penalty

# Fix p=30 and use different n
p=30
# the number of different sample sizes
n_num=5
# For each sample size, the number of times you generate new data
# Please make sure loop is large enough (>=100) since we will use z-score to exclude some outliers of running time. Small loop will lead to bad z-score estimate and wrongly exclude some correct running times, which may cause NA in our final result.
loop=100



# Define n_num null time matrices
Time=list()
for(j in 1:n_num)
{
  Time[[j]]=matrix(nrow=loop,ncol=2)
}

# Define n_num null mean time matrices
MTime=list()
for(j in 1:n_num)
{
  MTime[[j]]=matrix(nrow=1,ncol=2)
}


# the increment of the sample size
incre=500
for(j in 1:n_num)
{
  n=1000+incre*(j-1)

  for(i in 1:loop)
  {
    # Generate Data
    # Baseline probability of event from cause 1
    prob=0.5
    # set true beta values
    # for similicity of simulation design, 
    # the relevant covariates are not changing with time
    non_zero_coeff 		= c(0.5,0.5)
    beta1 				= c(non_zero_coeff, rep(0, (p - length(non_zero_coeff))))	
    beta2=rep(c(-.5,.5,0),c(2,2,p-4))
    
    ## Generate data
    # Each id associates with multiple from-to entries for time-dependent covariates
    # The coding of survival entries is tricky, so the data generation steps are shown here.
    #----------------------------------------------------------------------------------------
    Z1=cov.gen.stdN(n,length(beta1)) # Covariate at baseline
    Z2=cov.gen.stdN(n,length(beta1)) # Covariate after the change
    Z2[,1:4] = Z1[,1:4] # keep the relevant covariate
    phaz=drop(exp(Z1%*%beta1)) # relative subdistribution risk cause 1
    phaz2 = drop(exp(Z1%*%beta2)) # cause 2 parameter
    
    cause=1+rbinom(n,1,prob=(1-prob)^phaz)
    
    # Generate the final sojourn times
    ftime=rep(0,n)
    ftime[cause==1]=FG.invT1(phaz[cause==1],prob)
    ftime[cause==2]=rT2.exp(phaz2[cause==2])
    
    Ctime=rcensor.unif(n) # Censoring time
    cause[Ctime<ftime] = 0 
    surv=Surv(pmin(ftime,Ctime),ftime<=Ctime)
    ftime = pmin(ftime,Ctime)
    cov.change = runif(n,0,1) # covariate change time
    id=1:length(surv)
    
    ## Generate weighted pseudo entries with left truncation time (start) and event time(stop) for each event for glmnet and phnet respectively
    # for glmnet
    dat<-FG.Surv_timedependent(surv,cause, Z1, Z2, cov.change, id)
    
    # for phnet
    # Most entries follow conventional coding
    cause01 = which(cause<2)
    n01 = length(cause01)
    end.before.change = (ftime <= cov.change)[cause01]
    fstatus1 = end.before.change*cause[cause01]
    surv1 = Surv(rep(0,n01),pmin(cov.change,ftime)[cause01],fstatus1)
    sec.entries = cause01[!end.before.change]
    surv2 = Surv(cov.change[sec.entries],ftime[sec.entries],cause[sec.entries])
    fstatus2 = cause[sec.entries]
    # Cause 2 events are different!
    cause2 = which(cause == 2)
    n2 = length(cause2)
    e2.after.change = (ftime >= cov.change)[cause2]
    entries.before.type2 = cause2[e2.after.change]
    n2.0 = length(entries.before.type2)
    surv3 = Surv(rep(0,n2.0),
                 cov.change[entries.before.type2],
                 rep(0,n2.0))
    fstatus3 = rep(0,n2.0)
    e2.before.change = cause2[!e2.after.change]
    n2.2 = n2-n2.0
    surv4 = Surv(c(rep(0,n2.2),cov.change[cause2]),
                 c(cov.change[e2.before.change], 
                   rep(max(ftime[cause==1])+1,n2)),
                 rep(0,n2.2+n2)) #typo here? without IPW weight?
    fstatus4 = rep(2,n2.2+n2) 
    # Put everything together
    id = c(cause01,sec.entries,entries.before.type2,
           e2.before.change,cause2)
    Y = rbind(surv1,surv2,surv3,surv4)
    Z = rbind(Z1[cause01,],Z2[sec.entries,],
              Z1[entries.before.type2,],Z1[e2.before.change,],
              Z2[cause2,])
    ftime = ftime[id]
    fstatus = c(fstatus1,fstatus2,fstatus3,fstatus4)
    
    
    # cmprskHD
    # generate the largest lambda that makes all the betas shrink to 0
    fit1<-phnet(Z,Y,ftime=ftime,fstatus=fstatus,id=id,family="finegray",alpha=1,nlambda=100,lambda.min.ratio=0.0001,standardize=TRUE,maxit=10000)
    
    # generate the lambda track based on the largest lambda and the smallest lambda set to be 0.00001
    lambda_seq<-seq(max(fit1$lambda),(0.1)^4,-(max(fit1$lambda)-(0.1)^4)/100)
    lambda_seq<-rev(lambda_seq)
    tic.clearlog()
    tic()
    fit1<-phnet(Z,Y, ftime = ftime, fstatus = fstatus,id=id,family="finegray",alpha=1,lambda=lambda_seq,standardize=TRUE,maxit=10000)
    toc(log = TRUE, quiet = TRUE)
    a<-tic.log(FALSE)
    
    # store the running time for loop i
    Time1 <- unlist(lapply(a, function(x) x$toc - x$tic))
    Time[[j]][i,1]=Time1
    
    ## glmnet
    # Time1<-benchmark(glmnet=glmnet_FG(cov,Y,cause),replications=1)
    tic.clearlog()
    tic()
    fit2<-glmnet(dat$Z,dat$surv,family="cox",weights=dat$weight,alpha=1,lambda=lambda_seq,standardize=TRUE,maxit=10000)
    toc(log = TRUE, quiet = TRUE)
    a<-tic.log(FALSE)
    # store the running time for loop i
    Time1 <- unlist(lapply(a, function(x) x$toc - x$tic))
    Time[[j]][i,2]=Time1
    
  }
  
  # for sample j, calculate the mean running time for each package
  for(t in 1:2)
  {
    A=Time[[j]][,t]
    MTime[[j]][t]=mean(Time[[j]][((A-mean(A))/sd(A))<=1,t])
  }
  
}

# print out the average running time for each package under each sample size (row represents sample size and column represents package)
B=MTime[[1]]
for(i in 2:n_num)
{
  B=rbind(B,MTime[[i]])
}
B=rbind(c("cmprskHD","glmnet"),B)
print(B)
