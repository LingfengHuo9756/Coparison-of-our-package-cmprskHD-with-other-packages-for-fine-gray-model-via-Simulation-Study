# Input: observed data
# Output: weighted pseudo entries with start and stop time
FG.Surv_timedependent=function(surv,cause, Z1, Z2, cov.change, id = 1:nrow(surv))
{
  n=nrow(surv)
  # include left truncation column
  if(ncol(surv)==2)
  {
    surv=Surv(rep(0, n),surv[,1], surv[,2]) 
  }
  
  if(nrow(Z1) != n)
    stop("Z1 and surv should have the same number of rows. ")
  
  if(nrow(Z2) != n)
    stop("Z2 and surv should have the same number of rows. ")
  
  if(length(id) != n)
    stop("Length of id and surv do not match. ")
  
  id.unique = unique(id) # assume it is sorted
  
  # every element of id.list contains the indexes of the adjacent same ids
  id.list = list(1)
  for(i in 2:n)
    if(id[i]==id[i-1])
    {  id.list[[length(id.list)]] = c(id.list[[length(id.list)]],i)
    }else{
      id.list[[length(id.list)+1]] = i
    }
  last.row = sapply(id.list,max)
  
  if(length(cause)==length(id.unique))  ## why do this?
  {
    tmp = rep(0, n)
    tmp[last.row] = cause
    cause = tmp
  }
  
  if(length(cause) != n)
    stop("Length of cause does not match the number of rows in surv or the number of unique ids. ")
  
  if(any(cause * c(surv[,3]) != cause))
    stop("Cause should be coded as zero for censored subjects and non-zero for uncensored subjects. ")
  
  C.Surv=Surv(surv[last.row, 2], 1-surv[last.row, 3])
  C.km <- survfit(C.Surv ~ 1,  type="kaplan-meier")
  C.knot = sort(unique(C.Surv[C.Surv[,2]==1,1])) # find censoring times/ Then find the unique times/ Then sort it in ascending order
  t1.knot = sort(surv[,2][cause==1])
  C.surv = c(summary(C.km, times = C.knot, extend = T)$surv)
  if(min(C.surv) > 0)
  {
    C.knot = c(C.knot,max(C.Surv[,1])+1) # largest observation time +1
    C.surv = c(C.surv, 0) # G(largest censoring time +1)=0
  }
  Ghat=stepfun(C.knot,c(1,C.surv))
  
  FG.start = FG.end = FG.event = FG.id = FG.weight = FG.Z = FG.X = NULL
  
  
  
  for(i in 1:length(id.unique))
  {
    ni = length(id.list[[i]]) # 1
    id.surv = surv[id.list[[i]],] # surv i
    id.Z = id.list[[i]] # id of i
    X = id.surv[ni, 2] # observation time
    
    if( ni > 1)
    {
      if(any(c(id.surv[,3]) != c(rep(0, ni-1),1)))
        stop(paste( "Wrong event indication for id = ",
                    id.unique[i],
                    ". Only the last row for each id may have event. "
                    ,sep=''))
      
      if(any(c(id.surv[,1]) != c(0,id.surv[-ni,2])))
        stop(paste("Wrong coding for start-end time for id = ",
                   id.unique[i],
                   ". The first row of each id must start at 0, ",
                   "and the start of next row must match the end of the previous. ",
                   sep=''))
    }
    
    if(cause[last.row[i]] > 1)
    {
      t1k=t1.knot[t1.knot>X]
      nk=length(t1k)
      omega = Ghat(t1k)/Ghat(X)
      if(cov.change[last.row[i]]<X)
      {
        if(nk==0)
        {
          FG.start = c(FG.start, 0, cov.change[last.row[i]])
          FG.end = c(FG.end, cov.change[last.row[i]],X)
          FG.event = c(FG.event, rep(0, 2))
          FG.weight = c(FG.weight, rep(1,2))
          FG.id = c(FG.id, rep(id.unique[i],2)) # give id to pseudo entries
          Z.before.change=Z1[id.Z,]
          Z.after.change=do.call("rbind", replicate(1, Z2[id.Z,], simplify = FALSE))
          Z.i=rbind(Z.before.change,Z.after.change)
          FG.Z=rbind(FG.Z,Z.i)
          FG.X = c(FG.X, rep(X,2))
        }
        else{
          FG.start = c(FG.start, 0, cov.change[last.row[i]],X,t1k[-nk])
          FG.end = c(FG.end, cov.change[last.row[i]],X,t1k)
          FG.event = c(FG.event, rep(0, nk+2))
          FG.weight = c(FG.weight, rep(1,2), omega)
          FG.id = c(FG.id, rep(id.unique[i],nk+2)) # give id to pseudo entries
          Z.before.change=Z1[id.Z,]
          Z.after.change=do.call("rbind", replicate(nk+1, Z2[id.Z,], simplify = FALSE))
          Z.i=rbind(Z.before.change,Z.after.change)
          FG.Z=rbind(FG.Z,Z.i)
          FG.X = c(FG.X, rep(X,nk+2))
        }
      }
      else{
        if(nk==0)
        {
          FG.start = c(FG.start, 0)
          FG.end = c(FG.end, X)
          FG.event = c(FG.event, rep(0, nk+1))
          FG.weight = c(FG.weight, 1)
          FG.id = c(FG.id, rep(id.unique[i],nk+1)) # give id to pseudo entries
          Z.before.change=Z1[id.Z,]
          Z.after.change=NULL
          Z.i=rbind(Z.before.change,Z.after.change)
          FG.Z=rbind(FG.Z,Z.i)
          FG.X = c(FG.X, rep(X,nk+1))
        }
        else{
          # identify the smallest tk (cause-1 event time) that the covariate changing time is less than
          k=0
          flag=TRUE
          while(flag)
          {
            k=k+1
            if(t1k[k]>cov.change[last.row[i]])
            {
              flag=FALSE
            }
            if(k==nk && flag==TRUE)
            {
              k=k+1
              flag=FALSE
            }
          }
          
          FG.start = c(FG.start, 0, X,t1k[-nk])
          FG.end = c(FG.end, X,t1k)
          FG.event = c(FG.event, rep(0, nk+1))
          FG.weight = c(FG.weight, 1, omega)
          FG.id = c(FG.id, rep(id.unique[i],nk+1)) # give id to pseudo entries
          Z.before.change=do.call("rbind", replicate((k-1)+1, Z1[id.Z,], simplify = FALSE))
          Z.after.change=do.call("rbind", replicate(nk-(k-1), Z2[id.Z,], simplify = FALSE))
          Z.i=rbind(Z.before.change,Z.after.change)
          FG.Z=rbind(FG.Z,Z.i)
          FG.X = c(FG.X, rep(X,nk+1))
        }
      }
      
    }
    else{
      if(cov.change[last.row[i]]<X)
      {
        FG.start = c(FG.start, 0, cov.change[last.row[i]])
        FG.end = c(FG.end, cov.change[last.row[i]], X)
        if(cause[last.row[i]] == 1)
        {
          FG.event = c(FG.event, c(0,1)) 
        }
        else{
          FG.event = c(FG.event, rep(0,2)) 
        }
        FG.weight = c(FG.weight, rep(1,2))
        #   FG.omega = c(FG.omega, rep(1,ni))
        FG.id = c(FG.id, rep(id.unique[i],2))
        Z.before.change=Z1[id.Z,]
        Z.after.change=Z2[id.Z,]
        Zi=rbind(Z.before.change,Z.after.change)
        FG.Z = rbind(FG.Z, Zi)
        FG.X = c(FG.X, rep(X,2))
      }
      else{
        FG.start = c(FG.start, 0)
        FG.end = c(FG.end, X)
        if(cause[last.row[i]] == 1)
        {
          FG.event = c(FG.event, 1)
        }
        else{
          FG.event = c(FG.event, 0)
        }
        FG.weight = c(FG.weight, 1)
        #   FG.omega = c(FG.omega, rep(1,ni))
        FG.id = c(FG.id, id.unique[i])
        Z.before.change=Z1[id.Z,]
        Zi=Z.before.change
        FG.Z = rbind(FG.Z, Zi)
        FG.X = c(FG.X, X)
      }
      
    }
  }
  
  out.surv= Surv(FG.start, FG.end, FG.event)
  out.Z = FG.Z
  
  return(list(surv = Surv(FG.start, FG.end, FG.event), weight = FG.weight, Z = out.Z,
              X=FG.X, id = FG.id, C.surv = C.Surv))
}
