### Vaccination
# Parameters: 
# tV: Number of individuals needed to start vaccination
# D:  Number of doses needed  
# oV: Vaccination type 
#     1 - Uniform
#     2 - Neighbours
#     3 - Among most connected
#     4 - Most connected
#     5 - Among least connected
#     6 - Least connected 
# who: Population groups to vaccinate
#     1 - c('S') Susceptible 
#     2 - c('S','I','R') Susceptible, Recovered or Infected

Sample <- function(vec){
  if(length(vec) == 1){
    return(vec[1])}
  else{
    return(sample(vec,1))
  }
}

# Check if it is time to start the vaccination process

check_begin <- function(t,nI,dnI){
  return(nI > dnI)
}

# Set type label S,V,I,R,D

setV <- function(ix,typ,qM,dgi){
  return(sapply(1:length(ix),
                function(i) if(typ[ix[i]] == 'S' & rbinom(1,1,pE) == 1 ){'V'}else{typ[ix[i]]}))
}

# Vaccination strategy 1,3,4,5,6
#     1 - Uniform
#     3 - Among most connected
#     4 - Most connected
#     5 - Among least connected
#     6 - Least connected 
vacc1 <- function(df,nb,D,m,who,N,qM){
  vindx <- which(df$typ %in% who)
  if(length(vindx) < D){
    D <- length(vindx)}
  if(m == 1){
    m <- 1:length(vindx) 
  }
  else if(m > 1){
    m <- 1:ceiling(N/m)
  }
  else if(m < 0){
    m <- tail(1:length(vindx),ceiling(-N/m))
  }
  lnb <- vindx[sort(df$dg[vindx],decreasing = T,index.return = T)$ix][m]
  ix<-if(D==1){lnb[1]}else{sample(lnb,D)}
  df$v[ix] <- T 
  df$typ[ix] <- setV(ix,df$typ,qM,df$dg[ix])
  return(df)
}
                
# Vaccination strategy 2
#     2 - Neighbours
                
vacc2 <- function(df,nb,D,who,qM){
  sindx <- which(df$typ %in% who)
  if(length(sindx) < D){
    D <- length(sindx)}
  while(D>0){
    if(length(sindx)>0){
      aux <- sapply(sindx,function(i) sum(df$typ[nb[[i]]] %in% who  & df$v[nb[[i]]] == F))
      if(sum(aux) == 0){
        ix <- if(D==1){sindx[1]}else{sample(sindx,D)}
        df$v[ix]<- T
        df$typ[ix] <- setV(ix,df$typ,qM,df$dg[ix])
        break
      }
      ix <- Sample(sindx[which(aux>0)])
      ni <- nb[[ix]][which(df$typ[nb[[ix]]] %in% who & df$v[nb[[ix]]] == F)]
      if(length(ni)>0){
        ix <- Sample(ni)
        df$v[ix]<- T
        df$typ[ix] <- if(df$typ[ix] == 'S' & rbinom(1,1,pE) == 1 ){'V'}else{df$typ[ix]}
        D <- D-1
        sindx <- sindx[!sindx %in% ix]
      }
    }
  }
  return(df)
}

#N/m >= D

vaccinate <- function(df,nb,D,oV,who,N,qM){
  if(oV == 1){
    df <- vacc1(df,nb,D,1,who,N,qM)
  }
  else if(oV == 2){
    df <- vacc2(df,nb,D,who,qM)
  }
  else if(oV == 3){
    df <- vacc1(df,nb,D,2,who,N,qM)
  }
  else if(oV == 4){
    df <- vacc1(df,nb,D,N/D,who,N,qM)
  }
  else if(oV == 5){
    df <- vacc1(df,nb,D,-2,who,N,qM)
  }
  else if(oV == 6){
    df <- vacc1(df,nb,D,-N/D,who,N,qM)
  }
  return(df)
}

covid_vacc <- function(M,mc,N,d,lamb,pD,rR,pI,t,q,nIs,pre,tV,D,oV,who,ij){
  hiaux <- c()
  while (T) {
    if(lamb > 0){
      typ <- rep(0,N)
      g <- block_power_matrix(matrix(d/N),lamb,typ)
    }
    else{g <- igraph::erdos.renyi.game(N,d/N)}
    ns <- setdiff(1:N,igraph::V(g))
    if(length(ns)>0){
      g <- igraph::add_vertices(g,length(ns))
    }
    nb <- igraph::adjacent_vertices(g,1:N)
    dg <- igraph::degree(g)
    qM <- quantile(dg,c(q))
    df <- data.frame(ID = 1:N,
                     typ = rep('S',N),
                     tI = rep(NA,N),
                     tR = rep(NA,N),
                     f = rep(NA,N),
                     ch = rep(0,N),
                     lt = rep(NA,N),
                     dg = dg,
                     v = rep(F,N),
                     stringsAsFactors = F)
    k <- 1
    rI <- adjustrI(pI,rR)
    str <- ' '
    li = c(sample(which(dg >0),1))
    lrI = c(rI)
    ls = list(nb[[li]])
    lsl = c(length(nb[[li]]))
    df$f[li] <- 0
    df$typ[li] <- 'I'
    df$lt[li] <- 0
    df$tI[li] <- 0
    str <- paste0(' ',li,' ;')
    t <- 0
    nI <- 1
    vacc <- F
    while(T){
      if(!vacc && check_begin(t,nI,tV)){
        df <- vaccinate(df,nb,D,oV,who,N,qM)
        df$tI[df$typ == 'V'] <- t
        ls <- lapply(li, function(i) nb[[i]][df$typ[nb[[i]]] == 'S'])
        lsl <- if(length(ls)==0){0}else{sapply(ls, function(ni) length(ni))}
        vacc <- T
      }
      lenI <- length(lrI)
      rateInfect <- sum(lrI*lsl)
      rateReco <- rR*lenI
      ratemin <- rateInfect + rateReco
      if(ratemin == 0){
        break
      }
      else{
        texp <- rexp(1,rate = ratemin)
        t <- t+texp
        if(tStop < t){
          break
        }
        if(runif(1) < rateInfect/ratemin){
          nI <- nI+1
          ix <- if(length(lrI)>1){sample(1:lenI,1,prob = rI*lsl/rateInfect)}else{1}
          ni <- if(length(ls[[ix]]) == 1){ls[[ix]][1]}else{sample(ls[[ix]],1)}
          lni <- unlist(lapply(li,function(i){if(ni %in% nb[[i]]){i}}))
          df$f[ni] <- li[ix]
          df$typ[ni] <- 'I'
          df$tI[ni] <- t
          df$ch[df$f[ni]]<-df$ch[df$f[ni]]+1 
          str <- newick(df$f[ni],ni,t-df$lt[df$f[ni]],str)
          df$lt[ni] <- t
          df$lt[df$f[ni]] <-t 
          inix <-which(li %in% lni) 
          ls[inix] = lapply(inix, function(i){ls[[i]][ls[[i]] != ni]})
          lsl[inix] <- lsl[inix]-1
          li <- c(li,ni)
          lrI <- c(lrI,rI)
          ls = append(ls,list(nb[[ni]][df$typ[nb[[ni]]] == 'S']))
          lsl = c(lsl,length(nb[[ni]][df$typ[nb[[ni]]] == 'S']))
        }
        else{
          ix <- if(length(lrI)>1){sample(1:lenI,1)}else{1}
          df$typ[li[ix]] <-  if(rbinom(1,1, probOfDeath(pD,df$dg[li[ix]],qM) ) == 1){'D'}else{'R'}
          df$tR[li[ix]] <- t
          str <- newick(li[ix],F,t-df$lt[li[ix]],str)
          li <- li[-ix]
          ls <- ls[-ix]
          lsl <- lsl[-ix]
          lrI <- lrI[-ix]
        }
      }
      k <- k+1
    }
    hiaux <- c(hiaux,sum(df$typ != 'S' & df$typ != 'V'))
    if(nI > nIs | t >tStop){
      break
    }
  }
  write.table(str, 
              file = paste0(pre,N,'/VAC','_',paste0(who,collapse = ''),'_',oV,'_',tV,'_',D,'/','NW','/',ij,'.txt'),
              sep = "\t",row.names = FALSE,col.names = FALSE,append=F)
  write.table(data.frame(hi = hiaux,lvl = pre, N= N), 
              file = paste0(pre,N,'/VAC','_',paste0(who,collapse = '')
                            ,'_',oV,'_',tV,'_',D,'/','HI','/',ij,'.txt'),
              sep = "\t",row.names = FALSE,col.names = FALSE,append=F)
  write.table(df, 
              paste0(pre,N,'/VAC','_',paste0(who,collapse = '')
                     ,'_',oV,'_',tV,'_',D,'/','DF','/',ij,'.txt'),
              sep = "\t",row.names = FALSE,col.names = FALSE,append=F)
  write.table(which(df$typ == 'S'), 
              file = paste0(pre,N,'/VAC','_',paste0(who,collapse = ''),'_',oV,'_',tV,'_',D,'/','G','/S',ij,'.txt'),
              sep = "\t",row.names = FALSE,col.names = FALSE,append=F)
  write.table(data.frame(dg = dg,lvl = pre, N= N ),
              file = paste0(pre,N,'/VAC','_',paste0(who,collapse = ''),'_',oV,'_',tV,'_',D,'/','DG','/',ij,'.txt'),
              sep ='\t', row.names = F,col.names = F,append = F)
  write.table(igraph::as_edgelist(g),
              file = paste0(pre,N,'/VAC','_',paste0(who,collapse = ''),'_',oV,'_',tV,'_',D,'/','G','/',ij,'.txt'),sep='\t',
              row.names = F,col.names = F,append = F)
  write.table(dfSIR(dfpreSIR(df),N),
              file = paste0(pre,N,'/VAC','_',paste0(who,collapse = ''),'_',oV,'_',tV,'_',D,'/','SIR','/',ij,'.txt'),sep='\t',
              row.names = FALSE,col.names = FALSE,append=F)
}

COVID_vacc <- function(M,mc,N,d,lamb,pD,rR,pI,t,q,nIs,pre,tV,D,oV,who){
  dir.create(paste0(pre,N),showWarnings = F)
  dir.create(paste0(pre,N,'/VAC','_',paste0(who,collapse = ''),'_',oV,'_',tV,'_',D))
  dir.create(paste0(pre,N,'/VAC','_',paste0(who,collapse = ''),'_',oV,'_',tV,'_',D,'/','HI'))
  dir.create(paste0(pre,N,'/VAC','_',paste0(who,collapse = ''),'_',oV,'_',tV,'_',D,'/','DG'))
  dir.create(paste0(pre,N,'/VAC','_',paste0(who,collapse = ''),'_',oV,'_',tV,'_',D,'/','G'))
  dir.create(paste0(pre,N,'/VAC','_',paste0(who,collapse = ''),'_',oV,'_',tV,'_',D,'/','NW'))
  dir.create(paste0(pre,N,'/VAC','_',paste0(who,collapse = ''),'_',oV,'_',tV,'_',D,'/','DF'))
  dir.create(paste0(pre,N,'/VAC','_',paste0(who,collapse = ''),'_',oV,'_',tV,'_',D,'/','SIR'))
  if(mc == 0){
    invisible(lapply(1:M, function(i) 
      covid_vacc(M,mc,N,d,lamb,pD,rR,pI,t,q,nIs,pre,tV,D,oV,who,i)))
  }
  else(
    invisible(parallel::mclapply(1:M,function(i) 
      covid_vacc(M,mc,N,d,lamb,pD,rR,pI,t,q,nIs,pre,tV,D,oV,who,i),mc.cores = mc))
  )
}
