codedir <- '/N/project/Covid/scr/'
source(paste0(codedir,'covidcontrol.R'))
Rcpp::sourceCpp(paste0(codedir,"Block.cpp"))

act <- function(i,dfsir,df,ncl){
  if(df[i,1] == 'R'){
    indx <- 3
  }
  else if(df[i,1] == 'D'){
    indx <- 4}
  else if(df[i,1] == 'I'){
    dfsir$S[(df[i,2]+1):ncl] <- dfsir$S[(df[i,2]+1):ncl]-1
    dfsir$I[(df[i,2]+1):ncl] <- dfsir$I[(df[i,2]+1):ncl]+1
    dfsir$pI[(df[i,2]+1):ncl] <- dfsir$pI[(df[i,2]+1):ncl]+df[i,4]
    dfsir$pD[(df[i,2]+1):ncl] <- dfsir$pD[(df[i,2]+1):ncl]+df[i,5]
    return(dfsir)
  }
  else{
    dfsir$S[(df[i,2]+1):ncl] <- dfsir$S[(df[i,2]+1):ncl]-1
    return(dfsir)
  }
  dfsir$S[(df[i,2]+1):ncl] <- dfsir$S[(df[i,2]+1):ncl]-1
  dfsir$I[(df[i,2]+1):(df[i,3]+1)] <- dfsir$I[(df[i,2]+1):(df[i,3]+1)]+1
  dfsir$pI[(df[i,2]+1):(df[i,3]+1)] <- dfsir$pI[(df[i,2]+1):(df[i,3]+1)]+df[i,4]
  dfsir$pD[(df[i,2]+1):(df[i,3]+1)] <- dfsir$pD[(df[i,2]+1):(df[i,3]+1)]+df[i,5]
  dfsir[(df[i,3]+2):ncl,indx] <- dfsir[(df[i,3]+2):ncl,indx]+1
  return(dfsir)
}

dfSIR <- function(df,N){
  ncl <- max(c(df[,3],df[,2]),na.rm = T)+2
  dfsir <- data.frame(S = rep(N,ncl),I = rep(0,ncl),R = rep(0,ncl), 
                      D = rep(0,ncl),pI = rep(0,ncl),pD = rep(0,ncl))
  
  for(i in 1:nrow(df) ){dfsir <- act(i,dfsir,df,ncl)}
  dfsir$pI[1:(ncl-1)] <- dfsir$pI[1:(ncl-1)]/dfsir$I[1:(ncl-1)]
  dfsir$pD[1:(ncl-1)] <- dfsir$pD[1:(ncl-1)]/dfsir$I[1:(ncl-1)]
  return(dfsir)  
} 

dfpreSIR <- function(df){
  df1 <- df[which(df$typ != 'S'),c(2,3,4,9,10)]
  df1[,2] <- ceiling(df1[,2])
  df1[,3] <- ceiling(df1[,3])
  return(df1)
}


#pDtopIA <- function(ipI,ipD,pI,pD){
#  return(((0.9-pI)/(-pD))*(ipD-pD)+pI)
#}

#pItopDA <- function(ipI,ipD,pI,pD){
#  return(((-pD)/(0.9-pI))*(ipI-pI)+pD)
#}

#pDtopIM <- function(ipI,ipD,pI,pD){
#  return(exp((log(0.9)-log(pI))/(-pD))*(ipD-pD)+log(pI))
#}

#pItopDM <- function(ipI,ipD,pI,pD){
#  return(((-pD)/(log(0.9)-log(pI)))*(log(ipI)-log(pI))+pD)
#}

# muta<- function(ipI,ipD,pI,pD,rR,mfit,mfit2,pIM,pDM,pC,pIMB,pDMB,funmut){
#   if(rbinom(1,1,pC) == 1){
#     mutPI <- rbinom(2,1,c(pIM,pDM))
#     if(pIMB == 'B'){
#       aux <- rbinom(1,1,pDMB)
#       mutS <- c(aux,1-aux)
#     }
#     else{
#       mutS <- rep(rbinom(1,1,pIMB),2) 
#     }
#   }
#   else{
#     mutPI <- rbinom(2,1,c(pIM,pDM))
#     mutS <-  rbinom(2,1,c(pIMB,pDMB))
#   }
#   return(funmut(mutPI,mutS,ipI,ipD,pI,pD,rR,mfit,mfit2))
# }

muta<- function(ipI,ipD,pI,pD,rR,mfit,mfit2,pIM,pDM,pC,pIMB,pDMB,funmut){
  if(rbinom(1,1,pC) == 1){
    mutPI <- rep(rbinom(1,1,pIM),2)
    if(pIMB == 'B'){
      aux <- rbinom(1,1,pDMB)
      mutS <- c(aux,1-aux)
    }
    else{
      mutS <- rbinom(2,1,c(pIMB,pDMB)) 
    }
  }
  else{
    mutPI <- rbinom(2,1,c(pIM,pDM))
    mutS <-  rbinom(2,1,c(pIMB,pDMB))
  }
  return(funmut(mutPI,mutS,ipI,ipD,pI,pD,rR,mfit,mfit2))
}


funmutA<-function(mutPI,mutS,ipI,ipD,pI,pD,rR,mfit,mfit2){
  if(mutPI[1] == 1){
    ipI <- max(min(ipI+mfit*(2*mutS[1]-1),1-mfit),mfit)
  }
  if(mutPI[2] == 1){
    ipD <- max(min(ipD+mfit2*(2*mutS[2]-1),1),0)
  }
  return(c(ipI,ipD))
}


# funmutM<-function(mutPI,mutS,ipI,ipD,pI,pD,rR,mfit){
#   if(mutPI[1] == 1){
#     ipI <- max(min(ipI+mfit*(2*mutS[1]-1),1-mfit),mfit)
#   }
#   if(mutPI[2] == 1){
#     ipD <- max(min(ipD+mfit2*(2*mutS[2]-1),1-mfit2),mfit2)
#   }
#   return(c(ipI,ipD))
# }

covid_mut <- function(M,mc,N,d,lamb,pD,rR,pI,t,q,nIs,pre,ij,
                      mfit,pIM,pDM,pC,pIMB,pDMB,funmut,mut){
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
                     pI = rep(NA,N),
                     pD = rep(NA,N),
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
    df$pI[li] <- pI
    df$pD[li] <- pD
    str <- paste0(' ',li,' ;')
    t <- 0
    nI <-1
    while(T){
      lenI <- length(lrI)
      rateInfect <- sum(lrI*lsl)
      rateReco <- rR*lenI
      ratemin <- rateInfect + rateReco
      if(ratemin == 0){
        break
      }
      else{
        texp <- rexp(1,rate = ratemin)
        if(tStop < t){
          break
        }
        t <- t+texp
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
          pID <- muta(df$pI[df$f[ni]],df$pD[df$f[ni]],pI,pD,rR,mfit,mfit2,
                      pIM,pDM,pC,pIMB,pDMB,funmut)
          df$pI[ni] <- pID[1]
          df$pD[ni] <- pID[2]
          lrI <- c(lrI, adjustrI(df$pI[ni],rR) )
          ls = append(ls,list(nb[[ni]][df$typ[nb[[ni]]] == 'S']))
          lsl = c(lsl,length(nb[[ni]][df$typ[nb[[ni]]] == 'S']))
        }
        else{
          ix <- if(length(lrI)>1){sample(1:lenI,1)}else{1}
          df$typ[li[ix]] <-  if(rbinom(1,1, pD )==1){'D'}else{'R'}
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
    hiaux <- c(hiaux,sum(df$typ != 'S'))
    if(nI  > nIs || t>tStop){
      break
    }
  }
  write.table(str, 
              file = paste0(pre,N,'/M_',mut,'_',pC,'_',pIM,'_',pDM,'_',pIMB,'_',pDMB,'/','NW','/',ij,'.txt'),
              sep = "\t",row.names = FALSE,col.names = FALSE,append=F)
  write.table(data.frame(hi = hiaux,lvl = pre, N= N), 
              file = paste0(pre,N,'/M_',mut,'_',pC,'_',pIM,'_',pDM,'_',pIMB,'_',pDMB,'/','HI','/',ij,'.txt'),
              sep = "\t",row.names = FALSE,col.names = FALSE,append=F)
  write.table(df, 
              file = paste0(pre,N,'/M_',mut,'_',pC,'_',pIM,'_',pDM,'_',pIMB,'_',pDMB,'/','DF','/',ij,'.txt'),
              sep = "\t",row.names = FALSE,col.names = FALSE,append=F)
  write.table(igraph::as_edgelist(g),
              file = paste0(pre,N,'/M_',mut,'_',pC,'_',pIM,'_',pDM,'_',pIMB,'_',pDMB,'/G','/',ij,'.txt'),sep='\t',
              row.names = F,col.names = F,append = T)
  write.table(which(df$typ == 'S'), 
              file = paste0(pre,N,'/M_',mut,'_',pC,'_',pIM,'_',pDM,'_',pIMB,'_',pDMB,'/','G','/S',ij,'.txt'),
              sep = "\t",row.names = FALSE,col.names = FALSE,append=F)
  write.table(data.frame(dg = dg,lvl = pre, N= N ),
              file = paste0(pre,N,'/M_',mut,'_',pC,'_',pIM,'_',pDM,'_',pIMB,'_',pDMB,'/','DG','/',ij,'.txt'),
              sep ='\t', row.names = F,col.names = F,append = F)
  write.table(dfSIR(dfpreSIR(df),N),
               file = paste0(pre,N,'/M_',mut,'_',pC,'_',pIM,'_',pDM,'_',pIMB,'_',pDMB,'/','SIR','/',ij,'.txt'),
               sep = "\t",row.names = FALSE,col.names = FALSE,append=F)
}

COVID_muta <- function(M,mc,N,d,lamb,pD,rR,pI,t,q,nIs,pre,
                       mfit,pIM,pDM,pC,pIMB,pDMB,funmut,mut){
  dir.create(paste0(pre,N),showWarnings = F)
  dir.create(paste0(pre,N,'/M_',mut,'_',pC,'_',pIM,'_',pDM,'_',pIMB,'_',pDMB))
  dir.create(paste0(pre,N,'/M_',mut,'_',pC,'_',pIM,'_',pDM,'_',pIMB,'_',pDMB,'/','HI'))
  dir.create(paste0(pre,N,'/M_',mut,'_',pC,'_',pIM,'_',pDM,'_',pIMB,'_',pDMB,'/','DG'))
  dir.create(paste0(pre,N,'/M_',mut,'_',pC,'_',pIM,'_',pDM,'_',pIMB,'_',pDMB,'/','G'))
  dir.create(paste0(pre,N,'/M_',mut,'_',pC,'_',pIM,'_',pDM,'_',pIMB,'_',pDMB,'/','NW'))
  dir.create(paste0(pre,N,'/M_',mut,'_',pC,'_',pIM,'_',pDM,'_',pIMB,'_',pDMB,'/','DF'))
  dir.create(paste0(pre,N,'/M_',mut,'_',pC,'_',pIM,'_',pDM,'_',pIMB,'_',pDMB,'/','SIR'))
  if(mc == 0){
    invisible(lapply(1:M, function(i) covid_mut(M,mc,N,d,lamb,pD,rR,pI,t,q,nIs,pre,i,
                                                            mfit,pIM,pDM,pC,pIMB,pDMB,funmut,mut) ))
  }
  else(
    invisible(parallel::mclapply(1:M,function(i) covid_mut(M,mc,N,d,lamb,pD,rR,pI,t,q,nIs,pre,i,
                                                                       mfit,pIM,pDM,pC,pIMB,pDMB,funmut,mut),mc.cores = mc))
  )
}
