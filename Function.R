# A function to do wilcox.test
tFun <- function(prf,phe,id,p,ex.0=T){
  wi <- which(colnames(phe)==id)
  wp <- which(colnames(phe)==p)
  phe.s <- data.frame(id=phe[,wi],group=phe[,wp])
  phe.s <- phe.s[which(!is.na(phe.s$group)),]
  phe.s$group <- as.factor(phe.s$group)
  pName <- colnames(prf)
  tsum <- t(apply(prf,1,FUN=function(x,pName,phe){
    #pName <- rownames(x)
    pi <- data.frame(id=pName,val=as.numeric(x))
    dat <- merge(pi,phe,by="id")
    dat <- dat[which(!is.na(dat$val)),]
    if(ex.0==T){
      dat <- dat[dat$val>0,]
    }
    
    dat$group <- droplevels(dat$group) 
    NZ <- table(dat$group[which(dat$val>0)])
    md <- ddply(dat,"group",summarize,median=median(val))
    #length(which(dat$va>0))
    if(length(levels(dat$group))==2){
      tt <- wilcox.test(val~group,dat)
      return(c(nrow(dat),NZ[1],NZ[2],md$median[1],md$median[2],tt$p.value))
    }else {
      return(c(0,NZ[1],NZ[2],0,0,NA))
    }
  },pName,phe.s))
  wpn <- unlist(levels(phe[,wp]))
  res <- data.frame(phe=p,mgs.V2=rownames(tsum),tsum)
  colnames(res)[3:8] <- c("len",paste0("Num.",wpn),paste0("Median.",wpn),"p.value")
  res$prevalent <- (res$len)/ncol(prf)
  res$FDR <- NA
  res$FDR[which(!is.na(res$p.value))] <- p.adjust(res$p.value[which(!is.na(res$p.value))],method = "BH")
  if(grepl("MGS",rownames(prf)[1])){
    res <- merge(mgs.V2.anno2,res,by="mgs.V2",all.y=T)
  }
  return(res)
}

# A function to do kruskal.test
kFun <- function(prf,phe,id,p){
  wi <- which(colnames(phe)==id)
  wp <- which(colnames(phe)==p)
  phe.s <- data.frame(id=phe[,wi],group=phe[,wp])
  phe.s <- phe.s[which(!is.na(phe.s$group)),]
  phe.s$group <- as.factor(phe.s$group)
  pName <- colnames(prf)
  tsum <- t(apply(prf,1,FUN=function(x,pName,phe){
    pi <- data.frame(id=pName,val=as.numeric(x))
    dat <- merge(pi,phe,by="id")
    dat <- dat[which(!is.na(dat$val)),]
    #dat <- dat[dat$val>0,]
    NZ <- table(dat$group[which(dat$val>0)])
    md <- ddply(dat,"group",summarize,median=median(val))
    tt <- kruskal.test(val~group,dat)
    return(c(nrow(dat),NZ[1],NZ[2],NZ[3],md$median[1],md$median[2],md$median[3],tt$p.value))
  },pName,phe.s))
  wpn <- unlist(levels(phe.s$group))
  res <- data.frame(phe=p,mgs.V2=rownames(tsum),tsum)
  colnames(res)[3:10] <- c("len",paste0("Num.",wpn),paste0("Median.",wpn),"p.value")
  res$prevalent <- (res$len)/ncol(prf)
  res$FDR <- NA
  res$FDR[which(!is.na(res$p.value))] <- p.adjust(res$p.value[which(!is.na(res$p.value))],method = "BH")
  if(grepl("MGS",rownames(prf)[1])){
    res <- merge(mgs.V2.anno2,res,by="mgs.V2",all.y=T)
  }
  return(res)
}

#update the wilcox.test function to do test within divided sub-group

ttFun <- function(prf,phe,id,p,block,ex.0=T){
  wi <- which(colnames(phe)==id)
  wp <- which(colnames(phe)==p)
  wb <- which(colnames(phe)==block)
  bGrp <- levels(phe[,wb])
  phe.s <- data.frame(id=phe[,wi],group=phe[,wp],block=phe[,wb])
  pName <- colnames(prf)
  wilcox.fun <- function(x,pName,phe){
    pi <- data.frame(id=pName,val=as.numeric(x))
    dat <- merge(pi,phe,by="id")
    dat <- dat[which(!is.na(dat$group)),]
    Num <- nrow(dat)
    dat <- dat[which(!is.na(dat$val)),]
    if(ex.0){
      dat <- dat[dat$val>0,]
    }
    dat$group <- droplevels(dat$group) 
    NZ <- table(dat$group[which(dat$val>0)])
    md <- ddply(dat,"group",summarize,median=median(val))
    #length(which(dat$va>0))
    if(length(levels(dat$group))==2){
      tt <- wilcox.test(val~group,dat)
      return(c(Num,NZ[1],NZ[2],md$median[1],md$median[2],tt$p.value))
    }else{
      return(c(0,NZ[1],NZ[2],0,0,NA))
    }
  }
  
  curateFun <- function(prf,phe,block){
    gN <- unlist(levels(phe$group))
    tsum <- t(apply(prf,1,FUN=wilcox.fun,pName,phe))
    colnames(tsum) <- c("Num",gN,paste0("median@",gN),"p.value")
    
    res <- data.frame(mgs.V2=dimnames(tsum)[[1]],block=block,tsum)
    res$Occ. <- (res[,4]+res[,5])/res$Num
    res$Enrich <- ifelse(res[,7]-res[,6]>0,gN[2],gN[1])
    res$FDR <- NA
    res$FDR[which(!is.na(res$p.value))] <- p.adjust(res$p.value[which(!is.na(res$p.value))],method = "BH")
    if(grepl("MGS",rownames(prf)[1])){
      res <- merge(mgs.V2.anno2,res[,c(1:5,9,6,7,10,8,11)],by="mgs.V2",all.y=T)
    }
    return(res)
  }
  tsum0 <- curateFun(prf,phe.s,"NULL")
  tsum1 <- curateFun(prf,phe.s[which(phe.s$block==bGrp[1]),],bGrp[1])
  tsum2 <- curateFun(prf,phe.s[which(phe.s$block==bGrp[2]),],bGrp[2])
  
  res <- cbind(tsum0,tsum1[,-c(1:3)],tsum2[,-c(1:3)])
  return(res)
  
}


# update the kruskal.test function to do test within divided sub-group

tkFun <- function(prf,phe,id,p,block,ex.0=T){
  wi <- which(colnames(phe)==id)
  wp <- which(colnames(phe)==p)
  wb <- which(colnames(phe)==block)
  bGrp <- levels(phe[,wb])
  phe.s <- data.frame(id=phe[,wi],group=phe[,wp],block=phe[,wb])
  phe.s$group <- as.factor(phe.s$group)
  pName <- colnames(prf)
  kw.fun <- function(x,pName,phe){
    pi <- data.frame(id=pName,val=as.numeric(x))
    dat <- merge(pi,phe,by="id")
    dat <- dat[which(!is.na(dat$group)),]
    Num <- nrow(dat)
    dat <- dat[which(!is.na(dat$val)),]
    dat$group <- as.factor(dat$group)
    if(ex.0){
      dat <- dat[dat$val>0,]
    }
    #dat$group <- droplevels(dat$group) 
    NZ <- table(dat$group[which(dat$val>0)])
    md <- ddply(dat,"group",summarize,median=median(val))
    #length(which(dat$va>0))
    if(length(levels(dat$group))==3){
      tk <- kruskal.test(val~group,dat)
      return(c(Num,NZ[1],NZ[2],NZ[3],md$median[1],md$median[2],md$median[3],tk$p.value))
    }else{
      return(c(0,NZ[1],NZ[2],NZ[3],0,0,0,NA))
    }
  }
  
  curateFun <- function(prf,phe,block){
    gN <- unlist(levels(phe$group))
    tsum <- t(apply(prf,1,FUN=kw.fun,pName,phe))
    colnames(tsum) <- c("Num",gN,paste0("median@",gN),"p.value")
    
    res <- data.frame(mgs.V2=dimnames(tsum)[[1]],block=block,tsum)
    res$Occ. <- (res[,4]+res[,5]+res[,6])/res$Num
    
    res$Enrich <- gN[as.numeric(apply(res[,7:9],1, function(x) {
      b <- which(x==max(x[which(!is.na(x))]))
      if(length(b)>1){return(NA)}else{return(b)}
    }))]
    res$FDR <- NA
    res$FDR[which(!is.na(res$p.value))] <- p.adjust(res$p.value[which(!is.na(res$p.value))],method = "BH")
    if(grepl("MGS",rownames(prf)[1])){
      res <- merge(mgs.V2.anno2,res[,c(1:6,11,7:9,12,10,13)],by="mgs.V2",all.y=T)
    }
    return(res)
  }
  tsum0 <- curateFun(prf,phe.s,"NULL")
  tsum1 <- curateFun(prf,phe.s[which(phe.s$block==bGrp[1]),],bGrp[1])
  tsum2 <- curateFun(prf,phe.s[which(phe.s$block==bGrp[2]),],bGrp[2])
  
  res <- cbind(tsum0,tsum1[,-c(1:3)],tsum2[,-c(1:3)])
  return(res)
  
}


# Modify the wilcox.test function to do pair test

ptFun <- function(prf1,phe1,tag1,prf2,phe2,tag2,id,ID,blockA,blockB,ex0=T,ext=F){
  # find index in phe
  wi <- c(which(colnames(phe1)==id),which(colnames(phe2)==id))
  wI <- c(which(colnames(phe1)==ID),which(colnames(phe2)==ID))
  wbA <- c(which(colnames(phe1)==blockA),which(colnames(phe2)==blockA))
  wbB <- c(which(colnames(phe1)==blockB),which(colnames(phe2)==blockB))
  
  phe <- rbind(data.frame(id=phe1[,wi[1]],ID=phe1[,wI[1]],
                          blockA=phe1[,wbA[1]],blockB=phe1[,wbB[1]],group=tag1,
                          id.opt=paste0(phe1[,wi[1]],"_",tag1)),
               data.frame(id=phe2[,wi[2]],ID=phe2[,wI[2]],
                          blockA=phe2[,wbA[2]],blockB=phe2[,wbB[2]],group=tag2,
                          id.opt=paste0(phe2[,wi[2]],"_",tag2)))
  #prf <- merge(data.frame(id=rownames(prf1),prf1),
  #               data.frame(id=rownames(prf2),prf2),by="id")
  if(length(which(table(phe$id)>1))>0){
    phe$id<-phe$id.opt
    colnames(prf1) <- paste0(colnames(prf1),"_",tag1)
    colnames(prf2) <- paste0(colnames(prf2),"_",tag2)
  }
  prf <- cbind(prf1,prf2)
  bAgrp <- levels(phe$blockA)
  bBgrp <- levels(phe$blockB)
  pName <- colnames(prf)
  phe.n <- phe[phe$id%in%pName,]
  ID.s <- names(which(table(phe.n$ID)==2))
  phe.s <- phe[phe$ID%in%ID.s,]
  wilcox.fun <- function(x,pName,phe,b){
    pi <- data.frame(id=pName,val=as.numeric(x))
    dat <- merge(pi,phe,by="id")
    mdat <- dcast(dat,ID~group,value.var = "val")
    mdats <- mdat[which(!is.na(mdat[,2]+mdat[,3])),]
    if(ex0){
      mdats <- mdat[which(mdat[,2]+mdat[,3]>0),]
    }
    
    dat$group <- droplevels(dat$group) 
    NZ <- nrow(mdats)
    md <- c(median(mdats[,2]),median(mdats[,3]))
    #length(which(dat$va>0))
    if(NZ>0){
      tt <- wilcox.test(mdats[,2],mdats[,3],paired=T)
      return(c(nrow(mdat),NZ,md[1],md[2],tt$p.value))
    }else{
      return(c(nrow(mdat),NZ,md[1],md[2],NA))
    }
  }
  curateFun <- function(prf,phe,block){
    gN <- unlist(levels(phe$group))
    tsum <- t(apply(prf,1,FUN=wilcox.fun,pName,phe))
    colnames(tsum) <- c("pair","NZ",paste0("median@",gN),"p.value")
    
    res <- data.frame(mgs.V2=dimnames(tsum)[[1]],block=block,tsum)
    res$Occ. <- res$NZ/res$pair
    res$Enrich <- ifelse(res[,6]-res[,5]>0,gN[2],gN[1])
    res$FDR <- NA
    res$FDR[!is.na(res$p.value)] <- p.adjust(res$p.value[!is.na(res$p.value)],method = "BH")
    return(res[,c(1:4,8,5,6,9,7,10)])
  }
  tsum0 <- curateFun(prf,phe.s,"NULL")
  tsumA1 <- curateFun(prf,phe.s[which(phe.s$blockA==bAgrp[1]),],bAgrp[1])
  tsumA2 <- curateFun(prf,phe.s[which(phe.s$blockA==bAgrp[2]),],bAgrp[2])
  tsumB1 <- curateFun(prf,phe.s[which(phe.s$blockB==bBgrp[1]),],bBgrp[1])
  tsumB2 <- curateFun(prf,phe.s[which(phe.s$blockB==bBgrp[2]),],bBgrp[2])
  
  res <- cbind(tsum0,tsumA1[,-c(1)],tsumA2[,-c(1)],tsumB1[,-c(1)],tsumB2[,-c(1)])
  
  #IF Extend
  if(ext){
    tsumA1B1 <- curateFun(prf,phe.s[which(phe.s$blockA==bAgrp[1]&
                                            phe.s$blockB==bBgrp[1]),],
                          paste0(bAgrp[1],"&",bBgrp[1]))
    tsumA1B2 <- curateFun(prf,phe.s[which(phe.s$blockA==bAgrp[1]&
                                            phe.s$blockB==bBgrp[2]),],
                          paste0(bAgrp[1],"&",bBgrp[2]))
    tsumA2B1 <- curateFun(prf,phe.s[which(phe.s$blockA==bAgrp[2]&
                                            phe.s$blockB==bBgrp[1]),],
                          paste0(bAgrp[2],"&",bBgrp[1]))
    tsumA2B2 <- curateFun(prf,phe.s[which(phe.s$blockA==bAgrp[2]&
                                            phe.s$blockB==bBgrp[2]),],
                          paste0(bAgrp[2],"&",bBgrp[2]))
    res <- cbind(res,tsumA1B1[,-1],tsumA1B2[,-1],tsumA2B1[,-1],tsumA2B2[,-1])
  }
  if(grepl("MGS",rownames(res)[1])){
    cname <- colnames(res)[-1]
    res <- merge(mgs.V2.anno2,res,by="mgs.V2",all.y=T)
    colnames(res) <- c(colnames(mgs.V2.anno2),cname)
  }
  return(res)
}
