## --------------------------------randomForest------------------------##
train.module <- function(prf,phe,ID,grp,size=NULL,max.marker=30,max.cv=0.4,
                         thread=5,rept=5,steps=2,cv.fold=5,rrf=10,block=NULL,
                         candy=NULL,pick=NULL,dw=NULL,ex.PID=NULL,mtry=3,
                         strati=NULL,pct=NULL){
  #step1
  min.cv <- 1
  marker.num=0
  # parallel method to speed up the process
  cl <- makeCluster(thread)
  max.try <- mtry
  tried = 0
  tmp.train <- NULL
  tmp.cv <- 1
  while((marker.num>max.marker|marker.num==0|min.cv > max.cv) & (tried <= max.try)){
    if(tried == max.try){
      tried <- "done"
      train.cv <- tmp.train
      
      error.cv <- sapply(train.cv, "[[", "error.cv")
      error.cv.rm <- rowMeans(error.cv)
      error.cv.sd <- apply(error.cv,1,sd)
      id <- max(which(error.cv.rm < min(error.cv.sd + error.cv.rm)))
      min.cv <- error.cv.rm[id]
      
      marker.num <- min(as.numeric(names(error.cv.rm)[id]))
      
      mode <- data.frame(num=train.cv[[1]]$n.var,error.cv)
      mmode <- melt(mode,id.vars = "num",variable.name = "times",value.name = "Err.cv")
      mmode$num  <- as.factor(mmode$num)
      
      smode <- data.frame(num=train.cv[[1]]$n.var,times="rm",Err.cv=error.cv.rm)
      smode$num <- as.factor(smode$num)
    }else{
      tried <- tried + 1
      
      if(is.null(strati)){
        tt <- aveSam(prf,phe,ID,grp,NULL,candy,
                     size,pick,ex.PID=ex.PID)
      }else{
        tt <- stratiSam(prf,phe,ID,grp,NULL,candy,
                        pct=pct,strati=strati,pick,ex.PID=ex.PID)
      }
      
      if(class(tt$train.y)=="factor"){ tt$train.y <- droplevels(tt$train.y) }
      
      #clusterEvalQ(cl,source("rfcv1.R"))
      clusterEvalQ(cl,library(randomForest))
      clusterExport(cl,c("tt","rfcv1"))
      train.cv <- parSapply(cl, 1:rept, function(i,...) {
        set.seed(i)
        rfcv1(tt$train.x, tt$train.y, cv.fold = cv.fold, step = steps, ...)
      },simplify = F)
      #
      error.cv <- sapply(train.cv, "[[", "error.cv")
      error.cv.rm <- rowMeans(error.cv)
      error.cv.sd <- apply(error.cv,1,sd)
      id <- max(which(error.cv.rm < min(error.cv.sd + error.cv.rm)))
      min.cv <- error.cv.rm[id]
      if(min.cv < tmp.cv){ 
        tmp.train <- train.cv
        tmp.cv <- min.cv
      }
      
      marker.num <- min(as.numeric(names(error.cv.rm)[id]))
      
      mode <- data.frame(num=train.cv[[1]]$n.var,error.cv)
      mmode <- melt(mode,id.vars = "num",variable.name = "times",value.name = "Err.cv")
      mmode$num  <- as.factor(mmode$num)
      
      smode <- data.frame(num=train.cv[[1]]$n.var,times="rm",Err.cv=error.cv.rm)
      smode$num <- as.factor(smode$num)
    }
    
    print(paste0("marker.num = ",marker.num," | min.cv = ",round(min.cv,2)," | tried: ",tried))
  }
  stopCluster(cl)

  #step 2
  marker.t <- table(unlist(lapply(train.cv, function(x) {
    lapply(x$res, "[", 1:marker.num)
  })))
  marker.t <- sort(marker.t, d = T)

  names(marker.t) <- colnames(tt$train.x)[as.numeric(names(marker.t))]
  marker.p <- names(marker.t)[1:marker.num]
  if(marker.num > 1){
    marker.anno <- as.data.frame(marker.t[1:marker.num])
  }else{
    marker.anno <- data.frame(id=names(marker.t[1:marker.num]),
                              marker.t[1:marker.num])
  }
  if(exists('mgs.V2.anno2')){
    colnames(marker.anno) <- c("mgs.V2","Freq")
    marker.anno <- merge(mgs.V2.anno2,marker.anno,by="mgs.V2",all.y=T)
  }else{
    colnames(marker.anno) <- c("ID","Freq")
  }


  rf.dat <- cbind(tt$train.x[, marker.p],grp=tt$train.y)
  train.rf <- randomForest(grp~.,data=rf.dat, importance = T,proximity=T)
  train.p <- NULL
  if(class(tt$train.y)=="factor"){
    for(i in 2:rrf){
      tmp.rf <- randomForest(grp~.,data=rf.dat, importance = T,proximity=T)
      if( mean(tmp.rf$err.rate[,1]) > mean(train.rf$err.rate[,1]) ){
        train.rf <- tmp.rf
      }
    }
    train.p <- predict(train.rf, type = "prob")
    train.dat <- data.frame(grp=tt$train.y,prob=train.p[,2])
  }else{
    train.p <- predict(train.rf)
    train.dat <- data.frame(grp=tt$train.y,prob=train.p)
  }
  
  
  #

  return(list(set=tt,cv=train.cv,mf=mmode,sf=smode,marker.num=marker.num,
              marker=marker.anno,rf=train.rf,train.pred=train.dat))
}

train.module2 <- function(prf,phe,ID,grp,size,
                         max.marker=30,max.cv=0.7,thread=5,rept=5,steps=2,cv.fold=5,
                         sam=F,rrf=10,smote=F,maxtry=10,impGrp=NULL,trees=500,
                         block=NULL,candy=NULL,pick=NULL,dw=NULL,ex.PID=NULL){
  #step1
  PD.min.cv <- 1
  min.cv <- 1
  marker.num=0
  tmp.min.cv <- 1
  ntree = trees
  # parallel method to speed up the process
  cl <- makeCluster(thread)
  clusterEvalQ(cl,library(randomForest))
  clusterExport(cl,c("tt","rfcv1"))
  tried <- 0
  while((marker.num>max.marker|marker.num==0|PD.min.cv > max.cv) && tried <= maxtry){
    #aveSam(prf,phe,ID,grp,block=NULL,candy=NULL,size,pick=NULL,dw=NULL,ex.ATB=F)
    tried = tried + 1
    tt <- aveSam(prf,phe,ID,grp,NULL,candy,
                 size,pick,ex.PID=ex.PID,sam=sam,smote=smote)
    
    tt$train.y <-droplevels(tt$train.y)
    
    if(tried > maxtry){
      train.cv <- min.train.cv
      tried <- "end"
    }else{
      train.cv <- parSapply(cl, 0:rept, function(i,...) {
        set.seed(i)
        rfcv1(tt$train.x, tt$train.y, cv.fold = cv.fold, step = steps, impGrp = impGrp, ntree = trees)
      },simplify = F)
    }
    #
    if(is.null(impGrp)){
      error.cv <- sapply(train.cv, "[[", "error.cv")
    }else{
      imp.err.cv <- t(matrix(unlist(lapply(train.cv,function(x){
        sapply(x$predicted,function(y) {
          mean((y!=tt$train.y)[which(tt$train.y==impGrp)])
        })
      })),ncol=length(train.cv[[1]]$n.var),byrow=T))
      rownames(imp.err.cv) <- names(train.cv[[1]]$predicted)
      error.cv <- imp.err.cv
    }
    
    error.cv.rm <- rowMeans(error.cv)
    min.cv <- min(error.cv.rm)
    id <- which(error.cv.rm == min(error.cv.rm))
    
    marker.num <- min(as.numeric(names(error.cv.rm)[id]))
    
    mode <- data.frame(num=train.cv[[1]]$n.var,error.cv)
    mmode <- melt(mode,id.vars = "num",variable.name = "times",value.name = "Err.cv")
    mmode$num  <- as.factor(mmode$num)
    
    smode <- data.frame(num=train.cv[[1]]$n.var,times="rm",Err.cv=error.cv.rm)
    smode$num <- as.factor(smode$num)
    if(length(which(grepl("PD",levels(tt$train.y))))>0){
      PD.cv <- reCalCv(train=NULL,train.cv=train.cv,grp="PD",train.mf=mmode,train.sf=smode,train.y=tt$train.y,rep=NULL)
      PD.min.cv <- min(PD.cv$Err.cv[which(PD.cv$grp=="PD")])

      print(sprintf("rfcv (cv=%.2f | PD.cv=%.2f | mk=%3d | tried %3s)",min.cv,PD.min.cv,marker.num,tried))
    }else{
      print(sprintf("rfcv (cv=%.2f | mk=%3d | tried %3s)",min.cv,marker.num,tried))
    }
    
    if(min.cv < tmp.min.cv){
      min.train.cv <- train.cv
      tmp.min.cv <- min.cv
    }
  }
  stopCluster(cl)
  #step 2
  marker.t <- table(unlist(lapply(train.cv, function(x) {
    lapply(x$res, "[", 1:marker.num)
  })))
  marker.t <- sort(marker.t, d = T)
  
  names(marker.t) <- colnames(tt$train.x)[as.numeric(names(marker.t))]
  marker.p <- names(marker.t)[1:marker.num]
  if(marker.num > 1){
    marker.anno <- as.data.frame(marker.t[1:marker.num])
  }else{
    marker.anno <- data.frame(id=names(marker.t[1:marker.num]),
                              marker.t[1:marker.num])
  }
  if(exists('mgs.V2.anno2')){
    colnames(marker.anno) <- c("mgs.V2","Freq")
    marker.anno <- merge(mgs.V2.anno2,marker.anno,by="mgs.V2",all.y=T)
  }else{
    colnames(marker.anno) <- c("ID","Freq")
  }
  
  
  rf.dat <- cbind(tt$train.x[, marker.p,F],grp=tt$train.y)
  train.rf <- randomForest(grp~.,data=rf.dat, importance = T,proximity=T,ntree=ntree)
  for(i in 2:rrf){
    tmp.rf <- randomForest(grp~.,data=rf.dat, importance = T,proximity=T,ntree=ntree)
    if( mean(tmp.rf$err.rate[,1]) > mean(train.rf$err.rate[,1]) ){
      train.rf <- tmp.rf
    }
  }
  train.p <- predict(train.rf, type = "prob")
  train.dat <- data.frame(grp=tt$train.y,prob=train.p[,2])
  #
  
  return(list(set=tt,cv=train.cv,mf=mmode,sf=smode,marker.num=marker.num,
              marker=marker.anno,rf=train.rf,train.pred=train.dat))
}

#for more than 2 groups
train.module3 <- function(prf,phe,ID,grp,size,
                          max.marker=30,max.cv=0.7,thread=5,rept=5,steps=2,cv.fold=5,
                          sam=F,rrf=10,smote=F,maxtry=10,impGrp=NULL,trees=500,
                          block=NULL,candy=NULL,pick=NULL,dw=NULL,ex.PID=NULL,...){
  #step1
  PD.min.cv <- 1
  min.cv <- 1
  marker.num=0
  tmp.min.cv <- 1
  ntree = trees
  # parallel method to speed up the process
  cl <- makeCluster(thread,...)
  clusterEvalQ(cl,library(randomForest))
  clusterExport(cl,c("tt","rfcv1"))
  tried <- 0
  while((marker.num>max.marker|marker.num==0|min.cv > max.cv) && tried <= maxtry){
    #aveSam(prf,phe,ID,grp,block=NULL,candy=NULL,size,pick=NULL,dw=NULL,ex.ATB=F)
    tried = tried + 1
    tt <- aveSam(prf,phe,ID,grp,NULL,candy,
                 size,pick,ex.PID=ex.PID,sam=sam,smote=smote)

    tt$train.y <-droplevels(tt$train.y)
    
    if(tried > maxtry){
      train.cv <- min.train.cv
      tried <- "end"
    }else{
      train.cv <- parSapply(cl, 0:rept, function(i,...) {
        set.seed(i)
        rfcv1(tt$train.x, tt$train.y, cv.fold = cv.fold, step = steps, impGrp = impGrp, ntree = trees,...)
      },simplify = F)
    }
    #
    if(is.null(impGrp)){
      err.cv <- t(matrix(unlist(lapply(train.cv,function(x){
        lv <- levels(tt$train.y)
        res <- sapply(x$predicted,function(y) {
          #mlv <- unlist(lapply(lv,FUN=function(x) mean((y!=tt$train.y)[which(tt$train.y==x)])))
          #res <- c(mlv,mean(mlv))
          #names(res) <- c(lv,"mean")
          res <- mean(unlist(lapply(lv,FUN=function(x) mean((y!=tt$train.y)[which(tt$train.y==x)]))))
          return(res)
        })
        return(res)
      })),ncol=length(train.cv[[1]]$n.var),byrow=T))
      rownames(err.cv) <- names(train.cv[[1]]$predicted)
      error.cv <- err.cv
    }else{
      imp.err.cv <- t(matrix(unlist(lapply(train.cv,function(x){
        sapply(x$predicted,function(y) {
          mean((y!=tt$train.y)[which(tt$train.y==impGrp)])
        })
      })),ncol=length(train.cv[[1]]$n.var),byrow=T))
      rownames(imp.err.cv) <- names(train.cv[[1]]$predicted)
      error.cv <- imp.err.cv
    }
    
    error.cv.rm <- rowMeans(error.cv)
    min.cv <- min(error.cv.rm)
    id <- which(error.cv.rm == min(error.cv.rm))
    
    marker.num <- min(as.numeric(names(error.cv.rm)[id]))
    
    mode <- data.frame(num=train.cv[[1]]$n.var,error.cv)
    mmode <- melt(mode,id.vars = "num",variable.name = "times",value.name = "Err.cv")
    mmode$num  <- as.factor(mmode$num)
    
    smode <- data.frame(num=train.cv[[1]]$n.var,times="rm",Err.cv=error.cv.rm)
    smode$num <- as.factor(smode$num)
    if(length(which(grepl("PD",levels(tt$train.y))))>0&&length(levels(tt$train.y))==3){
      PD.cv <- reCalCv(train=NULL,train.cv=train.cv,grp="PD",train.mf=mmode,train.sf=smode,train.y=tt$train.y,rep=NULL)
      PD.min.cv <- min(PD.cv$Err.cv[which(PD.cv$grp=="PD")])
      SD.cv <- reCalCv(train=NULL,train.cv=train.cv,grp="SD",train.mf=mmode,train.sf=smode,train.y=tt$train.y,rep=NULL)
      SD.min.cv <- min(SD.cv$Err.cv[which(SD.cv$grp=="SD")])
      PR.cv <- reCalCv(train=NULL,train.cv=train.cv,grp="PR",train.mf=mmode,train.sf=smode,train.y=tt$train.y,rep=NULL)
      PR.min.cv <- min(PR.cv$Err.cv[which(PR.cv$grp=="PR")])
      
      print(sprintf("rfcv (cv=%.2f | PD(%.2f)|SD(%.2f)|PR(%.2f) | mk=%3d | tried %3s)",
                    min.cv,PD.min.cv,SD.min.cv,PR.min.cv,marker.num,tried))
    }else{
      print(sprintf("rfcv (cv=%.2f | mk=%3d | tried %3s)",min.cv,marker.num,tried))
    }
    
    if(min.cv < tmp.min.cv){
      min.train.cv <- train.cv
      tmp.min.cv <- min.cv
    }
  }
  stopCluster(cl)
  #step 2
  marker.t <- table(unlist(lapply(train.cv, function(x) {
    lapply(x$res, "[", 1:marker.num)
  })))
  marker.t <- sort(marker.t, d = T)
  
  names(marker.t) <- colnames(tt$train.x)[as.numeric(names(marker.t))]
  marker.p <- names(marker.t)[1:marker.num]
  if(marker.num > 1){
    marker.anno <- as.data.frame(marker.t[1:marker.num])
  }else{
    marker.anno <- data.frame(id=names(marker.t[1:marker.num]),
                              marker.t[1:marker.num])
  }
  if(exists('mgs.V2.anno2')){
    colnames(marker.anno) <- c("mgs.V2","Freq")
    marker.anno <- merge(mgs.V2.anno2,marker.anno,by="mgs.V2",all.y=T)
  }else{
    colnames(marker.anno) <- c("ID","Freq")
  }
  
  
  rf.dat <- cbind(tt$train.x[, marker.p,F],grp=tt$train.y)
  train.rf <- randomForest(grp~.,data=rf.dat, importance = T,proximity=T,...)
  for(i in 2:rrf){
    tmp.rf <- randomForest(grp~.,data=rf.dat, importance = T,proximity=T,...)
    if( mean(tmp.rf$err.rate[,1]) < mean(train.rf$err.rate[,1]) ){
      train.rf <- tmp.rf
    }
  }
  train.p <- predict(train.rf, type = "prob")
  train.dat <- data.frame(grp=tt$train.y,train.p)
  #
  
  return(list(set=tt,cv=train.cv,mf=mmode,sf=smode,marker.num=marker.num,
              marker=marker.anno,rf=train.rf,train.pred=train.dat))
}
