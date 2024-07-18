## atakrig
## Function: area-to-point/area Kriging
## Author: Maogui Hu.

# require(sp)
# require(gstat)
# require(FNN)


## ataKriging: Area-to-area ordinary Kriging ----
# Input:
#   x: discretized area, list(areaValues, discretePoints):
#       areaValues: sample values, data.frame(areaId,centx,centy,value).
#       discretePoints: discretized area-samples, data.frame(areaId,ptx,pty,weight), weight is normalized.
#   unknown: discreteArea object or data.frame(areaId,ptx,pty,weight), discretized destination areas, weight is normalized.
#   ptVgm: point scale variogram, ataKrigVgm.
#   nmax: max number of neighborhoods used for interpolation.
#	  longlat: coordinates are longitude/latitude or not.
#   showProgress: show progress bar for batch interpolation (multi destination areas).
#   nopar: for internal use. Disable parallel process in the function even if ataStartCluster() has been called.
# Output: estimated value of destination area and its variance
ataKriging <- function(x, unknown, ptVgm, nmax=10, longlat=FALSE, showProgress=FALSE, nopar=FALSE, clarkAntiLog = FALSE) {
  stopifnot(nmax > 0)
  if(nmax < Inf) { # local neigbourhood Kriging.
    return(ataKriging.local(x, unknown, ptVgm, nmax, longlat, showProgress, nopar, clarkAntiLog))
  }

  if(is(unknown, "discreteArea")) unknown <- unknown$discretePoints
  if(is(ptVgm, "ataKrigVgm")) ptVgm <- extractPointVgm(ptVgm)

  sampleIds <- x$areaValues[,1]
  nSamples <- length(sampleIds)		# number of samples

  ## Kriging system: C * wmu = D
  # C matrix
  C <- matrix(1, nrow=nSamples+1, ncol=nSamples+1)
  sampleIndex <- list()
  for(i in 1:nSamples) {
    if(length(sampleIndex) < i) sampleIndex[[i]] <- x$discretePoints[,1] == sampleIds[i]
    sampleI <- x$discretePoints[sampleIndex[[i]],]
    for(j in i:nSamples) {
      if(length(sampleIndex) < j) sampleIndex[[j]] <- x$discretePoints[,1] == sampleIds[j]
      sampleJ <- x$discretePoints[sampleIndex[[j]],]
      C[i,j] <- ataCov(sampleI[,2:4], sampleJ[,2:4], ptVgm, longlat = longlat)
      C[j,i] <- C[i,j]
    }
  }
  C[nSamples+1,nSamples+1] <- 0

  unknownAreaIds <- sort(unique(unknown[,1]))

  krigOnce <- function(k) {
    cur <- unknown[unknown[,1] == unknownAreaIds[k], 2:4, drop=FALSE]

    # D matrix
    D <- matrix(1, nrow=nSamples+1, ncol=1)
    for(i in 1:nSamples) {
      sampleI <- x$discretePoints[sampleIndex[[i]],]
      D[i] <- ataCov(sampleI[,2:4,drop=FALSE], cur, ptVgm, longlat = longlat)
    }

    # solving
    solvedByGInv <- FALSE
    wmu <- try(solve(C, D), TRUE)
    if(is(wmu, "try-error")) {
      wmu <- MASS::ginv(C) %*% D
      solvedByGInv <- TRUE
    }
    w <- wmu[1:nSamples]
    mu <- wmu[(nSamples+1):nrow(wmu)]

    # estimation
    yest <- sum(w*x$areaValues[,4])
    yvar <- ataCov(cur, cur, ptVgm, longlat = longlat) - sum(wmu * D)

    if(!clarkAntiLog)
      return(data.frame(areaId=unknownAreaIds[k], pred=yest, var=yvar))
    else
      return(data.frame(areaId=unknownAreaIds[k], pred=yest, var=yvar,
                        pred.Clark=yest + sum(wmu * D)))
  }

  hasCluster <- ataIsClusterEnabled()
  if(showProgress) pb <- txtProgressBar(min=0, max=length(unknownAreaIds), width = 50, style = 3)

  if(!hasCluster || nopar || length(unknownAreaIds) == 1) {
    estResults <- c()
    for (k in 1:length(unknownAreaIds)) {
      estResults <- rbind(estResults, krigOnce(k))
      if(showProgress) setTxtProgressBar(pb, k)
    }
  } else {
    progress <- function(k) if(showProgress) setTxtProgressBar(pb, k)
    estResults <-
      foreach(k = 1:length(unknownAreaIds), .combine = rbind, .options.snow=list(progress=progress),
              .export = c("ataCov","calcAreaCentroid"),
              .packages = c("sp","gstat")) %dopar% {
                krigOnce(k)
              }
    ataClusterClearObj()
  }

  if(showProgress) close(pb)

  unknownCenter <- calcAreaCentroid(unknown)
  estResults <- merge(unknownCenter, estResults)

  return(estResults)
}


## ataKriging.cv: ataKriging cross-validation ----
#   nfold: integer; n-fold cross validation.
ataKriging.cv <- function(x, nfold=10, ptVgm, nmax=10, longlat=FALSE, showProgress=FALSE, nopar=FALSE, clarkAntiLog=FALSE) {
  N <- nrow(x$areaValues)
  if(missing(nfold)) {
    nfold <- N
  }

  if(length(nfold) == 1) {
    if(nfold <= 1 || nfold > N) {
      nfold <- N
    }

    if(nfold == N) {
      # leave-one-out
      indexM <- matrix(sort(x$areaValues[,1]), ncol = 1)
    } else {
      # n fold
      rndIds <- sample(x$areaValues[,1], N)
      nsize <- ceiling(N/nfold)
      indexM <- matrix(NA, nrow=nfold, ncol=nsize)
      nfrom <- nto <- 0
      for (i in 1:nfold) {
        nfrom <- nto + 1
        nto <- min(nsize * i, N)
        indexM[i,1:(nto-nfrom+1)] <- sort(rndIds[nfrom:nto])
      }
    }
  } else {
    indexM <- matrix(nfold, nrow = 1)
  }

  hasCluster <- ataIsClusterEnabled()
  if(showProgress) pb <- txtProgressBar(min=0, max=nrow(indexM), width = 50, style = 3)


  if(!hasCluster || nopar || nrow(indexM) == 1) {
    estResults <- c()
    for (k in 1:nrow(indexM)) {
      xknown <- subsetDiscreteArea(x, indexM[k,], revSel = TRUE)
      unknown <- subsetDiscreteArea(x, indexM[k,])$discretePoints
      estResults <- rbind(estResults, ataKriging(xknown, unknown, ptVgm, nmax, longlat, showProgress = FALSE, nopar = TRUE, clarkAntiLog))
      if(showProgress) setTxtProgressBar(pb, k)
    }
  } else {
    bInnerParallel <- ncol(indexM) > 2*nrow(indexM)
    if(bInnerParallel) {
      estResults <- c()
      for (k in 1:nrow(indexM)) {
        xknown <- subsetDiscreteArea(x, indexM[k,], revSel = TRUE)
        unknown <- subsetDiscreteArea(x, indexM[k,])$discretePoints
        estResults <- rbind(estResults, ataKriging(xknown, unknown, ptVgm, nmax, longlat, showProgress = FALSE, nopar = FALSE, clarkAntiLog))
        if(showProgress) setTxtProgressBar(pb, k)
      }
    } else {
      progress <- function(k) if(showProgress) setTxtProgressBar(pb, k)
      estResults <-
        foreach(k = 1:nrow(indexM), .combine = rbind, .options.snow=list(progress=progress),
                .export = c("subsetDiscreteArea","ataCov","calcAreaCentroid","ataKriging","ataKriging.local"),
                .packages = c("sp","gstat","FNN")) %dopar% {
                  xknown <- subsetDiscreteArea(x, indexM[k,], revSel = TRUE)
                  unknown <- subsetDiscreteArea(x, indexM[k,])$discretePoints
                  ataKriging(xknown, unknown, ptVgm, nmax, longlat, showProgress = FALSE, nopar = TRUE, clarkAntiLog)
                }
      ataClusterClearObj()
    }
  }
  if(showProgress) close(pb)

  estResults <- estResults[order(estResults$areaId),]
  indx <- match(estResults[,1], x$areaValues[,1])
  # estResults$diff <- x$areaValues[indx,4] - estResults[,4]
  estResults$value <- x$areaValues[indx,4]

  return(estResults)
}


## ataKriging.local: [internal use only]. ----
ataKriging.local <- function(x, unknown, ptVgm, nmax=10, longlat=FALSE, showProgress=FALSE, nopar=FALSE, clarkAntiLog=FALSE) {
  if(is(unknown, "discreteArea")) unknown <- unknown$discretePoints

  unknown <- unknown[sort.int(unknown[,1], index.return = TRUE)$ix,]
  unknownCenter <- calcAreaCentroid(unknown)

  nb <- FNN::get.knnx(as.matrix(x$areaValues[,2:3,drop=FALSE]), as.matrix(unknownCenter[,2:3,drop=FALSE]), nmax)
  nb$nn.index <- matrix(x$areaValues[nb$nn.index,1], ncol = ncol(nb$nn.index))

  unknownAreaIds <- unknownCenter[,1]

  krigOnce <- function(k) {
    curUnknown <- unknown[unknown[,1] == unknownAreaIds[k],]
    curAreaPts <- x$discretePoints[x$discretePoints[,1] %in% nb$nn.index[k,],]
    curAreaVals <- x$areaValues[x$areaValues[,1] %in% nb$nn.index[k,],]

    estResult <- ataKriging(x=list(discretePoints=curAreaPts, areaValues = curAreaVals),
                            curUnknown, ptVgm, nmax=Inf, longlat=longlat, showProgress=FALSE, nopar=TRUE, clarkAntiLog)
    return(estResult)
  }

  hasCluster <- ataIsClusterEnabled()
  if(showProgress) pb <- txtProgressBar(min=0, max=length(unknownAreaIds), width = 50, style = 3)

  if(!hasCluster || nopar || length(unknownAreaIds) == 1) {
    estResults <- c()
    for (k in 1:length(unknownAreaIds)) {
      estResults <- rbind(estResults, krigOnce(k))
      if(showProgress) setTxtProgressBar(pb, k)
    }
  } else {
    progress <- function(k) if(showProgress) setTxtProgressBar(pb, k)
    estResults <-
      foreach(k = 1:length(unknownAreaIds), .combine = rbind, .options.snow=list(progress=progress),
              .export = c("ataKriging","ataCov","calcAreaCentroid","extractPointVgm"),
              .packages = c("sp","gstat")) %dopar% {
                krigOnce(k)
              }
    ataClusterClearObj()
  }
  if(showProgress) close(pb)

  return(estResults)
}


## ataCoKriging: Area-to-area ordinary CoKriging ----
# Input:
#   x: discretized areas, list(
#      `varId1`=list(areaValues=data.frame(areaId,centx,centy,value), discretePoints=data.frame(areaId,ptx,pty,weight)),
#      `varId2`=list(areaValues=data.frame(areaId,centx,centy,value), discretePoints=data.frame(areaId,ptx,pty,weight)),
#       ...,
#      `varIdn`=list(areaValues=data.frame(areaId,centx,centy,value), discretePoints=data.frame(areaId,ptx,pty,weight)))
#   unknownVarId: variable name (charaster) defined in x for prediction.
#   unknown: discretized destination area, data.frame(areaId,ptx,pty,weight).
#   ptVgms: point scale direct or cross variograms, ataKrigVgm.
#   nmax: max number of neighborhoods used for interpolation.
#	  longlat: coordinates are longitude/latitude or not.
#   oneCondition: oneCondition coKriging, assuming expected means of variables known and constant with the study area.
#   meanVal: expected means of variables for oneCondition coKriging, data.frame(varId,value). If missing, simple mean values of
#      areas from x will be used instead.
#   auxRatioAdj: for oneCondition Kriging, adjusting the auxiliary variable residue by a ratio between the primary variable mean
#      and auxiliary variable mean.
#   showProgress: show progress bar for batch interpolation (multi destination areas).
#   nopar: disable parallel process in the function even if ataStartCluster() has been called, mainly for  internal use.
# Output: estimated value of destination area and its variance
ataCoKriging_backup <- function(x, unknownVarId, unknown, ptVgms, nmax=10, longlat=FALSE, oneCondition=FALSE,
                         meanVal=NULL, auxRatioAdj=TRUE, showProgress=FALSE, nopar=FALSE, clarkAntiLog=FALSE) {
  
    
  stopifnot(nmax > 0)
  
  if(nmax < Inf) {
    cat("usign ataCoKriging local\n")
    return(ataCoKriging.local(x, unknownVarId, unknown, ptVgms, nmax, longlat, oneCondition,
                              meanVal, auxRatioAdj, showProgress, nopar, clarkAntiLog))
  }
  
  if(is(unknown, "discreteArea")) unknown <- unknown$discretePoints
  if(is(ptVgms, "ataKrigVgm")) ptVgms <- extractPointVgm(ptVgms)
  
  # sort areaId in ascending order.
  for (i in 1:length(x)) {
    x[[i]]$areaValues <- x[[i]]$areaValues[sort.int(x[[i]]$areaValues[,1], index.return = TRUE)$ix,]
  }
  
  # combine all data together.
  varIds <- sort(names(x))
  xAll <- list(areaValues=NULL, discretePoints=NULL)
  for (id in varIds) {
    if(!hasName(x[[id]], "discretePoints")) {
      x[[id]]$discretePoints <- cbind(x[[id]]$areaValues[,1:3], data.frame(weight=rep(1,nrow(x[[id]]$areaValues))))
      names(x[[id]]$discretePoints)[2:3] <- c("ptx","pty")
    }
    
    x[[id]]$areaValues$varId <- id
    x[[id]]$areaValues$var_areaId <- paste(id, x[[id]]$areaValues[,1], sep = "_")
    x[[id]]$discretePoints$varId <- id
    x[[id]]$discretePoints$var_areaId <- paste(id, x[[id]]$discretePoints[,1], sep = "_")
    
    xAll$areaValues <- rbind(xAll$areaValues, x[[id]]$areaValues)
    xAll$discretePoints <- rbind(xAll$discretePoints, x[[id]]$discretePoints)
  }
  
  sampleIds <- sort(unique(xAll$discretePoints$var_areaId))
  nSamples <- length(sampleIds)		# number of all samples
  nVars <- length(x) # number of variables
  
  ## Kriging system: C * wmu = D
  if(oneCondition) {
    C <- matrix(0, nrow=nSamples+1, ncol=nSamples+1)
    D <- matrix(0, nrow=nSamples+1, ncol=1)
  } else {
    C <- matrix(0, nrow=nSamples+nVars, ncol=nSamples+nVars)
    D <- matrix(0, nrow=nSamples+nVars, ncol=1)
  }
  
  # C matrix
  sampleIndex <- list()
  for(i in 1:nSamples) {
    if(length(sampleIndex) < i) sampleIndex[[i]] <- xAll$discretePoints$var_areaId == sampleIds[i]
    sampleI <- xAll$discretePoints[sampleIndex[[i]],]
    for(j in i:nSamples) {
      if(length(sampleIndex) < j) sampleIndex[[j]] <- xAll$discretePoints$var_areaId == sampleIds[j]
      sampleJ <- xAll$discretePoints[sampleIndex[[j]],]
      ptVgm <- ptVgms[[crossName(sampleI$varId[1], sampleJ$varId[1])]]
      C[i,j] <- ataCov(sampleI[,2:4], sampleJ[,2:4], ptVgm, longlat = longlat)
      C[j,i] <- C[i,j]
    }
  }
  
  if(oneCondition) {
    C[nSamples+1, ] <- 1
    C[, nSamples+1] <- 1
    C[nSamples+1, nSamples+1] <- 0
    D[nSamples+1] <- 1
  } else {
    for (i in 1:nVars) {
      indx <- xAll$areaValues$varId == varIds[i]
      C[nSamples+i, (1:nSamples)[indx]] <- 1
      C[(1:nSamples)[indx], nSamples+i] <- 1
    }
    D[nSamples + which(unknownVarId == varIds)] <- 1
  }
  
  
  unknownAreaIds <- sort(unique(unknown[,1]))
  
  krigOnce <- function(k) {
    curUnknown <- unknown[unknown[,1] == unknownAreaIds[k], 2:4]
    
    # D matrix
    for(i in 1:nSamples) {
      sampleI <- xAll$discretePoints[xAll$discretePoints$var_areaId == sampleIds[i],]
      ptVgm <- ptVgms[[crossName(sampleI$varId[1], unknownVarId)]]
      D[i] <- ataCov(sampleI[,2:4], curUnknown, ptVgm, longlat = longlat)
    }

    # Soglie di mal condizionamento
    cond_number_threshold <- 1e4
    determinant_threshold <- 1e-10
    eigenvalue_threshold <- 1e-10
  
    # Calcoli
    kappa <- calc_condition_number(C)
    det_C <- calc_determinant(C)
    eigenvalues <- calc_eigenvalues(C)
    
    # Verifica delle condizioni
    if (kappa > cond_number_threshold) {
      # cat("La matrice è mal condizionata: Numero di condizionamento elevato\n")
      # cat("Numero di condizionamento:", kappa, "\n")
      C <- C + 1e-5 * diag(nrow(C))
      kappa2 <- calc_condition_number(C)
      
      # Verifica delle condizioni
      if (kappa2 > cond_number_threshold) {
        # cat("La matrice è ancora mal condizionata: Numero di condizionamento elevato\n")
        # cat("Numero di condizionamento:", kappa, "\n")
        # cat("Numero di condizionamento post normalizzazione:", kappa2, "\n")
      }
    }
  
    if (abs(det_C) < determinant_threshold) {
      cat("La matrice è mal condizionata: Determinante molto piccolo\n")
      cat("Determinante della matrice:", det_C, "\n")
    }
  
    if (any(abs(eigenvalues) < eigenvalue_threshold)) {
      cat("La matrice è mal condizionata: Autovalori molto piccoli\n")
      cat("Autovalori della matrice:", eigenvalues, "\n")
    }

    
    # solving
    # cat("Solving system for unknown area ID: ", unknownAreaIds[k], "\n")
    # solvedByGInv <- FALSE
    # wmu <- try(solve(C, D), FALSE)
    # if(is(wmu, "try-error")) {
    #   cat("Solving by GInv\n")
    #   wmu <- MASS::ginv(C) %*% D
    #   solvedByGInv <- TRUE
    # }
    wmu <- MASS::ginv(C) %*% D
    # wmu <- solve_via_svd(C,D)
    
    # estimation
    if(oneCondition) {
      if(is.null(meanVal)) {
        for (id in varIds) {
          meanVal <- rbind(meanVal, data.frame(varId=id, value=mean(x[[id]]$areaValues[,4])))
        }
      }
      rownames(meanVal) <- meanVal$varId
      
      w <- wmu[1:nSamples]
      w1 <- w[unknownVarId == xAll$areaValues$varId]
      yest <- sum(w1 * x[[unknownVarId]]$areaValues[,4])
      for (id in varIds[varIds != unknownVarId]) {
        w2 <- w[id == xAll$areaValues$varId]
        if (auxRatioAdj && abs(meanVal[id,2]) > 1e-6) {
          yest <- yest + sum(w2 * ((x[[id]]$areaValues[,4] - meanVal[id, 2])*(meanVal[unknownVarId, 2]/meanVal[id,2]) + meanVal[unknownVarId, 2]))
        } else {
          yest <- yest + sum(w2 * (x[[id]]$areaValues[,4] - meanVal[id, 2] + meanVal[unknownVarId, 2]))
        }
      }
    } else {
      w <- wmu[1:nSamples][unknownVarId == xAll$areaValues$varId]
      yest <- sum(w * x[[unknownVarId]]$areaValues[,4])
    }
    yvar <- ataCov(curUnknown, curUnknown, ptVgms[[unknownVarId]], longlat = longlat) - sum(wmu * D)
    
    if(!clarkAntiLog)
      return(data.frame(areaId=unknownAreaIds[k], pred=yest, var=yvar))
    else
      return(data.frame(areaId=unknownAreaIds[k], pred=yest, var=yvar,
                        pred.Clark=yest + sum(wmu * D)))
  }
  
  hasCluster <- ataIsClusterEnabled()
  if(showProgress) pb <- txtProgressBar(min=0, max=length(unknownAreaIds), width = 50, style = 3)
  
  if(!hasCluster || nopar || length(unknownAreaIds) == 1) {
    estResults <- c()
    for (k in 1:length(unknownAreaIds)) {
      estResults <- rbind(estResults, krigOnce(k))
      if(showProgress) setTxtProgressBar(pb, k)
    }
  } else {
    progress <- function(k) if(showProgress) setTxtProgressBar(pb, k)
    estResults <-
      foreach(k = 1:length(unknownAreaIds), .combine = rbind, .options.snow=list(progress=progress),
              .export = c("D","meanVal","crossName","ataCov","calcAreaCentroid", "spDistsNN"),
              .packages = c("sp","gstat")) %dopar% {
                krigOnce(k)
              }
    ataClusterClearObj()
  }
  
  if(showProgress) close(pb)
  
  unknownCenter <- calcAreaCentroid(unknown)
  estResults <- merge(unknownCenter, estResults)

  # Stampa i valori anomali
  if(abs(estResults$pred) > 1000){
    cat("Valore di pred:", estResults$pred, "\n")
  }
  
  return(estResults)
}


## ataCoKriging.cv: ataCoKriging cross validation. ----
#   nfold: integer; n-fold cross validation.
ataCoKriging.cv <- function(x, unknownVarId, nfold=10, ptVgms, nmax=10, longlat=FALSE, oneCondition=FALSE,
                            meanVal=NULL, auxRatioAdj=TRUE, showProgress=FALSE, nopar=FALSE, clarkAntiLog=FALSE) {
  N <- nrow(x[[unknownVarId]]$areaValues)
  if(missing(nfold)) {
    nfold <- N
  }

  if(length(nfold) == 1) {
    if(nfold <= 1 || nfold > N) {
      nfold <- N
    }

    if(nfold == N) {
      # leave-one-out
      indexM <- matrix(1:N, ncol = 1)
    } else {
      # n fold
      rndIds <- sample(x[[unknownVarId]]$areaValues[,1], N)
      nsize <- ceiling(N/nfold)
      indexM <- matrix(NA, nrow=nfold, ncol=nsize)
      nfrom <- nto <- 0
      for (i in 1:nfold) {
        nfrom <- nto + 1
        nto <- min(nsize * i, N)
        indexM[i,1:(nto-nfrom+1)] <- sort(rndIds[nfrom:nto])
      }
    }
  } else {
    indexM <- matrix(nfold, nrow = 1)
  }

  hasCluster <- ataIsClusterEnabled()
  if(showProgress) pb <- txtProgressBar(min=0, max=nrow(indexM), width = 50, style = 3)

  xknown <- x

  if(!hasCluster || nopar || nrow(indexM) == 1) {
    estResults <- c()
    for (k in 1:nrow(indexM)) {
      xknown[[unknownVarId]] <- subsetDiscreteArea(x[[unknownVarId]], indexM[k,], revSel = TRUE)
      unknown <- subsetDiscreteArea(x[[unknownVarId]], indexM[k,])$discretePoints
      estResults <- rbind(estResults,
                          ataCoKriging(xknown, unknownVarId, unknown, ptVgms, nmax, longlat, oneCondition,
                                       meanVal, auxRatioAdj, showProgress = FALSE, nopar = TRUE, clarkAntiLog))
      if(showProgress) setTxtProgressBar(pb, k)
    }
  } else {
    bInnerParallel <- ncol(indexM) > 2*nrow(indexM)
    if(bInnerParallel) {
      estResults <- c()
      for (k in 1:nrow(indexM)) {
        xknown[[unknownVarId]] <- subsetDiscreteArea(x[[unknownVarId]], indexM[k,], revSel = TRUE)
        unknown <- subsetDiscreteArea(x[[unknownVarId]], indexM[k,])$discretePoints
        estResults <- rbind(estResults,
                            ataCoKriging(xknown, unknownVarId, unknown, ptVgms, nmax, longlat, oneCondition,
                                         meanVal, auxRatioAdj, showProgress = FALSE, nopar = FALSE, clarkAntiLog))
        if(showProgress) setTxtProgressBar(pb, k)
      }
    } else {
      progress <- function(k) if(showProgress) setTxtProgressBar(pb, k)
      estResults <-
        foreach(k = 1:nrow(indexM), .combine = rbind, .options.snow=list(progress=progress),
                .export = c("crossName","ataCov","calcAreaCentroid","subsetDiscreteArea","ataCoKriging","ataCoKriging.local"),
                .packages = c("sp","gstat","FNN")) %dopar% {
                  xknown[[unknownVarId]] <- subsetDiscreteArea(x[[unknownVarId]], indexM[k,], revSel = TRUE)
                  unknown <- subsetDiscreteArea(x[[unknownVarId]], indexM[k,])$discretePoints
                  ataCoKriging(xknown, unknownVarId, unknown, ptVgms, nmax, longlat, oneCondition, meanVal,
                               auxRatioAdj, showProgress = FALSE, nopar = TRUE, clarkAntiLog)
                }
      ataClusterClearObj()
    }
  }
  if(showProgress) close(pb)

  estResults <- estResults[order(estResults$areaId),]
  indx <- match(estResults[,1], x[[unknownVarId]]$areaValues[,1])
  # estResults$diff <- x[[unknownVarId]]$areaValues[indx,4] - estResults[,4]
  estResults$value <- x[[unknownVarId]]$areaValues[indx,4]

  return(estResults)
}


## ataCoKriging.local: [internal use only]. ----
ataCoKriging.local_backup <- function(x, unknownVarId, unknown, ptVgms, nmax=10, longlat=FALSE,
                               oneCondition=FALSE, meanVal=NULL, auxRatioAdj=TRUE,
                               showProgress=FALSE, nopar=FALSE, clarkAntiLog=FALSE) {

  
  if(is(unknown, "discreteArea")) {
    unknown <- unknown$discretePoints
    cat("Convertito 'unknown' da 'discreteArea' a 'discretePoints'\n")
  }
  
  # sort areaId in ascending order.
  unknown <- unknown[sort.int(unknown[,1], index.return = TRUE)$ix,]
  cat("Ordinato 'unknown' per 'areaId'\n")
  unknownCenter <- calcAreaCentroid(unknown)
  cat("Calcolato il centroide di 'unknown'\n")
  
  # neighbor indexes for each unknown point.
  varIds <- sort(names(x))
  nb <- list()
  for (id in varIds) {
    cat(sprintf("Calcolo dei vicini per la variabile %s\n", id))
    nb[[id]] <- FNN::get.knnx(as.matrix(x[[id]]$areaValues[,2:3,drop=FALSE]), as.matrix(unknownCenter[,2:3,drop=FALSE]), nmax)
    nb[[id]]$nn.index <- matrix(x[[id]]$areaValues[,1][nb[[id]]$nn.index], ncol = nmax)
  }
  
  # only consider covariables within the radius of unknownVarId
  for (id in varIds[varIds != unknownVarId]) {
    cat(sprintf("Filtraggio dei vicini per la variabile %s\n", id))
    indx <- nb[[id]]$nn.dist > matrix(rep(nb[[unknownVarId]]$nn.dist[,nmax] * 1.5, nmax), ncol = nmax)
    nb[[id]]$nn.dist[indx] <- NA
    nb[[id]]$nn.index[indx] <- NA
  }
  
  unknownAreaIds <- sort(unique(unknown[,1]))
  cat("Calcolati gli ID delle aree sconosciute\n")
  
  krigOnce <- function(k) {
    curUnknown <- unknown[unknown[,1] == unknownAreaIds[k], ]
    
    curx <- list()
    for (id in varIds) {
      if(!hasName(x[[id]], "discretePoints")) {
        x[[id]]$discretePoints <- cbind(x[[id]]$areaValues[,1:3], data.frame(weight=rep(1,nrow(x[[id]]$areaValues))))
        names(x[[id]]$discretePoints)[2:3] <- c("ptx","pty")
      }
      
      curVals <- x[[id]]$areaValues[x[[id]]$areaValues[,1] %in% nb[[id]]$nn.index[k,],]
      curPts <- x[[id]]$discretePoints[x[[id]]$discretePoints[,1] %in% nb[[id]]$nn.index[k,],]
      if(nrow(curVals) > 0) {
        curx[[id]] <- list(areaValues=curVals, discretePoints=curPts)
      }
    }
    
    estResult <- ataCoKriging(curx, unknownVarId, curUnknown, ptVgms, nmax=Inf, longlat, oneCondition,
                              meanVal, auxRatioAdj, showProgress=FALSE, nopar=TRUE, clarkAntiLog)
    return(estResult)
  }
  
  hasCluster <- ataIsClusterEnabled()
  if(showProgress) pb <- txtProgressBar(min=0, max=length(unknownAreaIds), width = 50, style = 3)
  
  if(!hasCluster || nopar || length(unknownAreaIds) == 1) {
    estResults <- c()
    for (k in 1:length(unknownAreaIds)) {
      estResults <- rbind(estResults, krigOnce(k))
      if(showProgress) setTxtProgressBar(pb, k)
    }
  } else {
    cat("Inizio stima parallela\n")
    progress <- function(k) if(showProgress) setTxtProgressBar(pb, k)
    estResults <-
      foreach(k = 1:length(unknownAreaIds), .combine = rbind, .options.snow=list(progress=progress),
              .export = c("x","ataCoKriging","crossName","ataCov","calcAreaCentroid", "spDistsNN"),
              .packages = c("sp","gstat")) %dopar% {
                krigOnce(k)
              }
    ataClusterClearObj()
  }
  if(showProgress) close(pb)
  
  cat("Fine funzione ataCoKriging.local\n")
  return(estResults)
}


## ataCov: [internal use only] Covariance between two discretized area-samples. ----
# Input:
#   areaPts1: first discretized area-sample, data.frame(ptx,pty,weight), weight is normalized.
#   areaPts2: second discretized area-sample, data.frame(ptx,pty,weight), weight is normalized.
#   ptVgm: point scale variogram (gstat vgm).
#	  longlat: indicator whether coordinates are longitude/latitude
ataCov <- function(areaPts1, areaPts2, ptVgm, longlat=FALSE) {
  # disM <- spDists(as.matrix(areaPts1[,1:2,drop=FALSE]), as.matrix(areaPts2[,1:2,drop=FALSE]), longlat=longlat)
  # mCov <- variogramLine(ptVgm, covariance=TRUE, dist_vector=disM)
  # return(sum(outer(areaPts1[,3], areaPts2[,3]) * mCov))

  disM <- spDistsNN(areaPts1[,1], areaPts1[,2], areaPts2[,1], areaPts2[,2], longlat=longlat)
  mCov <- variogramLineSimple(ptVgm, disM, bCov = TRUE)
  return(sum(outerProd(areaPts1[,3], areaPts2[,3]) * mCov))
}


## atpKriging: Area-to-point ordinary Kriging ----
# Input:
#   x: discretized areas, list(discretePoints, areaValues):
#       areaValues: values of areas, data.frame(areaId,centx,centy,value).
#       discretePoints: discretized points of areas, data.frame(areaId,ptx,pty,weight), the weight is normalized.
#   unknown0: single discretized destination area, data.frame(ptx,pty).
#   ptVgm: point scale variogram, ataKrigVgm.
#   nmax: max number of neighborhoods used for interpolation.
#	  longlat: coordinates are longitude/latitude or not.
#   showProgress: show progress bar for batch interpolation (multi destination areas).
#   nopar: for internal use. Disable parallel process in the function even if ataStartCluster() has been called.
# Output: estimated value of destination area and its variance.
atpKriging <- function(x, unknown0, ptVgm, nmax=10, longlat=FALSE, showProgress=FALSE, nopar=FALSE) {
  unknown <- cbind(areaId=1:nrow(unknown0), unknown0, weight=1)
  return(ataKriging(x, unknown, ptVgm, nmax, longlat, showProgress, nopar))
}


## atpCoKriging: Area-to-point ordinary CoKriging ----
# Input:
#   x: discretized areas, list(
#      `varId1`=list(areaValues=data.frame(areaId,centx,centy,value), discretePoints=data.frame(areaId,ptx,pty,weight)),
#      `varId2`=list(areaValues=data.frame(areaId,centx,centy,value), discretePoints=data.frame(areaId,ptx,pty,weight)),
#       ...,
#      `varIdn`=list(areaValues=data.frame(areaId,centx,centy,value), discretePoints=data.frame(areaId,ptx,pty,weight)))
#   unknown0: unknown points, data.frame(ptx,pty).
#   ptVgms: point scale direct or cross variograms, ataKrigVgm.
#   nmax: max number of neighborhoods used for interpolation for the main variable.
#	  longlat: indicator whether coordinates are longitude/latitude.
#   oneCondition: use oneCondition coKriging, assuming expected means of variables known and constant with the study area.
#   meanVal: expected means of variables for oneCondition coKriging, data.frame(varId,value). If missing, simple mean values of
#      areas from x will be used instead.
#   showProgress: show progress bar for batch interpolation (multi destination areas).
#   nopar: for internal use. Disable parallel process in the function even if ataStartCluster() has been called.
# Output: estimated value of destination area and its variance.
atpCoKriging <- function(x, unknownVarId, unknown0, ptVgms, nmax=10, longlat=FALSE, oneCondition=FALSE,
                         meanVal=NULL, auxRatioAdj=TRUE, showProgress=FALSE, nopar=FALSE) {
  unknown <- cbind(areaId=1:nrow(unknown0), unknown0, weight=1)
  return(ataCoKriging(x, unknownVarId, unknown, ptVgms, nmax, longlat, oneCondition, meanVal, auxRatioAdj, showProgress, nopar))
}

# Calcolo del numero di condizionamento usando i valori singolari
calc_condition_number <- function(C) {
  svd_C <- svd(C)
  kappa <- max(svd_C$d) / min(svd_C$d)
  return(kappa)
}

# Calcolo del determinante
calc_determinant <- function(C) {
  det_C <- det(C)
  return(det_C)
}

# Verifica degli autovalori
calc_eigenvalues <- function(C) {
  eigen_C <- eigen(C)
  return(eigen_C$values)
}

ataCoKriging.local <- function(x, unknownVarId, unknown, ptVgms, nmax=10, longlat=FALSE,
                               oneCondition=FALSE, meanVal=NULL, auxRatioAdj=TRUE,
                               showProgress=FALSE, nopar=FALSE, clarkAntiLog=FALSE) {
  
  cat("Inizio funzione ataCoKriging.local\n")
  
  if (is(unknown, "discreteArea")) {
    unknown <- unknown$discretePoints
    cat("Convertito 'unknown' da 'discreteArea' a 'discretePoints'\n")
  }
  
  # Sort areaId in ascending order.
  unknown <- unknown[sort.int(unknown[, 1], index.return = TRUE)$ix, ]
  cat("Ordinato 'unknown' per 'areaId'\n")
  unknownCenter <- calcAreaCentroid(unknown)
  cat("Calcolato il centroide di 'unknown'\n")
  
  # Neighbor indexes for each unknown point.
  varIds <- sort(names(x))
  nb <- list()
  for (id in varIds) {
    cat(sprintf("Calcolo dei vicini per la variabile %s\n", id))
    nb[[id]] <- FNN::get.knnx(as.matrix(x[[id]]$areaValues[, 2:3, drop = FALSE]), as.matrix(unknownCenter[, 2:3, drop = FALSE]), nmax)
    nb[[id]]$nn.index <- matrix(x[[id]]$areaValues[, 1][nb[[id]]$nn.index], ncol = nmax)
  }
  
  # Only consider covariables within the radius of unknownVarId
  for (id in varIds[varIds != unknownVarId]) {
    cat(sprintf("Filtraggio dei vicini per la variabile %s\n", id))
    indx <- nb[[id]]$nn.dist > matrix(rep(nb[[unknownVarId]]$nn.dist[, nmax] * 1.5, nmax), ncol = nmax)
    nb[[id]]$nn.dist[indx] <- NA
    nb[[id]]$nn.index[indx] <- NA
  }
  
  unknownAreaIds <- sort(unique(unknown[, 1]))
  cat("Calcolati gli ID delle aree sconosciute\n")
  
  krigOnce <- function(k) {
    cat(sprintf("Inizio krigOnce per areaId %d\n", unknownAreaIds[k]))
    curUnknown <- unknown[unknown[, 1] == unknownAreaIds[k], ]
    
    if (nrow(curUnknown) == 0) {
      cat("Nessun dato per questa area. Restituzione NA\n")
      return(list(estimates = NA, C = NA, D = NA))
    }

    curx <- list()
    for (id in varIds) {
      if (!hasName(x[[id]], "discretePoints")) {
        x[[id]]$discretePoints <- cbind(x[[id]]$areaValues[, 1:3], data.frame(weight = rep(1, nrow(x[[id]]$areaValues))))
        names(x[[id]]$discretePoints)[2:3] <- c("ptx", "pty")
      }
      
      curVals <- x[[id]]$areaValues[x[[id]]$areaValues[, 1] %in% nb[[id]]$nn.index[k, ], ]
      curPts <- x[[id]]$discretePoints[x[[id]]$discretePoints[, 1] %in% nb[[id]]$nn.index[k, ], ]
      if (nrow(curVals) > 0) {
        curx[[id]] <- list(areaValues = curVals, discretePoints = curPts)
      }
    }

    if (length(curx) == 0) {
      cat("Nessun valore per le variabili. Restituzione NA\n")
      return(list(estimates = NA, C = NA, D = NA))
    }
    
    result <- ataCoKriging(curx, unknownVarId, curUnknown, ptVgms, nmax = Inf, longlat, oneCondition,
                           meanVal, auxRatioAdj, showProgress = FALSE, nopar = TRUE, clarkAntiLog)
    estResult <- result$estimates
    C <- result$C
    D <- result$D
    cat("Fine krigOnce\n")
    
    return(list(estimates = estResult, C = C, D = D))
  }
  
  hasCluster <- ataIsClusterEnabled()
  if (showProgress) pb <- txtProgressBar(min = 0, max = length(unknownAreaIds), width = 50, style = 3)
  cat("Verificata presenza di cluster\n")
  
  estResults <- list(estimates = NULL, C = NULL, D = NULL, highEstimateAreas = list())
  if (!hasCluster || nopar || length(unknownAreaIds) == 1) {
    for (k in 1:length(unknownAreaIds)) {
      cat(sprintf("Inizio kriging per areaId %d\n", unknownAreaIds[k]))
      res <- krigOnce(k)
      if (!is.null(res$estimates) && all(!is.na(res$estimates))) {
        estResults$estimates <- rbind(estResults$estimates, res$estimates)
        if (is.null(estResults$C) && is.null(estResults$D)) {
          estResults$C <- res$C
          estResults$D <- res$D
        }
        if (any(res$estimates$pred > 100000)) {
          cat(sprintf("Aggiunta areaId %d alla lista highEstimateAreas\n", unknownAreaIds[k]))
          estResults$highEstimateAreas[[as.character(unknownAreaIds[k])]] <- list(C = res$C, D = res$D)
        }
      }
      if (showProgress) setTxtProgressBar(pb, k)
    }
  } else {
    cat("Inizio stima parallela\n")
    progress <- function(k) if (showProgress) setTxtProgressBar(pb, k)
    results <-
      foreach(k = 1:length(unknownAreaIds), .combine = list, .multicombine = TRUE, .options.snow = list(progress = progress),
              .export = c("x", "ataCoKriging", "crossName", "ataCov", "calcAreaCentroid", "spDistsNN"),
              .packages = c("sp", "gstat")) %dopar% {
                krigOnce(k)
              }
    for (res in results) {
      if (!is.null(res$estimates) && all(!is.na(res$estimates))) {
        estResults$estimates <- rbind(estResults$estimates, res$estimates)
        if (is.null(estResults$C) && is.null(estResults$D)) {
          estResults$C <- res$C
          estResults$D <- res$D
        }
        if (any(res$estimates$pred > 100000)) {
          cat(sprintf("Aggiunta areaId %s alla lista highEstimateAreas\n", names(res$estimates)))
          estResults$highEstimateAreas[[as.character(names(res$estimates))]] <- list(C = res$C, D = res$D)
        }
      }
    }
    ataClusterClearObj()
  }
  if (showProgress) close(pb)
  
  cat("Fine funzione ataCoKriging.local\n")
  return(estResults)
}



ataCoKriging <- function(x, unknownVarId, unknown, ptVgms, nmax=10, longlat=FALSE, oneCondition=FALSE,
                                 meanVal=NULL, auxRatioAdj=TRUE, showProgress=FALSE, nopar=FALSE, clarkAntiLog=FALSE) {

  stopifnot(nmax > 0)

  if(nmax < Inf) {
    cat("usign ataCoKriging local\n")
    return(ataCoKriging.local(x, unknownVarId, unknown, ptVgms, nmax, longlat, oneCondition,
                                      meanVal, auxRatioAdj, showProgress, nopar, clarkAntiLog))
  }

  if(is(unknown, "discreteArea")) unknown <- unknown$discretePoints
  if(is(ptVgms, "ataKrigVgm")) ptVgms <- extractPointVgm(ptVgms)

  # sort areaId in ascending order.
  for (i in 1:length(x)) {
    x[[i]]$areaValues <- x[[i]]$areaValues[sort.int(x[[i]]$areaValues[,1], index.return = TRUE)$ix,]
  }

  # combine all data together.
  varIds <- sort(names(x))
  xAll <- list(areaValues=NULL, discretePoints=NULL)
  for (id in varIds) {
    if(!hasName(x[[id]], "discretePoints")) {
      x[[id]]$discretePoints <- cbind(x[[id]]$areaValues[,1:3], data.frame(weight=rep(1,nrow(x[[id]]$areaValues))))
      names(x[[id]]$discretePoints)[2:3] <- c("ptx","pty")
    }

    x[[id]]$areaValues$varId <- id
    x[[id]]$areaValues$var_areaId <- paste(id, x[[id]]$areaValues[,1], sep = "_")
    x[[id]]$discretePoints$varId <- id
    x[[id]]$discretePoints$var_areaId <- paste(id, x[[id]]$discretePoints[,1], sep = "_")

    xAll$areaValues <- rbind(xAll$areaValues, x[[id]]$areaValues)
    xAll$discretePoints <- rbind(xAll$discretePoints, x[[id]]$discretePoints)
  }

  sampleIds <- sort(unique(xAll$discretePoints$var_areaId))
  nSamples <- length(sampleIds)		# number of all samples
  nVars <- length(x) # number of variables

  ## Kriging system: C * wmu = D
  if(oneCondition) {
    C <- matrix(0, nrow=nSamples+1, ncol=nSamples+1)
    D <- matrix(0, nrow=nSamples+1, ncol=1)
  } else {
    C <- matrix(0, nrow=nSamples+nVars, ncol=nSamples+nVars)
    D <- matrix(0, nrow=nSamples+nVars, ncol=1)
  }

  # C matrix
  sampleIndex <- list()
  for(i in 1:nSamples) {
    if(length(sampleIndex) < i) sampleIndex[[i]] <- xAll$discretePoints$var_areaId == sampleIds[i]
    sampleI <- xAll$discretePoints[sampleIndex[[i]],]
    for(j in i:nSamples) {
      if(length(sampleIndex) < j) sampleIndex[[j]] <- xAll$discretePoints$var_areaId == sampleIds[j]
      sampleJ <- xAll$discretePoints[sampleIndex[[j]],]
      ptVgm <- ptVgms[[crossName(sampleI$varId[1], sampleJ$varId[1])]]
      C[i,j] <- ataCov(sampleI[,2:4], sampleJ[,2:4], ptVgm, longlat = longlat)
      C[j,i] <- C[i,j]
    }
  }

  if(oneCondition) {
    C[nSamples+1, ] <- 1
    C[, nSamples+1] <- 1
    C[nSamples+1, nSamples+1] <- 0
    D[nSamples+1] <- 1
  } else {
    for (i in 1:nVars) {
      indx <- xAll$areaValues$varId == varIds[i]
      C[nSamples+i, (1:nSamples)[indx]] <- 1
      C[(1:nSamples)[indx], nSamples+i] <- 1
    }
    D[nSamples + which(unknownVarId == varIds)] <- 1
  }

  unknownAreaIds <- sort(unique(unknown[,1]))

  krigOnce <- function(k) {
    curUnknown <- unknown[unknown[,1] == unknownAreaIds[k], 2:4]

    # D matrix
    for(i in 1:nSamples) {
      sampleI <- xAll$discretePoints[xAll$discretePoints$var_areaId == sampleIds[i],]
      ptVgm <- ptVgms[[crossName(sampleI$varId[1], unknownVarId)]]
      D[i] <- ataCov(sampleI[,2:4], curUnknown, ptVgm, longlat = longlat)
    }

    # solving
    wmu <- MASS::ginv(C) %*% D

    # estimation
    if(oneCondition) {
      if(is.null(meanVal)) {
        for (id in varIds) {
          meanVal <- rbind(meanVal, data.frame(varId=id, value=mean(x[[id]]$areaValues[,4])))
        }
      }
      rownames(meanVal) <- meanVal$varId

      w <- wmu[1:nSamples]
      w1 <- w[unknownVarId == xAll$areaValues$varId]
      yest <- sum(w1 * x[[unknownVarId]]$areaValues[,4])
      for (id in varIds[varIds != unknownVarId]) {
        w2 <- w[id == xAll$areaValues$varId]
        if (auxRatioAdj && abs(meanVal[id,2]) > 1e-6) {
          yest <- yest + sum(w2 * ((x[[id]]$areaValues[,4] - meanVal[id, 2])*(meanVal[unknownVarId, 2]/meanVal[id,2]) + meanVal[unknownVarId, 2]))
        } else {
          yest <- yest + sum(w2 * (x[[id]]$areaValues[,4] - meanVal[id, 2] + meanVal[unknownVarId, 2]))
        }
      }
    } else {
      w <- wmu[1:nSamples][unknownVarId == xAll$areaValues$varId]
      yest <- sum(w * x[[unknownVarId]]$areaValues[,4])
    }
    yvar <- ataCov(curUnknown, curUnknown, ptVgms[[unknownVarId]], longlat = longlat) - sum(wmu * D)

    if(!clarkAntiLog)
      return(data.frame(areaId=unknownAreaIds[k], pred=yest, var=yvar))
    else
      return(data.frame(areaId=unknownAreaIds[k], pred=yest, var=yvar,
                        pred.Clark=yest + sum(wmu * D)))
  }

  hasCluster <- ataIsClusterEnabled()
  if(showProgress) pb <- txtProgressBar(min=0, max=length(unknownAreaIds), width = 50, style = 3)

  if(!hasCluster || nopar || length(unknownAreaIds) == 1) {
    estResults <- c()
    for (k in 1:length(unknownAreaIds)) {
      estResults <- rbind(estResults, krigOnce(k))
      if(showProgress) setTxtProgressBar(pb, k)
    }
  } else {
    progress <- function(k) if(showProgress) setTxtProgressBar(pb, k)
    estResults <-
      foreach(k = 1:length(unknownAreaIds), .combine = rbind, .options.snow=list(progress=progress),
              .export = c("D","meanVal","crossName","ataCov","calcAreaCentroid", "spDistsNN"),
              .packages = c("sp","gstat")) %dopar% {
                krigOnce(k)
              }
    ataClusterClearObj()
  }

  if(showProgress) close(pb)

  unknownCenter <- calcAreaCentroid(unknown)
  estResults <- merge(unknownCenter, estResults)

  return(list(estimates = estResults, C = C, D = D))
}
