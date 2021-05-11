#############
## fig_label
#############
##
##  Description: Add a label to a plot into the external area.
##
fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}

###########
## get.AF
###########
##
## Description: Get all the matrices (design and constraints)
##              for the additional factors of interest
##
get.AF <- function(yX){
  addFactors <- lapply(as.vector(unique(yX$CellMarker)), function(mrk){
    
    res.vID <- lapply(as.vector(unique(yX$VectorID)), function(vID){
      sp.ord <- 2
      nKTS <- 2
      p <- nKTS + sp.ord
      
      yX_vecMrk <- yX[which(yX$VectorID == vID & yX$CellMarker == mrk),]
      x <- yX_vecMrk$TimePoint
      kts <- tail(head(seq(from = min(x), to = max(x), length.out = nKTS + 2), -1), -1)
      bdKts <- range(x)
      ## design:
      H <- bs(x, intercept = F, degree = sp.ord, knots = kts, Boundary.knots = bdKts)  
      rownames(H) <- rownames(yX_vecMrk)
      colnames(H) <- rep(paste(vID, mrk, collapse = "-"), ncol(H))
      H <- as.matrix(as.data.frame(H))
      ## design for prediction:
      H.pred <- bs(seq(min(x), max(x), length.out = 100), intercept = F, degree = sp.ord, knots = kts, Boundary.knots = bdKts)  
      colnames(H.pred) <- rep(paste(vID, mrk, collapse = "-"), ncol(H.pred))
      H.pred <- as.matrix(as.data.frame(H.pred))
      
      ## derivative of pline basis:
      dM <- tail(dbs(sort(x, decreasing = F), intercept = F, degree = sp.ord, knots = kts, Boundary.knots = bdKts, derivs = sp.ord - 1), 1)
      ## linear equality constraints for natural behaviour:
      AT <- rbind(diag(1, p - 1),
                  c(rep(0,p-2), -dM[,p-1]/dM[,p]))
      ## linear inequality constraints for monotonicity:
      CC.mnt <- diag(1, p)
      CC.mnt[abs(row(CC.mnt) - col(CC.mnt)) == 1 & row(CC.mnt) > col(CC.mnt)] <- -1
      CC.mnt <- CC.mnt[-1,]
      
      return(list(H = H,
                  H.pred = H.pred,
                  AT = AT,
                  CC.mnt = CC.mnt,
                  rn = rownames(H),
                  cn = colnames(H)
      ))
    })
    
    
    H <- cbind(1, bdiag(lapply(res.vID, function(v){v$H})))
    rownames(H) <- as.vector(unlist(lapply(res.vID, function(v){rownames(v$H)})))
    # colnames(H) <- c(paste("Intercept", mrk, sep = "-"), as.vector(unlist(lapply(res.vID, function(v){colnames(v$H)}))))
    colnames(H) <- c("intercept", as.vector(unlist(lapply(res.vID, function(v){colnames(v$H)}))))
    
    H.pred <- cbind(1, bdiag(lapply(res.vID, function(v){v$H.pred})))
    # colnames(H.pred) <- c(paste("Intercept", mrk, sep = "-"), as.vector(unlist(lapply(res.vID, function(v){colnames(v$H.pred)}))))
    colnames(H.pred) <- c("intercept", as.vector(unlist(lapply(res.vID, function(v){colnames(v$H.pred)}))))
    
    AT <- bdiag(1, bdiag(lapply(res.vID, function(v){v$AT})))
    colnames(AT) <- rep("af", ncol(AT))
    
    CC.mnt <- bdiag(1, bdiag(lapply(res.vID, function(v){v$CC.mnt})))
    
    
    return(list(H = H,
                H.pred = H.pred,
                AT = AT,
                CC.mnt = CC.mnt,
                rn = rownames(H),
                cn = colnames(H)
    ))
  })
  # return(addFactors)
  
  AT.af <- bdiag(lapply(addFactors, function(v){v$AT}))
  # colnames(AT.af) <- rep("af", ncol(AT.af))
  colnames(AT.af) <- rep(c("intercept", rep("af", ncol(AT.af)/length(addFactors) - 1)), times = length(addFactors))
  # AT.af <- diag(1, ncol(H.af))
  H.af <- bdiag(lapply(addFactors, function(v){v$H}))
  rownames(H.af) <- as.vector(unlist(lapply(addFactors, function(v){rownames(v$H)})))
  colnames(H.af) <- as.vector(unlist(lapply(addFactors, function(v){colnames(v$H)})))
  H.af <- H.af[rownames(yX),]
  C.af <- bdiag(lapply(addFactors, function(v){v$CC.mnt}))
  
  H.pred.af <- bdiag(lapply(addFactors, function(v){v$H.pred}))
  colnames(H.pred.af) <- as.vector(unlist(lapply(addFactors, function(v){colnames(v$H.pred)})))
  
  return(list(H = H.af,
              H.pred = H.pred.af,
              AT = AT.af,
              CC.mnt = C.af,
              rn = rownames(H.af),
              cn = colnames(H.af)
  ))
}

###########
## get.CFDs
###########
##
## Description: Get all the matrices (design and constraints)
##              for the confounders
##
get.CFDs <- function(x.L, nKTS.L, sp.ord.L){
  H_CM.L <- lapply(1:length(x.L), function(v){
    
    x <- x.L[[v]]
    nKTS <- nKTS.L[[v]]
    sp.ord <- sp.ord.L[[v]]
    
    kts <- tail(head(seq(from = min(x), to = max(x), length.out = nKTS + 2), -1), -1)
    bdKts <- range(x)
    
    p <- nKTS + sp.ord
    H <- bs(x, intercept = F, degree = sp.ord, knots = kts, Boundary.knots = bdKts)
    colnames(H) <- rep(names(x.L)[v], ncol(H))
    
    H.sorted <- bs(sort(x, decreasing = F), intercept = F, degree = sp.ord, knots = kts, Boundary.knots = bdKts)
    colnames(H.sorted) <- rep(names(x.L)[v], ncol(H))
    
    if(v %% 2 == 1){
      H.pred <- bs(rep(seq(min(x), max(x), length.out = 100), each = 100), intercept = F, degree = sp.ord, knots = kts, Boundary.knots = bdKts)
    }else{
      H.pred <- bs(rep(seq(min(x), max(x), length.out = 100), times = 100), intercept = F, degree = sp.ord, knots = kts, Boundary.knots = bdKts)
    }
    
    dM <- tail(dbs(sort(x, decreasing = F), intercept = F, degree = sp.ord, knots = kts, Boundary.knots = bdKts, derivs = sp.ord - 1), 1)
    
    CC.mnt <- diag(1, p)
    CC.mnt[abs(row(CC.mnt) - col(CC.mnt)) == 1 & row(CC.mnt) > col(CC.mnt)] <- -1
    CC.mnt <- CC.mnt[-1,]
    
    CC.cncv <- diag(2, p)
    CC.cncv[abs(row(CC.cncv) - col(CC.cncv)) == 1] <- -1
    CC.cncv <- matrix(CC.cncv[-c(1,p),], nrow = p - 2, byrow = F)
    
    CM <- as.matrix(-1*rbind(CC.mnt, CC.cncv))
    
    AT <- rbind(diag(1, p - 1),
                c(rep(0,p-2), -dM[,p-1]/dM[,p]))
    colnames(AT) <- rep(names(x.L)[v], ncol(AT))
    
    return(list(H = H,
                H.sorted = H.sorted,
                H.pred = H.pred,
                CM = CM,
                AT = AT))
  })
  
  AT.cfd <- as.matrix(bdiag(1, bdiag(lapply(H_CM.L, function(v){v$AT}))))
  colnames(AT.cfd) <- c("intercept", as.vector(unlist(lapply(H_CM.L, function(v){colnames(v$AT)}))))
  H.cfd <- cbind(1, Reduce(cbind, lapply(H_CM.L, function(v){v$H})))
  CM.cfd <- as.matrix(cbind(0, bdiag(lapply(H_CM.L, function(v){v$CM}))))
  
  return(list(H = H.cfd,
              H.sorted = cbind(1, Reduce(cbind, lapply(H_CM.L, function(v){v$H.sorted}))),
              H.pred = cbind(1, Reduce(cbind, lapply(H_CM.L, function(v){v$H.pred}))),
              CM = CM.cfd,
              AT = AT.cfd))
}

#######
## scs
#######
##
## Description: Fit the SCS for a particular set of confounders and factors of interest.
##
scs <- function(x.L, y, nKTS.L, sp.ord.L, confounders.L, AF.mats, RELTOL = 1e-4, graphic = F){
  
  H.af <- AF.mats$H
  AT.af <- AF.mats$AT
  C.af <- AF.mats$CC.mnt
  
  CFD.mats <- get.CFDs(x.L, nKTS.L, sp.ord.L)
  
  AT.cfd <- CFD.mats$AT
  H.cfd <- CFD.mats$H
  CM.cfd <- CFD.mats$CM
  
  
  H <- as.matrix(cbind(H.cfd, H.af))
  CM <- as.matrix(bdiag(CM.cfd, C.af))
  AT <- as.matrix(bdiag(AT.cfd, AT.af))
  colnames(AT) <- c("intercept",
                    tail(colnames(AT.cfd), -1),
                    colnames(AT.af)) 
  
  
  local_opts <- list( "algorithm" = "NLOPT_LD_LBFGS",
                      "xtol_rel"  = RELTOL)
  
  res.opt <- nloptr(x0 = c(log(1:ncol(AT.cfd)), (-ncol(AT.af)):(-1)),
                  eval_f = f,
                  eval_grad_f = g,
                  eval_g_ineq = inConstrF,
                  eval_jac_g_ineq = inConstrJac,
                  lb = rep(-Inf, ncol(AT)),
                  ub = rep(Inf, ncol(AT)),
                  X = H, 
                  y = log(y),
                  CM = CM,
                  AT = AT,
                  opts = list("algorithm" = "NLOPT_LD_AUGLAG", # NLOPT_LD_AUGLAG, NLOPT_LD_MMA, NLOPT_LD_LBFGS
                              "xtol_rel"= RELTOL, 
                              "maxeval" = 10000,
                              "local_opts" = local_opts,
                              "check_derivatives" = TRUE,
                              "print_level" = 1))
  b_MLE <- res.opt$solution
  
  names(b_MLE) <- c("intercept",
                    tail(colnames(AT.cfd), -1),
                    colnames(AT.af)) 
  
  
  nDF <- ncol(AT) - sum(CM %*% AT %*% b_MLE == 0)
  s2_MLE <- as.numeric(t(log(y) - H%*%AT%*%b_MLE)%*%(log(y) - H%*%AT%*%b_MLE)/(nrow(H) - nDF))
  AIC_MLE <- 2*nDF - 2*sum(dnorm(x = log(y), mean = as.numeric(H%*%AT%*%b_MLE), sd = sqrt(s2_MLE), log = T)) + 2*nDF*(nDF+1)/(nrow(H) - nDF)
  BIC_MLE <- log(nrow(H))*nDF - 2*sum(dnorm(x = log(y), mean = as.numeric(H%*%AT%*%b_MLE), sd = sqrt(s2_MLE), log = T))
  
  b_MLE_cfd <- b_MLE[1:ncol(CFD.mats$AT)]
  
  if(graphic){
    H.pred <- CFD.mats$H.pred
    
    par(mar = c(5,2,2,2))
    sc3 <- scatterplot3d(x.L[[1]], 
                         x.L[[2]], 
                         y, pch = 20, 
                         cex = 3, 
                         color = alpha("black", .5),
                         main = paste("AICc = ", round(AIC_MLE,2), sep = ""), cex.main = 2,
                         xlab = confounders.L[1], 
                         ylab = confounders.L[2], 
                         zlab = "entropy", 
                         cex.axis = 1.5, 
                         cex.lab = 1.5)
    
    b_MLE_cfd.avgInt <- c(mean(head(b_MLE[names(b_MLE) %in% "intercept"], 1) + tail(b_MLE[names(b_MLE) %in% "intercept"], -1)), 
                          tail(b_MLE_cfd,-1))
    
    sc3$points3d(rep(seq(min(x.L[[1]]), max(x.L[[1]]), length.out = 100), each = 100), 
                 rep(seq(min(x.L[[2]]), max(x.L[[2]]), length.out = 100), times = 100),
                 exp(H.pred %*% AT.cfd %*% b_MLE_cfd.avgInt), col = alpha("green", .2), pch = 20, cex = .3)
  }
  
  ####
  res <- list()
  res$par <- b_MLE
  res$convergence <- res.opt$status
  res$message <- res.opt$message
  
  res$residuals <- exp(log(y) - H.cfd[,which(colnames(H.cfd) %in% confounders.L)]%*%(AT.cfd%*%b_MLE_cfd)[which(colnames(H.cfd) %in% confounders.L)])
  
  
  H.sorted <- CFD.mats$H.sorted
  H.diff <- matrix(data = tail(H.sorted,1), nrow = nrow(H.cfd), ncol = ncol(H.cfd), byrow = T) - H.cfd
  
  res$residualsFromMax <- exp(log(y) + H.diff[,which(colnames(H.cfd) %in% confounders.L)]%*%(AT.cfd%*%b_MLE_cfd)[which(colnames(H.cfd) %in% confounders.L)])
  
  res$AIC <- AIC_MLE
  res$BIC <- BIC_MLE
  res$nDF <- nDF
  
  return(res)
}

#######
## fma
#######
##
## Description: Compute the Frequentist Model Averaging (FMA) estimator 
##              and all the corresponding statistics.
##
fma <- function(yX, x.L, nKTS.L, sp.ord.L, AF.mats, RELTOL, graphic = F){
  
  allAICRes <- apply(head(expand.grid(c(T,F), c(T,F), c(T,F), c(T,F)), -1), 1, 
                     FUN = function(model){
                       scs(x.L[t(model)], 
                           y = as.numeric(yX$h.idx), 
                           nKTS.L = nKTS.L[t(model)], 
                           sp.ord.L = sp.ord.L[t(model)], 
                           confounders.L = names(x.L[t(model)]), 
                           AF.mats = AF.mats, 
                           RELTOL = RELTOL, 
                           graphic = graphic)
                     })
  
  all_pars <- matrix(data = 0, nrow = length(allAICRes), ncol = length(allAICRes[[1]]$par))
  colnames(all_pars) <- names(allAICRes[[1]]$par)
  apply(matrix(1:length(allAICRes)), 1, function(model){all_pars[model, which(colnames(all_pars) %in% names(allAICRes[[model]]$par))] <<- allAICRes[[model]]$par})
  
  probs <- exp(-as.vector(unlist(lapply(allAICRes, function(l){l$BIC})))/2)/sum(exp(-as.vector(unlist(lapply(allAICRes, function(l){l$BIC})))/2))
  varProbs <- apply(head(expand.grid(c(T,F), c(T,F), c(T,F), c(T,F)), -1), 2, function(v){sum(probs[v])})
  names(varProbs) <- c("DNA", "VCN", "PS", "SD")
  b_FMA <- colSums(all_pars * t(t(probs))[,rep(1, ncol(all_pars))])
  
  H.af <- AF.mats$H
  AT.af <- AF.mats$AT
  C.af <- AF.mats$CC.mnt
  
  CFD.mats <- get.CFDs(x.L, nKTS.L, sp.ord.L)
  
  AT.cfd <- CFD.mats$AT
  H.cfd <- CFD.mats$H
  CM.cfd <- CFD.mats$CM
  
  H <- as.matrix(cbind(H.cfd, H.af))
  CM <- as.matrix(bdiag(CM.cfd, C.af))
  AT <- as.matrix(bdiag(AT.cfd, AT.af))
  colnames(AT) <- c("intercept",
                    tail(colnames(AT.cfd), -1),
                    colnames(AT.af)) 
  
  nDF <- ncol(AT) - sum(CM %*% AT %*% b_FMA == 0)
  s2_MLE <- as.numeric(t(log(yX$h.idx) - H%*%AT%*%b_FMA)%*%(log(yX$h.idx) - H%*%AT%*%b_FMA)/(nrow(H) - nDF))
  
  res <- list()
  res$par <- b_FMA
  res$s2_MLE <- s2_MLE
  
  H.cfd <- CFD.mats$H
  AT.cfd <- CFD.mats$AT
  b_MLE_cfd <- b_FMA[1:ncol(AT.cfd)]
  res$par.cfd <- b_MLE_cfd
  
  res$residuals <- exp(log(yX$h.idx) - H.cfd[,which(colnames(H.cfd) %in% confounders.L)]%*%(AT.cfd%*%b_MLE_cfd)[which(colnames(H.cfd) %in% confounders.L)])
  
  
  H.sorted <- get.CFDs(x.L, nKTS.L, sp.ord.L)$H.sorted
  H.diff <- matrix(data = tail(H.sorted,1), nrow = nrow(H.cfd), ncol = ncol(H.cfd), byrow = T) - H.cfd
  
  res$residualsFromMax <- exp(log(yX$h.idx) + H.diff[,which(colnames(H.cfd) %in% confounders.L)]%*%(AT.cfd%*%b_MLE_cfd)[which(colnames(H.cfd) %in% confounders.L)])
  
  res$probs <- probs
  res$varProbs <- varProbs
  
  return(res)
}

#######
## f
#######
##
## Description: negative log-likelihood function
##
f <- function(b, X, y, CM, AT){
  
  return(as.numeric(-2*t(b)%*%t(X%*%AT)%*%y + t(b)%*%t(X%*%AT)%*%X%*%AT%*%b))
}

#######
## g
#######
##
## Description: gradient of the negative log-likelihood function
##
g <- function(b, X, y, CM, AT){
  
  g <- as.numeric(-2*t(X%*%AT)%*%y + 2*t(X%*%AT)%*%X%*%AT%*%b)
  return(g)
}

##############
## inConstrF
##############
##
## Description: inequality constraints
##
inConstrF <- function(b, X, y, CM, AT){
  
  return(as.numeric(CM %*% AT %*% b))
}

##############
## inConstrJac
##############
##
## Description: Jacobian of the inequality constraints
##
inConstrJac <- function(b, X, y, CM, AT){
  
  gCM <- CM %*%AT
  return(gCM)
}
