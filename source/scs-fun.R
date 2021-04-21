################################
## SCS - supporting functions ##
################################
##
## Author: Luca Del Core

scs <- function(x.L, y, nKTS.L, sp.ord.L, confounders.L, RELTOL = 1e-4, graphic = F){
  
  H_CM.L <- lapply(1:length(x.L), function(v){
    
    x <- x.L[[v]]
    nKTS <- nKTS.L[[v]]
    sp.ord <- sp.ord.L[[v]]
    
    kts <- tail(head(seq(from = min(x), to = max(x), length.out = nKTS + 2), -1), -1)
    bdKts <- range(x)
    
    p <- nKTS + sp.ord
    H <- bs(x, 
            intercept = F, 
            degree = sp.ord, 
            knots = kts, 
            Boundary.knots = bdKts)
    colnames(H) <- rep(names(x.L)[v], ncol(H))
    
    H.sorted <- bs(x = sort(x, decreasing = F), 
                   intercept = F, 
                   degree = sp.ord, 
                   knots = kts, 
                   Boundary.knots = bdKts)
    colnames(H.sorted) <- rep(names(x.L)[v], ncol(H))
    
    if(v %% 2 == 1){
      H.pred <- bs(x = rep(seq(min(x), max(x), length.out = 100), each = 100), 
                   intercept = F, 
                   degree = sp.ord, 
                   knots = kts, 
                   Boundary.knots = bdKts)
    }else{
      H.pred <- bs(x = rep(seq(min(x), max(x), length.out = 100), times = 100), 
                   intercept = F, 
                   degree = sp.ord, 
                   knots = kts, 
                   Boundary.knots = bdKts)
    }
    
    dM <- tail(dbs(x = sort(x, decreasing = F), 
                   intercept = F, 
                   degree = sp.ord, 
                   knots = kts, 
                   Boundary.knots = bdKts, 
                   derivs = 1), 1)
    
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
  
  
  H <- as.matrix(cbind(H.cfd))
  CM <- cbind(CM.cfd)
  AT <- as.matrix(AT.cfd)
  colnames(AT) <- c("intercept",
                    tail(colnames(AT.cfd), -1)) 
  
  
  local_opts <- list( "algorithm" = "NLOPT_LD_LBFGS",
                      "xtol_rel"  = RELTOL)
  
  b_MLE <- nloptr(x0 = log(1:ncol(AT)),
                  eval_f = f,
                  eval_grad_f = g,
                  eval_g_ineq = inConstrF,
                  eval_jac_g_ineq = inConstrJac,
                  lb = c(-Inf, rep(1e-16, ncol(AT) - 1)), # c(-Inf, -Inf, rep(0, d-2))
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
                              "print_level" = 1))$solution
  
  names(b_MLE) <- c("intercept",
                    tail(colnames(AT.cfd), -1)) 
  
  
  nDF <- ncol(AT) - sum(CM %*% AT %*% b_MLE == 0)
  s2_MLE <- as.numeric(t(log(y) - H%*%AT%*%b_MLE)%*%(log(y) - H%*%AT%*%b_MLE)/(nrow(H) - nDF))
  AIC_MLE <- 2*nDF - 2*sum(dnorm(x = log(y), mean = as.numeric(H%*%AT%*%b_MLE), sd = sqrt(s2_MLE), log = T)) + 2*nDF*(nDF+1)/(nrow(H) - nDF)
  BIC_MLE <- log(nrow(H))*nDF - 2*sum(dnorm(x = log(y), mean = as.numeric(H%*%AT%*%b_MLE), sd = sqrt(s2_MLE), log = T))
  
  b_MLE_cfd <- b_MLE[1:ncol(as.matrix(bdiag(1, bdiag(lapply(H_CM.L, function(v){v$AT})))))]
  
  if(graphic){
    H.pred <- cbind(1, Reduce(cbind, lapply(H_CM.L, function(v){v$H.pred})))
    
    pdf(file = paste(currFigPath, "h_confounders.pdf", sep = ""), width = 7, height = 5)
    par(mar = c(5,2,2,2))
    sc3 <- scatterplot3d(x.L[[confounders.L[1]]], 
                         x.L[[confounders.L[2]]], 
                         y, pch = 20, 
                         main = paste("AICc = ", round(AIC_MLE,2), sep = ""), 
                         cex.main = 2, 
                         cex = 4, 
                         color = alpha("black", .5), # zlim = c(1, 6),
                         xlab = confounders.L[1], 
                         ylab = confounders.L[2], 
                         zlab = "entropy", 
                         cex.axis = 1.5, 
                         cex.lab = 1.5)
    
    sc3$points3d(rep(seq(min(x.L[[confounders.L[1]]]), max(x.L[[confounders.L[1]]]), length.out = 100), each = 100), 
                 rep(seq(min(x.L[[confounders.L[2]]]), max(x.L[[confounders.L[2]]]), length.out = 100), times = 100),
                 exp(H.pred %*% AT.cfd %*% b_MLE_cfd), col = alpha("green", .5), pch = 20, cex = .3)
    dev.off()
  }
  
  ####
  res <- list()
  res$par <- b_MLE
  
  w.conf <- which(colnames(H.cfd) %in% confounders.L)
  res$residuals <- exp(log(y) - H.cfd[,w.conf]%*%(AT.cfd%*%b_MLE_cfd)[w.conf])
  
  
  H.sorted <- cbind(1, Reduce(cbind, lapply(H_CM.L, function(v){v$H.sorted})))
  H.diff <- matrix(data = tail(H.sorted,1), nrow = nrow(H.cfd), ncol = ncol(H.cfd), byrow = T) - H.cfd
  res$residualsFromMax <- exp(log(y) + H.diff[,w.conf]%*%(AT.cfd%*%b_MLE_cfd)[w.conf])
  
  res$AIC <- AIC_MLE
  res$BIC <- BIC_MLE
  res$AT <- AT
  res$CM <- CM
  
  return(res)
}

f <- function(b, X, y, CM, AT){
  
  return(as.numeric(-2*t(b)%*%t(X%*%AT)%*%y + t(b)%*%t(X%*%AT)%*%X%*%AT%*%b))
}


g <- function(b, X, y, CM, AT){
  
  g <- as.numeric(-2*t(X%*%AT)%*%y + 2*t(X%*%AT)%*%X%*%AT%*%b)
  return(g)
}

inConstrF <- function(b, X, y, CM, AT){
  
  return(as.numeric(CM %*% AT %*% b))
}

inConstrJac <- function(b, X, y, CM, AT){
  
  gCM <- CM %*%AT
  return(gCM)
}

strsplits <- function(x, splits, ...)
{
  for (split in splits)
  {
    x <- unlist(strsplit(x, split, ...))
  }
  return(x[!x == ""])
}
