################################
## SCS - supporting functions ##
################################
##
## Author: Luca Del Core

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

#######
## scs
#######
##
## Description: Fit the SCS for a particular set of confounders and factors of interest.
##
scs <- function(x.L, y, nKTS.L, sp.ord.L, confounders.L, RELTOL = 1e-4, graphic = F, label = ""){
  
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
  

  res.opt <- nloptr(x0 = log(1:ncol(AT)),
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
                                "print_level" = 1))
  
  b_MLE <- res.opt$solution
  names(b_MLE) <- c("intercept",
                    tail(colnames(AT.cfd), -1)) 
  
  
  nDF <- ncol(AT) - sum(CM %*% AT %*% b_MLE == 0)
  s2_MLE <- as.numeric(t(log(y) - H%*%AT%*%b_MLE)%*%(log(y) - H%*%AT%*%b_MLE)/(nrow(H) - nDF))
  AIC_MLE <- 2*nDF - 2*sum(dnorm(x = log(y), mean = as.numeric(H%*%AT%*%b_MLE), sd = sqrt(s2_MLE), log = T)) + 2*nDF*(nDF+1)/(nrow(H) - nDF)
  BIC_MLE <- log(nrow(H))*nDF - 2*sum(dnorm(x = log(y), mean = as.numeric(H%*%AT%*%b_MLE), sd = sqrt(s2_MLE), log = T))
  
  b_MLE_cfd <- b_MLE[1:ncol(as.matrix(bdiag(1, bdiag(lapply(H_CM.L, function(v){v$AT})))))]
  
  zLim <- list(original = range(y),
               rescaled = c(1,6))
  
  if(graphic != FALSE){
    H.pred <- cbind(1, Reduce(cbind, lapply(H_CM.L, function(v){v$H.pred})))
    
    pdf(file = paste(currFigPath, "h_", graphic, "_confounders.pdf", sep = ""), width = 7, height = 5)
    par(mar = c(5,2,0,2))
    sc3 <- scatterplot3d(x.L[[confounders.L[1]]], 
                         x.L[[confounders.L[2]]], 
                         y, pch = 20,
                         zlim = zLim[[graphic]],
                         main = "", # paste("AICc = ", round(AIC_MLE,2), sep = ""), 
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
    fig_label(label, region="figure", pos="topleft", cex = 3, font = 2)
    dev.off()
  }
  
  ####
  res <- list()
  res$par <- b_MLE
  res$convergence <- res.opt$status
  res$message <- res.opt$message
  
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

##############
## strsplits
##############
##
## Description: Generalization of strsplit for multiple string splits.
##
strsplits <- function(x, splits, ...)
{
  for (split in splits)
  {
    x <- unlist(strsplit(x, split, ...))
  }
  return(x[!x == ""])
}

#########
## SRS_2
#########
##
## Description: reimplementation of SRS's SRS() function which allows the use of sparse matrices.
##
SRS_2 <- function (data, Cmin) 
{
  if (Cmin > min(colSums(data))) {
    print("ERROR: Cmin > minimum library size. Please select a Cmin that is <= the minimum library size of the dataset.")
  }
  else {
    if (Cmin < 0) {
      print("ERROR: Cmin < 0. Please select a Cmin >= 0.")
    }
    else {
      if (Cmin%%1 > 0) {
        print("ERROR: Please select a Cmin without decimal places.")
      }
      else {
        counter = 1
        for (i in seq(1, ncol(data), 1)) {
          if (i == 1) {
            fixed_factor <- (data[, i]/(sum(data[, i])/Cmin))
            assign(paste(colnames(data)[i], sep = ""), fixed_factor)
            fixed_factor_1 <- Matrix(get(colnames(data)[i]))
            colnames(fixed_factor_1)[i] <- colnames(data)[i]
          }
          else {
            fixed_factor <- (data[, i]/(sum(data[, i])/Cmin))
            assign(paste(colnames(data)[i], sep = ""), fixed_factor)
            fixed_factor_1 <- cbind(fixed_factor_1, fixed_factor)
            colnames(fixed_factor_1)[{
              counter = counter + 1
            }] <- colnames(data)[i]
          }
        }
        
        revtrunc_fixed_factor_1 <- floor(fixed_factor_1)
        
        diff_counts <- Cmin - colSums(revtrunc_fixed_factor_1)
        
        revtrunc <- function(x) {
          sign(x) * (x - floor(x))
        }
        revtrunc_fixed_factor <- (round(revtrunc(fixed_factor_1), 1e+07))
        
        x <- revtrunc_fixed_factor
        counter = 1
        for (i in seq(1, ncol(x), 1)) {
          if (i == 1) {
            if (diff_counts[i] == 0) {
              fixed_factor <- revtrunc_fixed_factor_1[, i]
              assign(paste(colnames(data)[i], sep = ""), fixed_factor)
              fixed_factor_1 <- Matrix(get(colnames(data)[i]))
              colnames(fixed_factor_1)[i] <- colnames(data)[i]
            }
            else {
              maxN <- function(x, N = diff_counts[i]) {
                len <- length(x)
                if (N > len) {
                  warning("N greater than length(x).  Setting N=length(x)")
                  N <- length(x)
                }
                sort(x, partial = len - N + 1)[len - N + 1]
              }
              max <- which(x[, i] == maxN(x[, i]), arr.ind = TRUE)
              
              normalization_value <- diff_counts[i] - sum(x[, i] > unique(x[, i][max]))
              
              lowest_level_choise <- matrix(which(x[, i] == unique(maxN(x[, i]))))
              
              if (sum(revtrunc_fixed_factor_1[, i][lowest_level_choise[, 1]]) == 0) {
                lowest_level <- as.numeric(as.vector(sample(as.factor(lowest_level_choise[, 1]), normalization_value, replace = F)))
                y <- t(Matrix(rep(0, length(x[, 1]))))
                y[lowest_level] = 1
                
              }
              else {
                sub_int <- matrix(subset(lowest_level_choise, (revtrunc_fixed_factor_1[, i][lowest_level_choise[, 1]] >= 1) == TRUE))
                sub_int_bind <- matrix(cbind(sub_int,  revtrunc_fixed_factor_1[, i][sub_int[,  1]]), ncol = 2)
                
                colnames(sub_int_bind) <- c("V1", "V2")
                
                sub_int_bind_ordered <- matrix(sub_int_bind[order(sub_int_bind[, "V2"], decreasing = TRUE), ], ncol = 2)
                
                colnames(sub_int_bind_ordered) <- c("V1", "V2")
                
                sub_int_bind_ordered_V1 <- sub_int_bind_ordered[, "V1"]
                sub_int_bind_ordered_V2 <- sub_int_bind_ordered[, "V2"]
                if ((length(unique(sub_int_bind_ordered_V2)) == 
                     1 & length(sub_int_bind_ordered_V2) > 
                     as.vector(normalization_value))) {
                  lowest_level <- as.numeric(as.vector(sample(as.factor(sub_int_bind_ordered_V1), normalization_value, replace = F)))
                  y <- t(Matrix(rep(0, length(x[, 1]))))
                  y[lowest_level] = 1
                  
                }
                else {
                  if (length(sub_int_bind_ordered_V1) > 
                      normalization_value) {
                    maxN_1 <- function(x, N = normalization_value) {
                      len <- length(x)
                      if (N > len) {
                        warning("N greater than length(x).  Setting N=length(x)")
                        N <- length(x)
                      }
                      sort(x, partial = len - N + 1)[len - N + 1]
                    }
                    max_1 <- which(Matrix(sub_int_bind_ordered_V2)[, 1] == maxN_1(Matrix(sub_int_bind_ordered_V2)[, 1]), arr.ind = TRUE)
                    
                    normalization_value_1 <- normalization_value - sum(Matrix(sub_int_bind_ordered_V2)[,  1] > unique(Matrix(sub_int_bind_ordered_V2)[, 1][max_1]))
                    
                    lowest_level_choise_1 <- Matrix(which(Matrix(sub_int_bind_ordered_V2)[, 1] == unique(maxN_1(Matrix(sub_int_bind_ordered_V2)[, 1]))))
                    
                    lowest_level <- as.numeric(as.vector(sample(as.factor(lowest_level_choise_1[, 1]), normalization_value_1, replace = F)))
                    lowest_level <- sub_int_bind_ordered_V1[lowest_level]
                    
                    lowest_level_1 <- sub_int_bind_ordered_V1[(Matrix(sub_int_bind_ordered_V2)[, 1] > unique(Matrix(sub_int_bind_ordered_V2)[, 1][max_1]))]
                    lowest_level <- c(lowest_level_1, lowest_level)
                    y <- t(Matrix(rep(0, length(x[, 1]))))
                    y[lowest_level] = 1
                  }
                  else {
                    if (length(sub_int_bind_ordered_V1) < 
                        normalization_value) {
                      sub_int_zeros <- Matrix(subset(lowest_level_choise, (revtrunc_fixed_factor_1[, i][lowest_level_choise[, 1]] < 1) == TRUE))
                      length(t(sub_int_zeros))
                      lowest_level_2 <- as.numeric(as.vector(sample(as.factor(sub_int_zeros[, 1]), (normalization_value - length(sub_int_bind_ordered_V1)), replace = F)))
                      
                      lowest_level_3 <- c(sub_int_bind_ordered_V1, lowest_level_2)
                      y <- t(Matrix(rep(0, length(x[, 1]))))
                      y[lowest_level_3] = 1
                      
                    }
                    else {
                      y <- t(Matrix(rep(0, length(x[, 1]))))
                      y[sub_int_bind_ordered_V1] = 1
                      
                    }
                  }
                }
              }
              SRS <- t(revtrunc_fixed_factor_1[, i] + ceiling(x[, i] > unique(x[, i][max])) + y)
              
              assign(paste(colnames(data)[i], sep = ""), SRS)
              fixed_factor_1 <- get(colnames(data)[i])
              colnames(fixed_factor_1)[i] <- colnames(data)[i]
            }
          }
          else {
            if (diff_counts[i] == 0) {
              fixed_factor <- revtrunc_fixed_factor_1[, i]
              assign(paste(colnames(data)[i], sep = ""), fixed_factor)
              fixed_factor_1 <- cbind(fixed_factor_1, fixed_factor)
              colnames(fixed_factor_1)[{
                counter = counter + 1
              }] <- colnames(data)[i]
            }
            else {
              maxN <- function(x, N = diff_counts[i]) {
                len <- length(x)
                if (N > len) {
                  warning("N greater than length(x).  Setting N=length(x)")
                  N <- length(x)
                }
                sort(x, partial = len - N + 1)[len - N + 1]
              }
              max <- which(x[, i] == maxN(x[, i]), arr.ind = TRUE)
              
              normalization_value <- diff_counts[i] - sum(x[, i] > unique(x[, i][max]))
              
              lowest_level_choise <- matrix(which(x[, i] == unique(maxN(x[, i]))))
              
              length(t(lowest_level_choise))
              if (sum(revtrunc_fixed_factor_1[, i][lowest_level_choise[, 1]]) == 0) {
                lowest_level <- as.numeric(as.vector(sample(as.factor(lowest_level_choise[, 1]), normalization_value, replace = F)))
                y <- t(Matrix(rep(0, length(x[, 1]))))
                y[lowest_level] = 1
                
              }
              else {
                sub_int <- matrix(subset(lowest_level_choise, (revtrunc_fixed_factor_1[, i][lowest_level_choise[, 1]] >= 1) == TRUE))
                sub_int_bind <- Matrix(cbind(sub_int, revtrunc_fixed_factor_1[, i][sub_int[, 1]]))
                colnames(sub_int_bind) <- c("V1", "V2")
                sub_int_bind_ordered <- matrix(sub_int_bind[order(sub_int_bind[, "V2"], decreasing = TRUE), ], ncol = 2)
                colnames(sub_int_bind_ordered) <- c("V1", "V2")
                
                sub_int_bind_ordered_V1 <- sub_int_bind_ordered[, "V1"]
                sub_int_bind_ordered_V2 <- sub_int_bind_ordered[, "V2"]
                if ((length(unique(sub_int_bind_ordered_V2)) == 1 & length(sub_int_bind_ordered_V2) > as.vector(normalization_value))) {
                  lowest_level <- as.numeric(as.vector(sample(as.factor(sub_int_bind_ordered_V1), normalization_value, replace = F)))
                  y <- t(Matrix(rep(0, length(x[, 1]))))
                  y[lowest_level] = 1
                  
                }
                else {
                  if (length(sub_int_bind_ordered_V1) > normalization_value) {
                    maxN_1 <- function(x, N = normalization_value) {
                      len <- length(x)
                      if (N > len) {
                        warning("N greater than length(x).  Setting N=length(x)")
                        N <- length(x)
                      }
                      sort(x, partial = len - N + 1)[len - N + 1]
                    }
                    max_1 <- which(Matrix(sub_int_bind_ordered_V2)[, 1] == maxN_1(Matrix(sub_int_bind_ordered_V2)[, 1]), arr.ind = TRUE)
                    
                    normalization_value_1 <- normalization_value - sum(Matrix(sub_int_bind_ordered_V2)[, 1] > unique(Matrix(sub_int_bind_ordered_V2)[, 1][max_1]))
                    
                    lowest_level_choise_1 <- Matrix(which(Matrix(sub_int_bind_ordered_V2)[, 1] == unique(maxN_1(Matrix(sub_int_bind_ordered_V2)[, 1]))))

                    lowest_level <- as.numeric(as.vector(sample(as.factor(lowest_level_choise_1[, 1]), normalization_value_1, replace = F)))
                    lowest_level <- sub_int_bind_ordered_V1[lowest_level]

                    lowest_level_1 <- sub_int_bind_ordered_V1[(Matrix(sub_int_bind_ordered_V2)[, 1] > unique(Matrix(sub_int_bind_ordered_V2)[, 1][max_1]))]
                    lowest_level <- c(lowest_level_1, lowest_level)
                    y <- t(Matrix(rep(0, length(x[, 1]))))
                    y[lowest_level] = 1
                  }
                  else {
                    if (length(sub_int_bind_ordered_V1) < normalization_value) {
                      sub_int_zeros <- Matrix(subset(lowest_level_choise, (revtrunc_fixed_factor_1[, i][lowest_level_choise[, 1]] < 1) == TRUE))
                      length(t(sub_int_zeros))
                      lowest_level_2 <- as.numeric(as.vector(sample(as.factor(sub_int_zeros[, 1]), (normalization_value - length(sub_int_bind_ordered_V1)), replace = F)))
                      lowest_level_3 <- c(sub_int_bind_ordered_V1, lowest_level_2)
                      y <- t(Matrix(rep(0, length(x[, 1]))))
                      y[lowest_level_3] = 1
                    }
                    else {
                      y <- t(Matrix(rep(0, length(x[, 1]))))
                      y[sub_int_bind_ordered_V1] = 1
                    }
                  }
                }
              }
              SRS <- t(revtrunc_fixed_factor_1[, i] + ceiling(x[, i] > unique(x[, i][max])) + y)
              assign(paste(colnames(data)[i], sep = ""), SRS)
              fixed_factor_1 <- cbind(fixed_factor_1, SRS)
              colnames(fixed_factor_1)[{
                counter = counter + 1
              }] <- colnames(data)[i]
            }
          }
        }
        SRS_output <- fixed_factor_1
        gc()
        return(SRS_output)
      }
    }
  }
}

############
## rrarefy2
############
##
## Description: reimplementation of VEGAN's rrarefy() function which allows the use of sparse matrices.
##
rrarefy2 <- function (x, sample) 
{
  if (!identical(all.equal(x, round(x)), TRUE)) 
    stop("function is meaningful only for integers (counts)")
  if (!is.integer(x)) 
    x <- round(x)
  if (ncol(x) == 1) 
    x <- t(x)
  if (length(sample) > 1 && length(sample) != nrow(x)) 
    stop(gettextf("length of 'sample' and number of rows of 'x' do not match"))
  sample <- rep(sample, length = nrow(x))
  for (i in 1:nrow(x)) {
    x[i, ] <- .Call(do_rrarefy, x[i, ], sample[i])
  }
  return(Matrix(x,sparse = T))
}
environment(rrarefy2) <- environment(rrarefy)
