#####################
## SCS - main code ##
#####################
##
## Author: Luca Del Core

cat("\nInstall/load packages")
inst.pkgs <- installed.packages()

## required packages:
l.pkgs <- c("Matrix", 
            "splines", 
            "splines2", 
            "nloptr", 
            "scatterplot3d", 
            "scales")
## check if packages are installed
lapply(l.pkgs, function(pkg){
  if(!(pkg %in% rownames(inst.pkgs))){
    install.packages(pkg)
  }
})

lapply(l.pkgs, function(pkg){library(pkg,character.only=TRUE)})


##############
## currTime ##
##############

currTime <- Sys.time()
currTime <- gsub(" ", "_", currTime)
currTime <- gsub("-", "", currTime)
currTime <- gsub(":", "", currTime)
currTime <- substr(currTime, start = 1, stop = 13)

#############
## folders ##
#############

sourcePath <- paste("./source/")
isMatrices <- "./input/MS/matrices/"

StatsFolder <- "./input/MS/stats/"
allResPath <- "./results/"
resPath <- "./results/MS/"

currAllResPath <- paste(allResPath, "/", currTime, sep = "")
currResPath <- paste(resPath, "/", currTime, sep = "")
currOutPath <- paste(currResPath, "/output/", sep = "")
currFigPath <- paste(currResPath, "/figures/", sep = "")

ifelse(!dir.exists(file.path(allResPath)), dir.create(file.path(allResPath)), FALSE)
ifelse(!dir.exists(file.path(resPath)), dir.create(file.path(resPath)), FALSE)
ifelse(!dir.exists(file.path(currResPath)), dir.create(file.path(currResPath)), FALSE)
ifelse(!dir.exists(file.path(currOutPath)), dir.create(file.path(currOutPath)), FALSE)
ifelse(!dir.exists(file.path(currFigPath)), dir.create(file.path(currFigPath)), FALSE)


cat("\nload functions")
## load functions:
source(paste(sourcePath,"functions.R" , sep = ""))

cat("\nload data")
metadata <- read.csv(paste(StatsFolder, "/metadata_MS.csv", sep = ""), sep = ",", dec = ".", header = T, row.names = 1, stringsAsFactors = FALSE)
ISmat <- read.csv(paste(isMatrices, "/ISmat_MS.tsv", sep = ""), sep = "\t", header = T, row.names = NULL, check.names = FALSE) 
ISmat <- Matrix(as.matrix(ISmat), sparse = TRUE)

## compute Shannon entropy:
metadata$h.idx <- NA
h.idx <- apply(ISmat, 2, function(s){
  s <- s[s>0];
  p <- s/sum(s);
  h <- -sum(p*log(p));
  return(h)
})
## compute sequencing depth:
metadata$seqDepth <- NA
seq.depth <- apply(ISmat, 2, function(s){
  s <- s[s>0];
  sd <- sum(s)
  return(sd)
})

metadata[names(h.idx), "h.idx"] <- h.idx
metadata[names(seq.depth), "seqDepth"] <- seq.depth 
yX <- metadata


cat("\nSCS rescaling...")
## confounding factors:
x.L <- list(DNA = as.numeric(yX[,"DNAngUsed"]),
            VCN = as.numeric(yX[,"VCN"]),
            PS = as.numeric(yX[,"N.mice.pool"]),
            SD = as.numeric(yX[,"seqDepth"]))
nKTS.L <- list(kts.1 = 2, kts.2 = 2, kts.3 = 2, kts.4 = 2)
sp.ord.L <- list(sp.ord.1 = 2, sp.ord.2 = 2, sp.ord.3 = 2, sp.ord.4 = 2)
confounders.L <- names(x.L)

## Generating matrix of additional factors of interest:
vecMrkComb <- as.matrix(unique.matrix(expand.grid(yX$VectorID, as.character(yX$CellMarker)), MARGIN = 1))
rownames(vecMrkComb) <- 1:nrow(vecMrkComb)
AF.mats <- get.AF(yX)


## Frequentist Model Averaging (FMA):
res.FMA <- fma(yX, x.L, nKTS.L, sp.ord.L, AF.mats, RELTOL = 1e-4, graphic = F)
yX$h.idx.res <- NA
yX$h.idx.res <- as.numeric(res.FMA$residuals)
b_MLE <- res.FMA$par
b_MLE_cfd <- res.FMA$par.cfd


###########
## Plots ##
###########

pdf(file = paste(currFigPath, "h_confounders.pdf", sep = ""), width = 12, height = 8)
par(mar = c(5,5,2,2), mfrow = c(2,2))
plot(yX$DNAngUsed, yX$h.idx, cex.axis = 2, cex.lab = 2, xlab = "DNA", ylab = "h.idx", pch = 20, cex = 3, col = alpha("black", .5))
fig_label("a", region="figure", pos="topleft", cex = 3, font = 2)
plot(yX$VCN, yX$h.idx, cex.axis = 2, cex.lab = 2, xlab = "VCN", ylab = "h.idx", pch = 20, cex = 3, col = alpha("black", .5))
fig_label("b", region="figure", pos="topleft", cex = 3, font = 2)
plot(yX$N.mice.pool, yX$h.idx, cex.axis = 2, cex.lab = 2, xlab = "pool size", ylab = "h.idx", pch = 20, cex = 3, col = alpha("black", .5))
fig_label("c", region="figure", pos="topleft", cex = 3, font = 2)
plot(yX$seqDepth, yX$h.idx, cex.axis = 2, cex.lab = 2, xlab = "seq. depth", ylab = "h.idx", pch = 20, cex = 3, col = alpha("black", .5))
fig_label("d", region="figure", pos="topleft", cex = 3, font = 2)
dev.off()

pdf(file = paste(currFigPath, "frequentistPosterior.pdf", sep = ""), width = 7, height = 5)
par(mar = c(5,5,2,2), mfrow = c(1,1))
barplot(res.FMA$probs, ylab = "posterior probability", xlab = "model", 
        cex.axis = 2, cex.lab = 2, col = "black", font = 2, cex.names = 2,
        names.arg = 1:15)
fig_label("a", region="figure", pos="topleft", cex = 3, font = 2)
dev.off()


pdf(file = paste(currFigPath, "frequentistPosterior_VarProbs.pdf", sep = ""), width = 7, height = 5)
par(mar = c(5,5,2,2), mfrow = c(1,1))
barplot(res.FMA$varProbs, ylab = "posterior probability", xlab = "variable", 
        cex.axis = 2, cex.lab = 2, col = "black", font = 2, cex.names = 2,
        names.arg = names(res.FMA$varProbs))
fig_label("b", region="figure", pos="topleft", cex = 3, font = 2)
dev.off()

#################
## Fitting SCS ##
#################
##

mrk_lab <- c("a", "b", "c", "d")
names(mrk_lab) <- as.vector(unique(yX$CellMarker))
pdf(file = paste(currFigPath, "rescaledFitted_FMA_smoothing.pdf", sep = ""), width = 12, height = 8)
par(mar = c(5,5,2,2), mfrow = c(2,2))
for (mrk in as.vector(unique(yX$CellMarker))) {
  
  xLim <- range(yX[, "TimePoint"])
  yLim <- range(yX[, "h.idx.res"])
  
  vecMrk <- which(vecMrkComb[,2] == mrk)[1]
  
  yX_vecMrk <- yX[which(yX$VectorID == vecMrkComb[vecMrk,1] & yX$CellMarker == vecMrkComb[vecMrk,2]),]
  x <- yX_vecMrk$TimePoint
  y.pred <- exp((cbind(1, AF.mats$H.pred) %*% bdiag(1, AF.mats$AT) %*% b_MLE[names(b_MLE) %in% c("intercept","af")])[((vecMrk - 1)*100 + 1):((vecMrk - 1)*100 + 100)])
  
  q.t <- qt(p = 1 - .05/2, df = ncol(AF.mats$AT) - sum(AF.mats$CC.mnt %*% AF.mats$AT %*% tail(b_MLE, -length(b_MLE_cfd)) == 0))
  
  y.pred.sd <- apply(matrix(1:nrow(AF.mats$H.pred)), 1, function(r){
    sqrt(res.FMA$s2_MLE * as.numeric((AF.mats$H.pred[r,]%*%AF.mats$AT) %*% solve(t(AF.mats$H.pred%*%AF.mats$AT)%*%(AF.mats$H.pred%*%AF.mats$AT)) %*% t(AF.mats$H.pred[r,]%*%AF.mats$AT))) * q.t
  })[((vecMrk - 1)*100 + 1):((vecMrk - 1)*100 + 100)]
  
  y.pred.sd.minus <- exp((cbind(1, AF.mats$H.pred) %*% bdiag(1, AF.mats$AT) %*% b_MLE[names(b_MLE) %in% c("intercept","af")])[((vecMrk - 1)*100 + 1):((vecMrk - 1)*100 + 100)] - y.pred.sd)
  y.pred.sd.plus <- exp((cbind(1, AF.mats$H.pred) %*% bdiag(1, AF.mats$AT) %*% b_MLE[names(b_MLE) %in% c("intercept","af")])[((vecMrk - 1)*100 + 1):((vecMrk - 1)*100 + 100)] + y.pred.sd)
  
  plot(yX_vecMrk$TimePoint, yX_vecMrk$h.idx.res, xlim = xLim, ylim = yLim,
       pch = 19, xlab = "t", ylab = "h idx.", cex.axis = 2, cex.lab = 2, # pch = pchPCR[yX_vecMrk$PCRMethod]
       main = mrk, cex.main = 2, col = alpha("red", .5), cex = 2)
  lines(x = seq(min(x), max(x), length.out = 100),
        y = y.pred, lwd = 3, col = "red")
  lines(x = seq(min(x), max(x), length.out = 100),
        y = y.pred.sd.minus, lwd = 3, col = "red", lty = 2)
  lines(x = seq(min(x), max(x), length.out = 100),
        y = y.pred.sd.plus, lwd = 3, col = "red", lty = 2)
  
  vecMrk <- which(vecMrkComb[,2] == mrk)[2]
  
  yX_vecMrk <- yX[which(yX$VectorID == vecMrkComb[vecMrk,1] & yX$CellMarker == vecMrkComb[vecMrk,2]),]
  x <- yX_vecMrk$TimePoint
  y.pred <- exp((cbind(1, AF.mats$H.pred) %*% bdiag(1, AF.mats$AT) %*% b_MLE[names(b_MLE) %in% c("intercept", "af")])[((vecMrk - 1)*100 + 1):((vecMrk - 1)*100 + 100)])
  
  y.pred.sd <- apply(matrix(1:nrow(AF.mats$H.pred)), 1, function(r){
    sqrt(res.FMA$s2_MLE * as.numeric((AF.mats$H.pred[r,]%*%AF.mats$AT) %*% solve(t(AF.mats$H.pred%*%AF.mats$AT)%*%(AF.mats$H.pred%*%AF.mats$AT)) %*% t(AF.mats$H.pred[r,]%*%AF.mats$AT))) * q.t
  })[((vecMrk - 1)*100 + 1):((vecMrk - 1)*100 + 100)]
  
  y.pred.sd.minus <- exp((cbind(1, AF.mats$H.pred) %*% bdiag(1, AF.mats$AT) %*% b_MLE[names(b_MLE) %in% c("intercept","af")])[((vecMrk - 1)*100 + 1):((vecMrk - 1)*100 + 100)] - y.pred.sd)
  y.pred.sd.plus <- exp((cbind(1, AF.mats$H.pred) %*% bdiag(1, AF.mats$AT) %*% b_MLE[names(b_MLE) %in% c("intercept","af")])[((vecMrk - 1)*100 + 1):((vecMrk - 1)*100 + 100)] + y.pred.sd)
  
  points(yX_vecMrk$TimePoint, yX_vecMrk$h.idx.res, 
         pch = 17, col = alpha("blue", .5), cex = 2)
  lines(x = seq(min(x), max(x), length.out = 100),
        y = y.pred, lwd = 3, col = "blue")
  lines(x = seq(min(x), max(x), length.out = 100),
        y = y.pred.sd.minus, lwd = 3, col = "blue", lty = 2)
  lines(x = seq(min(x), max(x), length.out = 100),
        y = y.pred.sd.plus, lwd = 3, col = "blue", lty = 2)
  
  fig_label(mrk_lab[mrk], region="figure", pos="topleft", cex = 3, font = 2)
}
dev.off()



## Fit standard LM:
fit <- lm(formula = "log(h.idx) ~ bs(x = TimePoint, knots = quantile(TimePoint, probs = c(.33,.66)), degree = 2)*VectorID*CellMarker", data = yX)

mrks <- as.vector(unique(yX$CellMarker))
vcts <- as.vector(unique(yX$VectorID))
newdat_all <- do.call(rbind,lapply(mrks, function(mrk){
  do.call(rbind, lapply(vcts, function(vct){
    return(data.frame(TimePoint = seq(min(yX[which(yX$CellMarker == mrk & yX$VectorID == vct),]$TimePoint), 
                                      max(yX[which(yX$CellMarker == mrk & yX$VectorID == vct),]$TimePoint), length.out = 100),
                      CellMarker = mrk,
                      VectorID = vct))
  }))
}))

conf_interval_all <- exp(predict(fit, 
                                 newdata=newdat_all, 
                                 interval="confidence",
                                 level = 0.95))

## Plot results of standard LM:
pdf(file = paste(currFigPath, "standardLM.pdf", sep = ""), width = 12, height = 8)
par(mar = c(5,5,2,2), mfrow = c(2,2))
for (mrk in unique(yX$CellMarker)) {
  yX_mrk <- yX[which(yX$CellMarker == mrk),]
  
  
  yX_mrk_vID <- yX[which(yX$CellMarker == mrk & yX$VectorID == "LV.SF.LTR"),]
  newdat = data.frame(TimePoint = seq(min(yX_mrk_vID$TimePoint), max(yX_mrk_vID$TimePoint), length.out = 100),
                      CellMarker = mrk,
                      VectorID = "LV.SF.LTR")
  newdat$pred = exp(predict(fit, newdata = newdat))
  conf_interval <- exp(predict(fit, 
                               newdata=newdat, 
                               interval="confidence",
                               level = 0.95))
  
  plot(h.idx ~ TimePoint, data = yX_mrk_vID, cex = 2,
       main = mrk, cex.main = 2,
       pch = 19, cex.axis = 2, cex.lab = 2, col = alpha("red", .5),
       xlim = range(yX$TimePoint),
       ylim = range(yX_mrk$h.idx, conf_interval_all[,2:3]), xlab = "t", ylab = "h idx.")
  with(newdat, lines(x = TimePoint, y = pred, col = "red", lwd = 3))
  
  lines(newdat$TimePoint, conf_interval[,2], col="red", lty=2, lwd = 3)
  lines(newdat$TimePoint, conf_interval[,3], col="red", lty=2, lwd = 3)
  
  ####
  yX_mrk_vID <- yX[which(yX$CellMarker == mrk & yX$VectorID == "PGK"),]
  newdat = data.frame(TimePoint = seq(min(yX_mrk_vID$TimePoint), max(yX_mrk_vID$TimePoint), length.out = 100),
                      CellMarker = mrk,
                      VectorID = "PGK")
  newdat$pred = exp(predict(fit, newdata = newdat))
  conf_interval <- exp(predict(fit, 
                               newdata=newdat, 
                               interval="confidence",
                               level = 0.95))
  points(h.idx ~ TimePoint, data = yX_mrk_vID, cex = 2,
         pch = 17, col = alpha("blue", .5))
  with(newdat, lines(x = TimePoint, y = pred, col = "blue", lwd = 3))
  
  lines(newdat$TimePoint, conf_interval[,2], col="blue", lty=2, lwd = 3)
  lines(newdat$TimePoint, conf_interval[,3], col="blue", lty=2, lwd = 3)
  fig_label(mrk_lab[mrk], region="figure", pos="topleft", cex = 3, font = 2)
  
}
dev.off()

## Plot legends:
pdf(file = paste(currFigPath, "vectorLegends.pdf", sep = ""), width = 16, height = 4)
par(mar = c(5,5,2,2), mfrow = c(2,2))
plot.new()
legend(x = "center", legend = c("PGK", "LTR"), col = c("blue", "red"), pch = c(17, 20), cex = 4, lwd = 7, lty = 1, horiz = T)
dev.off()

