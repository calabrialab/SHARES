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
            "scales",
            "parallel",
            "vegan",
            "beeswarm")
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
isMatrices <- "./input/VA/matrices/"

StatsFolder <- "./input/VA/stats/"
allResPath <- "./results/"
resPath <- "./results/VA/"

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
metadata <- read.csv(paste(StatsFolder, "/metadata_VA.csv", sep = ""), sep = ",", dec = ".", header = T, row.names = 1, stringsAsFactors = FALSE)
ISmat <- read.csv(paste(isMatrices, "/ISmat_VA.tsv", sep = ""), sep = "\t", header = T, row.names = NULL, check.names = FALSE) 
ISmat <- Matrix(as.matrix(ISmat), sparse = TRUE)

## compute Shannon entropy:
metadata$h.idx<- NA
h.idx <- apply(ISmat, 2, function(s){
  s <- s[s>0];
  p <- s/sum(s);
  h <- -sum(p*log(p));
  return(h)
})
## compute sequencing depth:
metadata$seqDepth<- NA
seq.depth <- apply(ISmat, 2, function(s){
  s <- s[s>0];
  sd <- sum(s)
  return(sd)
})

metadata[names(h.idx), "h.idx"] <- h.idx
metadata[names(seq.depth), "seqDepth"] <- seq.depth 

yX_JY <- metadata

cat("\nSCS rescaling...")
## confounding factors:
x.L <- list(DNA = as.numeric(yX_JY[,"DNAngUsed"]),
            VCN = as.numeric(yX_JY[,"VCN"]),
            seqDepth = as.numeric(yX_JY[,"seqDepth"]))
nKTS.L <- list(kts.1 = 2, kts.2 = 2, kts.3 = 2)
sp.ord.L <- list(sp.ord.1 = 2, sp.ord.2 = 2, sp.ord.3 = 2)
confounders.L <- names(x.L)

## Apply SCS on VA data:
res <- scs(x.L = x.L, 
           y = as.numeric(yX_JY$h.idx), 
           nKTS.L = nKTS.L, 
           sp.ord.L = sp.ord.L, 
           confounders.L = confounders.L, 
           RELTOL = 1e-8, 
           graphic = "original",
           label = "a")
yX_JY$h.idx.res <- NA
yX_JY$h.idx.res <- as.numeric(res$residuals)


#################################################################################################
############################# COMPARISON WITH OTHER METHODS #####################################
#################################################################################################

cMin <- min(colSums(ISmat)) ## select the rarefaction level
set.seed(1) ## set seed of replicability
nRep <- 1000 ## set the number of random sub-samplings
cl <- makeCluster(getOption("cl.cores", round(detectCores()/2))) ## select the number of cores
clusterExport(cl = cl, varlist = list("ISmat", "cMin", "rrarefy2", "SRS_2")) # export variables to the cluster
# export packages to the cluster:
clusterEvalQ(cl, library("Matrix")) 
clusterEvalQ(cl, library("vegan"))

######## RAREFACTION ########
h.idx.RAR_allRep <- parApply(cl = cl, X = matrix(1:nRep), MARGIN = 1, FUN = function(r){
  message(r)
  ISmat_JY_RAR <- t(Matrix(rrarefy2(x = t(ISmat), cMin)));
  h.idx.RAR <- apply(ISmat_JY_RAR, 2, function(s){
    s <- s[s>0];
    p <- s/sum(s);
    h <- -sum(p*log(p));
    return(h)
  });
  return(h.idx.RAR)
})
h.idx.RAR_avg <- apply(h.idx.RAR_allRep, 1, mean)
yX_JY$h.idx.res.RAR <- NA
yX_JY[names(h.idx.RAR_avg), "h.idx.res.RAR"] <- h.idx.RAR_avg
h.idx.RAR_sd <- apply(h.idx.RAR_allRep, 1, sd)
yX_JY$h.idx.res.RAR.sd <- NA
yX_JY[names(h.idx.RAR_sd), "h.idx.res.RAR.sd"] <- h.idx.RAR_sd

######### SRS #########
h.idx.SRS_allRep <- parApply(cl = cl, X = matrix(1:nRep), MARGIN = 1, FUN = function(r){
  message(r)
  ISmat_JY_SRS <- SRS_2(ISmat, cMin)
  
  h.idx.SRS <- apply(ISmat_JY_SRS, 2, function(s){
    s <- s[s>0];
    p <- s/sum(s);
    h <- -sum(p*log(p));
    return(h)
  });
  
  return(h.idx.SRS)
})

h.idx.SRS_avg <- apply(h.idx.SRS_allRep, 1, mean)
yX_JY$h.idx.res.SRS <- NA
yX_JY[names(h.idx.SRS_avg), "h.idx.res.SRS"] <- h.idx.SRS_avg
h.idx.SRS_sd <- apply(h.idx.SRS_allRep, 1, sd)
yX_JY$h.idx.res.SRS.sd <- NA
yX_JY[names(h.idx.SRS_sd), "h.idx.res.SRS.sd"] <- h.idx.SRS_sd


###########
## Plots ##
###########

cors <- cor(yX_JY[,c("DNAngUsed", "VCN", "seqDepth", "h.idx", "h.idx.res","h.idx.res.RAR","h.idx.res.SRS")], method = "spearman")[4:7,1:3]
corPvals <- t(apply(matrix(rownames(cors)), 1, function(r){apply(matrix(colnames(cors)), 1, function(c){as.numeric(cor.test(x = yX_JY[,c], y = yX_JY[,r], method = "spearman")$p.value)})}))
dimnames(corPvals) <- dimnames(cors)

rownames(cors) <- c("OBS","SCS","RAR","SRS")
rownames(corPvals) <- rownames(cors)
cors <- cors[c("OBS", "RAR", "SRS", "SCS"),]
corPvals <- corPvals[c("OBS", "RAR", "SRS", "SCS"),]

pdf(file = paste(currFigPath, "Comparisons_correlations.pdf", sep = ""), width = 12, height = 4)
par(mar = c(5,6,2,2), mfrow = c(1,3))
bp <- barplot(abs(cors[,1]), cex.axis = 2, cex.names = 2, names.arg = rownames(cors[,]), col = "black", ylab = expression(rho(h,DNA)), cex.lab = 2)
lab <- attr(regexpr("(?<=\\.)0+", corPvals[,1], perl = TRUE), "match.length")
lab[lab < 0] <- 0
lab <- apply(matrix(lab), 1, FUN = function(l){if(l > 0){
  return(paste(rep("*", l), collapse = ""))
}else{
  return("")
}})
text(bp, abs(cors[,1])*.9, 
     labels = lab, col = "white", font = 2, srt = 90, cex = 4 ,pos = 1)
fig_label("a", region="figure", pos="topleft", cex = 3, font = 2)
bp <- barplot(abs(cors[,2]), cex.axis = 2, cex.names = 2, names.arg = rownames(cors[,]), col = "black", ylab = expression(rho(h,VCN)), cex.lab = 2)
lab <- attr(regexpr("(?<=\\.)0+", formatC(corPvals[,2], format = "f", digits = 10), perl = TRUE), "match.length")
lab[lab < 0] <- 0
lab <- apply(matrix(lab), 1, FUN = function(l){if(l > 0){
  return(paste(rep("*", l), collapse = ""))
}else{
  return("")
}})
text(bp, abs(cors[,2])*.8, 
     labels = lab, col = "white", font = 2, srt = 90, cex = 4, pos = 1)
fig_label("b", region="figure", pos="topleft", cex = 3, font = 2)
bp <- barplot(abs(cors[,3]), cex.axis = 2, cex.names = 2, names.arg = rownames(cors[,]), col = "black", ylab = expression(rho(h,SD)), cex.lab = 2)
lab <- attr(regexpr("(?<=\\.)0+", formatC(corPvals[,3], format = "f", digits = 10), perl = TRUE), "match.length")
lab[lab < 0] <- 0
lab <- apply(matrix(lab), 1, FUN = function(l){if(l > 0){
  return(paste(rep("*", l), collapse = ""))
}else{
  return("")
}})
text(bp, abs(cors[,3])*.7, 
     labels = lab, col = "white", font = 2, srt = 90, cex = 4, pos = 1)
fig_label("c", region="figure", pos="topleft", cex = 3, font = 2)
dev.off()


pdf(file = paste(currFigPath, "Comparisons_boxplots.pdf", sep = ""), width = 7, height = 5)
boxplot(yX_JY[, c("h.idx", "h.idx.res")], lwd = 3, cex.axis = 2, cex.lab = 2, cex = 2, names = c("OBS", "SCS"), outline = F)
beeswarm(yX_JY[, c("h.idx", "h.idx.res")], add=T, pch = 20, cex = 1.5, col = "red")
dev.off()

pdf(file = paste(currFigPath, "Comparisons_boxplots_allMethods.pdf", sep = ""), width = 15, height = 5)
boxplot(yX_JY[, c("h.idx", "h.idx.res", "h.idx.res.RAR", "h.idx.res.SRS")], lwd = 3, cex.axis = 2, cex.lab = 2, cex = 2, names = c("OBS", "SCS" ,"RAR", "SRS"), outline = F)
beeswarm(yX_JY[, c("h.idx", "h.idx.res", "h.idx.res.RAR", "h.idx.res.SRS")], add=T, pch = 20, cex = 1.5, col = "red")
dev.off()


scs(x.L = x.L, 
    y = as.numeric(yX_JY$h.idx.res), 
    nKTS.L = nKTS.L, 
    sp.ord.L = sp.ord.L, 
    confounders.L = confounders.L, 
    RELTOL = 1e-8, 
    graphic = "rescaled",
    label = "b")

pdf(file = paste(currFigPath, "h_confounders_VA.pdf", sep = ""), width = 10, height = 3)
par(mar = c(5,5,2,2), mfrow = c(1,3))
plot(yX_JY$DNAngUsed, yX_JY$h.idx, cex.axis = 2, cex.lab = 2, xlab = "DNA", ylab = "entropy", pch = 20, cex = 3, col = alpha("black", .5))
fig_label("a", region="figure", pos="topleft", cex = 3, font = 2)
plot(yX_JY$VCN, yX_JY$h.idx, cex.axis = 2, cex.lab = 2, xlab = "VCN", ylab = "entropy", pch = 20, cex = 3, col = alpha("black", .5))
fig_label("b", region="figure", pos="topleft", cex = 3, font = 2)
plot(yX_JY$seqDepth, yX_JY$h.idx, cex.axis = 2, cex.lab = 2, xlab = "SD", ylab = "entropy", pch = 20, cex = 3, col = alpha("black", .5))
fig_label("c", region="figure", pos="topleft", cex = 3, font = 2)
dev.off()







