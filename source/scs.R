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
source(paste(sourcePath,"scs-fun.R" , sep = ""))

cat("\nload data")
metadata <- read.csv(paste(StatsFolder, "/metadata.csv", sep = ""), sep = ";", dec = ".", header = T, row.names = NULL, stringsAsFactors = FALSE)
rownames(metadata) <- metadata$CompleteAmplificationID
ISmat <- read.csv(paste(isMatrices, "/ISmat_noMate_FE.tsv", sep = ""), sep = "\t", header = T, row.names = NULL, check.names = FALSE) 

rownames(ISmat) <- apply(ISmat[,1:5], 1, function(r){paste(r, collapse = "_")})
ISmat[is.na(ISmat)] <- 0
ISmat <- ISmat[, -c(1:5)]
ISmat <- round(as.matrix(ISmat))
ISmat <- Matrix(as.matrix(ISmat), sparse = TRUE)
ISmat <- ISmat[rowSums(ISmat) > 0,]
matrix(colnames(ISmat))
sum(!(colnames(ISmat) %in% metadata$CompleteAmplificationID))

metadata$h.idx<- NA
h.idx <- apply(ISmat, 2, function(s){
  s <- s[s>0];
  p <- s/sum(s);
  h <- -sum(p*log(p));
  return(h)
})
metadata$seqDepth<- NA
seq.depth <- apply(ISmat, 2, function(s){
  s <- s[s>0];
  sd <- sum(s)
  return(sd)
})
metadata$nIS<- NA
nIS <- apply(ISmat, 2, function(s){
  nis <- sum(s>0);
  return(nis)
})


metadata[names(h.idx), "h.idx"] <- h.idx
metadata[names(seq.depth), "seqDepth"] <- seq.depth 
metadata[names(nIS), "nIS"] <- nIS 

metadata <- metadata[intersect(rownames(metadata), colnames(ISmat)),]
metadata$mix <- as.vector(unlist(lapply(rownames(metadata), function(r){as.vector(unlist(strsplits(r, splits = c("-","_"), fixed = T)))[9]})))
JY_mixes <- c("jy10", "jy0", "jy1")

yX <- metadata[, c("DNAngUsed", "VCN", "seqDepth", "h.idx", "mix", "Kapa")]
yX_JY <- yX[which(yX$mix %in% JY_mixes),]

cat("\nSCS rescaling...")
## confounding factors:
x.L <- list(DNAngUsed = as.numeric(yX_JY[,"DNAngUsed"]),
            VCN = as.numeric(yX_JY[,"VCN"]),
            seqDepth = as.numeric(yX_JY[,"seqDepth"]))
nKTS.L <- list(kts.1 = 2, kts.2 = 2, kts.3 = 2)
sp.ord.L <- list(sp.ord.1 = 2, sp.ord.2 = 2, sp.ord.3 = 2)
confounders.L <- names(x.L)


res <- scs(x.L = x.L, 
           y = as.numeric(yX_JY$h.idx), 
           nKTS.L = nKTS.L, 
           sp.ord.L = sp.ord.L, 
           confounders.L = confounders.L, 
           RELTOL = 1e-8, 
           graphic = T)


yX_JY$h.idx.res <- NA
yX_JY$h.idx.res <- as.numeric(res$residuals)

write.csv(yX_JY, paste(currOutPath, "/VA-JYs.csv", sep = ""), row.names = TRUE)

# allAICRes2 <- apply(head(expand.grid(c(T,F), c(T,F), c(T,F)), -1), 1, 
#                     FUN = function(model){
#                       scs(x.L = x.L[t(model)], 
#                           y = as.numeric(yX_JY$h.idx), 
#                           nKTS.L = nKTS.L[t(model)], 
#                           sp.ord.L = sp.ord.L[t(model)], 
#                           confounders.L = confounders.L, 
#                           RELTOL = 1e-8, 
#                           graphic = F)
#                     })
# bestAIC <- order(as.vector(unlist(lapply(allAICRes2, function(l){l$AIC}))), decreasing = F)[1]
# bestBIC <- order(as.vector(unlist(lapply(allAICRes2, function(l){l$BIC}))), decreasing = F)[1]




