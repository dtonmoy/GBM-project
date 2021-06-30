
library(org.Hs.eg.db)
library(openxlsx)
library(pheatmap)
library(GO.db)
library(fgsea)
library(GeneAnswers)
library(viridis)
library(gridExtra)
library(msigdbr)
library(ggplot2)
library(ungeviz)
library(limma)
library(tidyr)


############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################

getKbLength <- function(entrez, organism = "hsa"){
  require(EDASeq)
  lgth <- getGeneLengthAndGCContent(entrez, organism, mode = "org.db")
  return(lgth[,1] / 1000)
}

getTPM <- function(countMat, organism = "hsa"){
  
  geneLength <- getKbLength(rownames(countMat), organism)
  
  # remove NA gene length
  toRemove <- is.na(geneLength)
  geneLength <- geneLength[!toRemove]
  countMat <-countMat[!toRemove, ]
  
  # get TPM
  x <- countMat / geneLength
  tpmMat <- t( t(x) * 1e6 / colSums(x) )
  
  return(tpmMat)
}

my.read_xlsx <- function(inFile)
{
  mysheets <- getSheetNames(inFile)
  mList <- lapply(mysheets, read.xlsx, xlsxFile = inFile)
  names(mList) <- mysheets
  return(mList)
}


fgsea_pipeline <- function(rank, pathways, minSize=15, maxSize=500, nperm=10000){
  fg <- fgsea(rank, pathways = pathways, minSize=minSize, maxSize=maxSize, nperm=nperm)
  fg <- as.data.frame(fg)
  fg.entrez <- fg$leadingEdge
  fg.symbol <- lapply(fg.entrez, entrez2symbol)
  fg.symbol <- unlist(lapply(fg.symbol, function(i) paste(i, collapse = ",")))
  avg.delta <- unlist(lapply(fg.entrez, function(i) mean(rank[i])))
  fg <- data.frame(fg, leadingEdge.symbol = fg.symbol, avgRank = avg.delta)
  fg <- fg[order(fg$pval), ]
  return(fg)
}


symbol2entrez <- function(symbol)
{
  entrez <- mget(as.character(symbol), org.Hs.egSYMBOL2EG, ifnotfound=NA)
  entrez <- unique(unlist(lapply(entrez, function(i) return(i[1]))))
  entrez <- entrez[!is.na(entrez)]
  return(entrez)
}

entrez2symbol <- function(entrez)
{
  symbol <- mget(as.character(entrez), org.Hs.egSYMBOL, ifnotfound=NA)
  symbol <- unique(unlist(lapply(symbol, function(i) return(i[1]))))
  return(symbol)
}

ensembl2entrez <- function(ensembl)
{
  entrez <- mget(as.character(ensembl), org.Hs.egENSEMBL2EG, ifnotfound=NA)
  entrez <- lapply(entrez, function(i) return(i[1]))
  return(unlist(entrez))
}

refseq2entrez <- function(refseq){
  entrez <- mget(as.character(refseq), org.Hs.egREFSEQ2EG, ifnotfound=NA)
  entrez <- unique(unlist(lapply(entrez, function(i) return(i[1]))))
  return(entrez)
}



getHeatmapMatrix <- function(gsea_list, scoreName)
{
  trms <- sort(unique(unlist(lapply(gsea_list, function(i) i[,1]))))
  m <- matrix(NA, nrow = length(trms), ncol = length(gsea_list))
  
  for(i in 1:length(gsea_list))
  {
    gs <- gsea_list[[i]]
    idx <- match(unlist(gs[,1]), trms)
    m[idx, i] <- as.numeric(gs[, scoreName])
  }
  m[is.na(m)] <- 0	
  rownames(m) <- trms
  colnames(m) <- names(gsea_list)
  return(m)
}

getHeatmapMatrix.pv <- function(gsea_list, adjusted = FALSE)
{
  trms <- sort(unique(unlist(lapply(gsea_list, function(i) i[,1]))))
  m <- matrix(NA, nrow = length(trms), ncol = length(gsea_list))
  
  for(i in 1:length(gsea_list))
  {
    gs <- gsea_list[[i]]
    idx <- match(unlist(gs[,1]), trms)
    if(adjusted)m[idx, i] <- as.numeric(gs$"pval")		
    else(m[idx, i] <- as.numeric(gs$"padj"))
  }
  m[is.na(m)] <- 1	
  rownames(m) <- trms
  colnames(m) <- names(gsea_list)
  return(m)
}










############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################

# BUILD GENE-SETS FROM FILE
mygs.raw <- read.xlsx(file.path("~/Research/Sajib/P10_GBM/doc/Gene_Set_ALA_POS.xlsx"), sheet = 1)
mygs <- lapply(1:ncol(mygs.raw), function(i) as.character(mygs.raw[, i]))
mygs <- lapply(mygs, function(i) i[!is.na(i)])
mygs <- lapply(mygs, symbol2entrez)
names(mygs) <- colnames(mygs.raw)

dbList <- list()
dbList[["leadingEdge"]] <- mygs


# PARAMETERS
fgseaDir <- file.path("../fgsea_GLASS_single")
dir.create(fgseaDir, recursive = TRUE, showWarnings = FALSE)

rankedGenesList <- readRDS(file.path("../data/survival/GLASS_single_rankedGenesList.rds"))

ann.sample <- readRDS(file.path("../data/survival/GLASS_clinical.rds"))


#########################
# Perform fgsea analysis

dbList <- dbList[c("leadingEdge")]

# FGSEA
setwd(fgseaDir)
lapply(1:length(dbList), function(i){
  mydb <- dbList[[i]]
  mydb.name <- names(dbList)[i]
  
  plotDir <- fgseaDir   
  
  dir.create(file.path(plotDir, mydb.name), showWarnings = FALSE)
  setwd(file.path(plotDir, mydb.name))   
  
  #########################
  # Perform fgsea analysis
  
  fgseaResList <- lapply(rankedGenesList,
                         fgsea, pathways = mydb, minSize=5, maxSize=1000
  )
  fgseaResList <- lapply(fgseaResList, function(i) i[order(i$pval),])	
  
  # Add symbols
  fgseaResList <- lapply(fgseaResList, as.data.frame)
  
  fgseaResList <- lapply(fgseaResList, function(k){
    k.symbol <- lapply(k$leadingEdge, entrez2symbol)
    k.symbol <- lapply(k.symbol, function(j) unique(j[!is.na(j)]))
    k.symbol <- lapply(k.symbol, toString)
    k$leadingEdge.symbol <- unlist(k.symbol)
    k
  })
  
  fgseaResList <- lapply(fgseaResList, function(j){
    j$NES[is.na(j$NES)] <- 0
    return(j)
  })
  
  # Divide UP and DOWN
  fgseaResList.UP <- lapply(fgseaResList, function(j) j[j$NES > 0, ])
  names(fgseaResList.UP) <- paste0(names(fgseaResList), ".UP")
  fgseaResList.DOWN <- lapply(fgseaResList, function(j) j[j$NES < 0, ])
  names(fgseaResList.DOWN) <- paste0(names(fgseaResList), ".DOWN")
  fgseaResList <- c(fgseaResList.UP, fgseaResList.DOWN)
  
  # Save fgsea output   
  lapply(1:length(rankedGenesList), function(j)
    write.xlsx(list(UP =  fgseaResList.UP[[j]], DOWN = fgseaResList.DOWN[[j]]),
               paste(names(rankedGenesList)[j], "_", mydb.name, "_fgsea.xlsx", sep = ""),
               row.names = FALSE, firstRow = T, headerStyle = createStyle(textDecoration = 'bold'))
  )	
  
})


#############
# NES HEATMAP

setwd(fgseaDir)

comp <- names(rankedGenesList)

lapply(1:length(dbList), function(i){
  mydb.name <- names(dbList)[i]
  
  # load fisher results
  fhFiles <- file.path(mydb.name, paste0(comp, "_", mydb.name, "_fgsea.xlsx"))
  
  fhList <- lapply(fhFiles, function(j){
    up <- read.xlsx(j, sheet = "UP")
    down <- read.xlsx(j, sheet = "DOWN")
    return(rbind(up, down))
  })
  
  names(fhList) <- gsub(paste0("_", mydb.name, "_fgsea.xlsx"), "", fhFiles)
  names(fhList) <- gsub(paste0(mydb.name, "\\/"), "", names(fhList))
  #fhList <- do.call(c, fhList)
  
  nesMat <- getHeatmapMatrix(fhList, "NES")
  write.xlsx(nesMat, file.path(mydb.name, paste0(mydb.name, "_NES.xlsx")), row.names = TRUE)
})

# AVERAGE NES PER SAMPLE
nesMat <- read.xlsx(file.path(fgseaDir, "leadingEdge/leadingEdge_NES.xlsx"), rowNames = TRUE)
nesAvg <- colMeans(nesMat[c("Inf..wound.Res.", "Mesenchymal"), ])
nesSum <- colSums(nesMat[c("Inf..wound.Res.", "Mesenchymal"), ])

##############
# Vs. survival

keep <- !(is.na(ann.sample$"Overall.survival.(months)"))
ann.sample <- ann.sample[keep, ]
nesSum <- nesSum[keep]


# PRIMARY
keep <- ann.sample$Status == "Primary tumor"

#cor.test(nesAvg[keep], ann.clin$days_to_death[keep], method = "spearman")
mycor <- cor.test(nesSum[keep], ann.sample$"Overall.survival.(months)"[keep], method = "spearman")

# PLOT
ggmat <- data.frame(Surv = ann.sample$"Overall.survival.(months)"[keep],
                    Score = nesSum[keep])
p <- ggplot(ggmat, aes(Surv, Score))
p <- p + geom_point(colour = "grey", size = 3)
p <- p + geom_smooth(method = "lm", se=F)
p <- p + xlab("Overall survival (month)") + ylab("Score")
p <- p + ggtitle(paste0("R= ", round(mycor$estimate, digits = 2), "; pvalue= ", signif(mycor$p.value, digits = 2)))
p <- p + theme_bw(base_size = 16)

ggsave(plot = p, filename = "Primary_score_VS_survival.pdf", width = 6, height = 6)


# RECURRENT
keep <- ann.sample$Status == "Recurrent"

#cor.test(nesAvg[keep], ann.clin$days_to_death[keep], method = "spearman")
mycor <- cor.test(nesSum[keep], ann.sample$"Overall.survival.(months)"[keep], method = "spearman")

# PLOT
ggmat <- data.frame(Surv = ann.sample$"Overall.survival.(months)"[keep],
                    Score = nesSum[keep])
p <- ggplot(ggmat, aes(Surv, Score))
p <- p + geom_point(colour = "grey", size = 3)
p <- p + geom_smooth(method = "lm", se=F)
p <- p + xlab("Overall survival (month)") + ylab("Score")
p <- p + ggtitle(paste0("R= ", round(mycor$estimate, digits = 2), "; pvalue= ", signif(mycor$p.value, digits = 2)))
p <- p + theme_bw(base_size = 16)

ggsave(plot = p, filename = "Recurrent_score_VS_survival.pdf", width = 6, height = 6)


# EXCLUDE OUTLIER
keep <- ann.sample$"Overall.survival.(months)" < 100
ann.sample <- ann.sample[keep, ]
nesSum <- nesSum[keep]



# PRIMARY
keep <- ann.sample$Status == "Primary tumor"

#cor.test(nesAvg[keep], ann.clin$days_to_death[keep], method = "spearman")
mycor <- cor.test(nesSum[keep], ann.sample$"Overall.survival.(months)"[keep], method = "spearman")

# PLOT
ggmat <- data.frame(Surv = ann.sample$"Overall.survival.(months)"[keep],
                    Score = nesSum[keep])
p <- ggplot(ggmat, aes(Surv, Score))
p <- p + geom_point(colour = "grey", size = 3)
p <- p + geom_smooth(method = "lm", se=F)
p <- p + xlab("Overall survival (month)") + ylab("Score")
p <- p + ggtitle(paste0("R= ", round(mycor$estimate, digits = 2), "; pvalue= ", signif(mycor$p.value, digits = 2)))
p <- p + theme_bw(base_size = 16)

ggsave(plot = p, filename = "Primary_score_VS_survival_woOutlier.pdf", width = 6, height = 6)


# RECURRENT
keep <- ann.sample$Status == "Recurrent"

#cor.test(nesAvg[keep], ann.clin$days_to_death[keep], method = "spearman")
mycor <- cor.test(nesSum[keep], ann.sample$"Overall.survival.(months)"[keep], method = "spearman")

# PLOT
ggmat <- data.frame(Surv = ann.sample$"Overall.survival.(months)"[keep],
                    Score = nesSum[keep])
p <- ggplot(ggmat, aes(Surv, Score))
p <- p + geom_point(colour = "grey", size = 3)
p <- p + geom_smooth(method = "lm", se=F)
p <- p + xlab("Overall survival (month)") + ylab("Score")
p <- p + ggtitle(paste0("R= ", round(mycor$estimate, digits = 2), "; pvalue= ", signif(mycor$p.value, digits = 2)))
p <- p + theme_bw(base_size = 16)

ggsave(plot = p, filename = "Recurrent_score_VS_survival_woOutlier.pdf", width = 6, height = 6)


