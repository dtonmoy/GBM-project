
library(ggplot2)
library(genefilter)
library(data.table)
library(pheatmap)
library(GO.db)
library(org.Hs.eg.db)
library(openxlsx)
library(viridis)
library(gridExtra)
library(limma)
library(msigdbr)

# BUILD DB
m_df <- msigdbr(species = "Homo sapiens")
m_df$gs_subcat <- gsub("\\:", "_", m_df$gs_subcat)
m_df$gs_subcat <- paste0(m_df$gs_cat, "_", m_df$gs_subcat)
m_df$gs_subcat <- gsub("_$", "", m_df$gs_subcat)

subcat.unique <- sort(unique(m_df$gs_subcat))

dbList <- lapply(subcat.unique, function(i){
  m_df.sub <- m_df[m_df$gs_subcat == i, ]
  terms <- unique(m_df.sub$gs_name)
  db <- mclapply(terms, function(j){
    return(m_df.sub$entrez_gene[m_df.sub$gs_name == j])
  }, mc.cores = 8)
  names(db) <- terms
  return(db)
})
names(dbList) <- subcat.unique



############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################

entrez2symbol <- function(entrez)
{
  symbol <- mget(as.character(entrez), org.Hs.egSYMBOL, ifnotfound=NA)
  symbol <- unlist(lapply(symbol, function(i) return(i[1])))
  return(symbol)
}

symbol2entrez <- function(symbol)
{
  entrez <- mget(as.character(symbol), org.Hs.egSYMBOL2EG, ifnotfound=NA)
  entrez <- unlist(lapply(entrez, function(i) return(i[1])))
  entrez <- unique(entrez[!is.na(entrez)])
  return(entrez)
}

ensembl2entrez <- function(ensembl)
{
  entrez <- mget(as.character(ensembl), org.Hs.egENSEMBL2EG, ifnotfound=NA)
  entrez <- lapply(entrez, function(i) return(i[1]))
  return(unlist(entrez))
}

hyperG <- function(geneSets,DEgenes,universe, org.library, cutoff=0.1, mincount=2, parallel=T,adj.P.Val = F){
  library(org.library, character.only = T)
  library(foreach)
  library(doMC)
  if(parallel){
    registerDoMC(cores=detectCores())
    cores=detectCores()
  }else{
    cores=1
  }
  results <- mclapply(1:length(geneSets), function(i){
    results <- matrix(data=NA,ncol=9,nrow = 1)
    colnames(results) <- c('Term','Count','Size.univ', 'Size.tot','p-value','adj.P.Val','odds ratio','Entrez','Symbol')
    geneSet <- intersect(universe, geneSets[[i]])
    a <- length(intersect(DEgenes,geneSet))
    b <- length(setdiff(DEgenes,intersect(DEgenes,geneSet)))
    c <- length(setdiff(geneSet,intersect(DEgenes,geneSet)))
    d <- length(setdiff(universe,DEgenes)) - c
    contigency.matrix <- cbind(c(a,b),c(c,d))
    res <- fisher.test(contigency.matrix,alternative = 'greater')
    results[1,'Term'] <- names(geneSets)[i]
    results[1,'Count'] <- a
    results[1,'Size.univ'] <- length(intersect(geneSets[[i]], universe))
    results[1,'Size.tot'] <- length(geneSets[[i]])
    results[1,'p-value'] <- res$p.value
    results[1,'odds ratio'] <- res$estimate[[1]]
    # find genes annotated in the consensus term
    if(a > 0){
      genes <- intersect(DEgenes,geneSet)
      eid <- genes
      eid <- eid[order(eid)]
      results[1,'Entrez'] <- paste(eid,collapse="|")
    }
    return(results)
  }, mc.cores=cores)
  
  results <- as.data.frame(do.call(rbind, results))
  for(i in c(2, 3, 4, 5)){
    results[, i] <- as.numeric(as.character(results[, i]))
  }
  
  if(nrow(results) != 1){
    results <- results[order(results[,'p-value'],decreasing = FALSE),]
    results[,'adj.P.Val'] <- p.adjust(results[,'p-value'], 'BH')
    if(adj.P.Val){
      results <- as.data.frame(subset(results,results[,'adj.P.Val']<=cutoff))
    }else{
      results <- as.data.frame(subset(results,results[,'p-value']<=cutoff))
    }
    results <- as.data.frame(subset(results,results[,'Count']>=mincount))
  }else results <- as.data.frame(results)
  
  org.symb <- gsub(".db", "SYMBOL", org.library)
  # find genes 
  results$Symbol <- sapply(results$Entrez, function(x){
    y <- unlist(strsplit(as.character(x), "|", fixed=T))
    syms <- paste(unlist(mget(y, eval(parse(text=org.symb)),ifnotfound = NA)), collapse="|")
  })
  return(results)
}


my.read_xlsx <- function(inFile)
{
  mysheets <- getSheetNames(inFile)
  mList <- lapply(mysheets, read.xlsx, xlsxFile = inFile)
  names(mList) <- mysheets
  return(mList)
}


getHeatmapMatrix <- function(gsea_list, adjusted = FALSE)
{
  trms <- sort(unique(unlist(lapply(gsea_list, function(i) i$Term))))
  m <- matrix(NA, nrow = length(trms), ncol = length(gsea_list))
  
  for(i in 1:length(gsea_list))
  {
    gs <- gsea_list[[i]]
    idx <- match(gs$Term, trms)
    if(adjusted)m[idx, i] <- as.numeric(gs$"adj.P.Val")		
    else(m[idx, i] <- as.numeric(gs$"p-value"))
  }
  m[is.na(m)] <- 1	
  rownames(m) <- trms
  colnames(m) <- names(gsea_list)
  return(m)
}



fisherHeatmapSingle <- function(fh_list, outFile, nb, adjusted = FALSE, kw = NULL, gs = NULL)
{
  gsMat <- getHeatmapMatrix(fh_list, adjusted)	
  colnames(gsMat) <- names(fh_list)	
  
  if(!is.null(kw)) gsMat <- gsMat[grepl(kw, rownames(gsMat), ignore.case = TRUE), ]
  if(!is.null(gs)) gsMat <- gsMat[match(intersect(rownames(gsMat), gs), rownames(gsMat)), ]
  
  # select top X gene-sets per column
  idxMat <- apply(gsMat, 2, order)
  if(nrow(idxMat) > nb) idxMat <- idxMat[1:nb, ]
  gsMat <- gsMat[unique(as.numeric(idxMat)), ]	
  
  # Merge UP and DOWN columns
  idx.up <- grep("UP$", colnames(gsMat))
  idx.down <- grep("DOWN$", colnames(gsMat))
  gsMat.up <- gsMat[, idx.up]
  gsMat.down <- gsMat[, idx.down]
  
  gsMat <- matrix(NA, nrow = nrow(gsMat.up), ncol = ncol(gsMat.up))
  for(rowIdx in 1:nrow(gsMat)){
    for(colIdx in 1:ncol(gsMat)){
      pv.up <- gsMat.up[rowIdx, colIdx] 
      pv.down <- gsMat.down[rowIdx, colIdx]
      if(pv.up <= pv.down) gsMat[rowIdx, colIdx] <- -log10(pv.up)
      if(pv.up > pv.down) gsMat[rowIdx, colIdx] <- log10(pv.down)
    }
  }
  rownames(gsMat) <- rownames(gsMat.up)
  colnames(gsMat) <- gsub("\\.UP", "", colnames(gsMat.up))
  
  gsMat[gsMat > 6] <- 6 # set the limit to 5
  gsMat[gsMat < -6] <- -6
  
  paletteLength <- 25
  myMax <- ceiling(max(abs(gsMat)))
  
  myBreaks <- seq(-myMax , myMax, length.out=paletteLength)
  myBreaks <- myBreaks[myBreaks != 0]
  myColor <- colorRampPalette(c("deepskyblue", "snow", "orangered"))(paletteLength-2)
  
  doClust <- ifelse(nrow(gsMat) > 1, TRUE, FALSE)
  
  mysize <- nrow(gsMat) * 15 / 40
  mysize <- max(c(mysize, 10))
  if(mysize > 30) mysize <- 30
  
  pheatmap::pheatmap(gsMat, color = myColor, breaks = myBreaks, filename = outFile,
                     annotation_row = NULL, annotation_col = NULL,
                     cluster_cols = TRUE, cluster_rows = doClust, show_rownames = TRUE,
                     cellwidth = 10, cellheight = 10, fontsize_row = 8, fontsize_col = 8,
                     width = mysize, height = mysize)
}







############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################


# PARAMETERS
limmaDir <- file.path("../eisaR")

fisherDir <- file.path("../Fisher/")
dir.create(fisherDir, recursive = TRUE, showWarnings = FALSE)


###################
# FISHER FROM LIMMA

# PLEASE RUN eisa.R BEFORE !

setwd(limmaDir)
limmaFiles <- list.files(pattern = "_edgeR.xlsx")
names(limmaFiles) <- gsub("_edgeR.xlsx", "", limmaFiles)
limmaList <- lapply(limmaFiles, read.xlsx, sheet = 1)

genes <- limmaList[[1]]$entrez
genes <- unique(genes[!is.na(genes)])

pvCutoff <- 0.05
fcCutoff <- 0
limmaList.UP <- lapply(limmaList, function(i) i$entrez[i$FDR < pvCutoff & i$logFC > fcCutoff])
limmaList.UP <- lapply(limmaList.UP, function(i) unique(i[!is.na(i)]))
limmaList.DOWN <- lapply(limmaList, function(i) i$entrez[i$FDR < pvCutoff & i$logFC < -fcCutoff])
limmaList.DOWN <- lapply(limmaList.DOWN, function(i) unique(i[!is.na(i)]))
limmaList.DEG <- lapply(limmaList, function(i) i$entrez[i$FDR < pvCutoff & abs(i$logFC) > fcCutoff])
limmaList.DEG <- lapply(limmaList.DEG, function(i) unique(i[!is.na(i)]))


# Fisher
setwd(fisherDir)
lapply(1:length(dbList), function(i){
  mydb <- dbList[[i]]
  mydb.name <- names(dbList)[i]
  dir.create(mydb.name, showWarnings = FALSE)
  
  fh.UP <- lapply(limmaList.UP, hyperG, geneSets = mydb, universe = genes, org.library = "org.Hs.eg.db", cutoff = 1, mincount = 2)
  fh.DOWN <- lapply(limmaList.DOWN, hyperG, geneSets = mydb, universe = genes, org.library = "org.Hs.eg.db", cutoff = 1, mincount = 2)
  fh.DEG <- lapply(limmaList.DEG, hyperG, geneSets = mydb, universe = genes, org.library = "org.Hs.eg.db", cutoff = 1, mincount = 2)
  
  # Save
  for(i in 1:length(fh.UP))
  {
    fhName <- paste(names(fh.UP)[i], "_pv", pvCutoff, "_fc", fcCutoff, "_", mydb.name, "_Fisher.xlsx", sep = "")
    write.xlsx(list(UP = fh.UP[[i]], DOWN = fh.DOWN[[i]], DEG = fh.DEG[[i]]), file.path(mydb.name, fhName),
               row.names = FALSE, firstRow = T, headerStyle = createStyle(textDecoration = 'bold'))	
  }
})



##########################################################################################
# HEATMAPS

comp <- names(limmaList)

setwd(fisherDir)

lapply(1:length(dbList), function(i){
  mydb.name <- names(dbList)[i]
  
  # load fisher results
  fhFiles <- file.path(mydb.name, paste0(comp, "_pv", pvCutoff, "_fc", fcCutoff, "_", mydb.name, "_Fisher.xlsx"))
  
  fhList <- lapply(fhFiles, my.read_xlsx)
  names(fhList) <- gsub(paste0("_pv", pvCutoff, "_fc", fcCutoff, "_", mydb.name, "_Fisher.xlsx"), "", fhFiles)
  fhList <- do.call(c, fhList)
  fhList <- fhList[-grep("DEG", names(fhList))]
  
  # plot heatmap
  fisherHeatmapSingle(fhList, file.path(mydb.name, paste0(mydb.name, "_Fisher_heatmap.pdf")), nb = 5, adjusted = FALSE)
})











