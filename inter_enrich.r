

load("tmp1.RData")




library(DOSE)
library(org.Mm.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)

  allg <- rownames(dat_i$x)
  
  Allgenes = bitr(allg, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
  
  #GO
  DE_CU_ALL <- list()
  for(i in 1:length(table(res2_i$cluster))){
    CU_genes <- bitr(allg[which(res2_i$cluster==i)], fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
    DE_CU_ALL[[i]] <- enrichGO(gene = CU_genes$ENTREZID,OrgDb = org.Hs.eg.db,keyType = 'ENTREZID',
                          ont = "ALL",universe=Allgenes$ENTREZID,pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,
                          readable = TRUE)
    #DE_CU_ALL1 <- setReadable(DE_CU_ALL, OrgDb = org.Hs.eg.db)
    #DE_CU_ALL2 <- summary(DE_CU_ALL1)
    
    #write.csv(DE_CU_ALL2,file=paste("CU_",i,"_GO_enrich_ALL.csv",sep=""))
    cat("CU:",i,"\n")
    #write.csv(DE_CU_ALL2[which(DE_CU_ALL2$pvalue<0.05),],file=paste("CU_",i,"_GO_enrich_PVALUE_0.05.csv",sep=""))
    
  }
  DE_CU_KEGG <- list()
for(i in 1:length(table(res2_i$cluster))){
  
  CU_genes <- bitr(allg[which(res2_i$cluster==i)], fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
  DE_CU_KEGG[[i]] <- enrichKEGG(gene=CU_genes$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 1,
                           pAdjustMethod = "BH",universe = Allgenes$ENTREZID,minGSSize = 10, maxGSSize = 500,
                           qvalueCutoff = 1, use_internal_data = TRUE)
  cat("CU:",i,"\n")
  #write.csv(DE_CU_KEGG@result,file=paste("CU_",i,"_KEGG_enrich_ALL.csv",sep=""))
  
  #write.csv(DE_CU_KEGG@result[which(DE_CU_KEGG@result$pvalue<0.05),],file=paste("CU_",i,"_KEGG_enrich_PVALUE_0.05.csv",sep=""))
}



  GO_extract <- function(GO1,thre=0.0001,count=15){
    
    
    n <- length(GO1)
    
    allk <-c()
    for(i in 1:n){
      
      ex1 <- GO1[[i]]@result
      indd <- (ex1[,7]<thre)+(ex1$Count>count)
      ex11 <- ex1[which(indd==2),3]
      allk <- c(allk,ex11)
      
    }
    allk1 <- unique(allk)
    
    alld <- c()
    for(i in 1:n){
      ex1 <- GO1[[i]]@result
      alld <- cbind(alld,ex1[match(allk1,ex1[,3]),7])
    }
    
    alld[which(is.na(alld))] <- 1
    rownames(alld) <- allk1
    colnames(alld) <- 1:n
    
    return(alld)
    
  }
  
  KEGG_extract <- function(KEGG){
    
    
    n <- length(KEGG)
    
    allk <-c()
    for(i in 1:n){
      
      ex1 <- KEGG[[i]]@result
      if(is.null(ex1))
         ex11 <- ex1[which(ex1[,6]<0.05),2]
         allk <- c(allk,ex11)
         
    }
    allk1 <- unique(allk)
    
    alld <- c()
    for(i in 1:n){
      ex1 <- KEGG[[i]]@result
      alld <- cbind(alld,ex1[match(allk1,ex1[,2]),6])
    }
    
    alld[which(is.na(alld))] <- 1
    rownames(alld) <- allk1
    colnames(alld) <- 1:n
    
    return(alld)
    
  } 
  
  
  #ID to GENE SYMBOL
  
  GO <- list();
  for(i in 1:length(table(CU))){
    
    GO[[i]] <- setReadable(DE_CU_ALL[[i]], OrgDb = org.Hs.eg.db)
    
  }
  
  
  go1 <- GO_extract(GO,thre=0.001,count=20)
  
  
  
  
  library(pheatmap)
  library(RColorBrewer)
  
  
  pheatmap(go1,treeheight_row=0, treeheight_col=0,cluster_row = TRUE, cluster_col = FALSE,
           display_numbers = FALSE,colorRampPalette((brewer.pal(n=7,name ="RdYlBu")))(100),
           angle_col=0,filename = "GO.pdf",width=16,height=20,fontsize = 13,fontsize_row=17)
  
  
  
  
  
  
  
  
  
  
