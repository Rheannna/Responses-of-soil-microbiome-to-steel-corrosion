#Random Forest for OTU and weight_loss
rm(list=ls())
options(stringsAsFactors = F)  
library(rfPermute)
library(randomForest)

load(file = "./different_otu.Rdata")
load(file = "../design_subset.Rdata")

#
if(F) {
  se <-read.table("../Fasting_map_all_beta_0223.txt",sep="\t",header=T,row.names=1)
  
  env <- se[,c("loss_ratio","redox")]
  
  design <- read.table("./design31.txt", header=T, row.names= 1, sep="\t")
  index_all = cbind(design[,c(3,5)], env[match(rownames(design), rownames(env)), ])
  data_all = cbind(design[,c(3,5)],t(otu_table_subset)[match(rownames(design), rownames(t(otu_table_subset))), ])
  idx1 <- design$Description%in% c("Q235s2","X80Cus2","Ps","Q235s3","X80Cus3")
  # 
  env1 <- index_all[idx1, ]
  
  data <- data_all[idx1,]
  data <- as.data.frame(t(data[,-c(1,2)]))
  
  #
  list = Reduce(intersect,  list(v1 = row.names(X80Cus3_enriched),
                                 v2 = row.names(X80Cus2_enriched),
                                 v3 = row.names(Q235s2_enriched),
                                 v4 = row.names(Q235s3_enriched),
                                 v5 = row.names(Ps_enriched)) )
  #
  union = Reduce(union,  list(v1 = row.names(X80Cus3_enriched),
                              v2 = row.names(X80Cus2_enriched),
                              v3 = row.names(Q235s2_enriched),
                              v4 = row.names(Q235s3_enriched)))
  
  dat <- data[union, ]
  #k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Caulobacterales;f__Caulobacteraceae;g__Phenylobacterium
  dat <- dat[-23, ]
 # dat <- log2(dat+1)
}
#
Dat.fill <- cbind(env1,t(dat)[match(row.names(env1),row.names(t(dat))), ])
####random forest for relation bwteen OTU and weight_loss

# #
# Dat.fill<-na.roughfix(index_abd_enrich)

##################RandomForest#####################
#rf1$importance
#rf1$predicted
#Mean Decrease Accuracy

# 

# loss_ratio
if(F){
  
  n <- length(names(Dat.fill[,c(5:ncol(Dat.fill))]))
  
  v <- c(0)
  set.seed(100)
  for (i in 1:(n-1)){
    mtry_fit <- randomForest(Dat.fill[,c(5:ncol(Dat.fill))],Dat.fill$loss_ratio,mtry = i)
    error <- mean(mtry_fit$mse)
    v[i] = error
    
    print(error)
  }
  
  #n = 29 for weight_loss 
  #run weight_loss
  set.seed(100)
  rf1 = randomForest(Dat.fill[,c(5:ncol(Dat.fill))],Dat.fill$loss_ratio,importance=TRUE, nodesize=5,proximity= TRUE,ntree=299, mtry= 29)
  plot(rf1) #ntree
  print(rf1)
  #
  round(importance(rf1), 2)
  
  varImpPlot(rf1,main = "Feature importance measured by Random Forest",n.var = 14, bg = par("bg"),
             color = par("fg"), gcolor = par("fg"), lcolor = "gray",cex.lab=1,cex=0.8)
  
  str(rf1$proximity)
  
  a <- rf1$importance
  
  
  
  #0223
  rf_CRs <- rfPermute(
    Dat.fill[,c(5:ncol(Dat.fill))],Dat.fill$loss_ratio,ntree = 100, 
    na.action = na.omit, mtry = 29, nrep = 50, num.cores = 1)
  rp.importance(rf_CRs, scale = TRUE)
  plot(rp.importance(rf_CRs, scale = TRUE))
  
  #
  rf_otu <- as.data.frame(rp.importance(rf_CRs,scale = TRUE))
  
  names(rf_otu)
  
  #write.table(file="randomforest_loss_ratio_enrich_otu_in_Ms.txt", rf_otu,sep="\t")
  
  #
  rf_otu$pval1 <- ifelse(rf_otu$`%IncMSE.pval`<0.05, "A", "B")
  rf_otu$pval2 <- ifelse(rf_otu$`IncNodePurity.pval`<0.05, "A", "B")
  
  #
  rf_otu <- subset(rf_otu,rf_otu$`%IncMSE`>0 & rf_otu$`IncNodePurity`>0 )
  
  
  
  rf_otu$name <- row.names(rf_otu)
  #
  library(stringr)
  name=str_split(rf_otu$name,';',simplify = T)[,c(3,4,5,6)] 
  name <- as.data.frame(name)
  name$names <- str_c(name[,1], name[,4],seq(1:nrow(rf_otu)),sep = ";")
  rf_otu$name <- name$names
  #
  rf_otu <- rf_otu[order(rf_otu$`%IncMSE`,rf_otu$IncNodePurity,decreasing = TRUE),]
  names(rf_otu)[1] <- "IncMSE" 
  
  # importance1 <- subset(rf_otu, rf_otu$`%IncMSE.pval`<0.05)
  # importance2 <- subset(rf_otu, rf_otu$`IncNodePurity.pval`<0.05)
  #IncMSE
  p1 <- ggplot(rf_otu, aes(x=reorder(name,IncMSE), y= IncMSE,fill = pval1)) + 
    geom_bar(stat = 'identity',color= "black") + 
    coord_flip() + labs(x="genus", y="IncMSE") + 
    ggtitle("IncMSE") + 
    theme(plot.title = element_text(size=18,face="bold"))
  p1
  
  ggsave(filename = "randomForest_loss_ratio_incMSE.pdf",p1,  width = 6, height = 5)
  
  
  #IncNodePurity
  p1 <- ggplot(rf_otu, aes(x=reorder(name,IncNodePurity), y= IncNodePurity,fill = pval2)) + 
    geom_bar(stat = 'identity',color= "black") + 
    coord_flip() + labs(x="genus", y="IncNodePurity") + 
    ggtitle("IncNodePurity") + 
    theme(plot.title = element_text(size=18,face="bold"))
  p1
  
  ggsave(filename = "randomForest_loss_ratio_IncNodePurity.pdf",p1,  width = 6, height = 5)
  
  
  #IncMSE
  p2 <- ggplot(data=rf_otu, aes(x=reorder(name,IncMSE), y=IncMSE, fill=pval1))+
    geom_bar(stat = "identity",position = "identity",width=0.9)+
    coord_flip()+
    scale_fill_manual(values = c("steelblue",'grey', "darkgreen"), guide=FALSE)+
    xlab("genus")+
    geom_abline(linetype="dashed",intercept = 50, slope = 0,size=1,colour='gray')+
    geom_abline(linetype="dashed",intercept = -50, slope = 0,size=1,colour='gray')+
    theme(panel.grid=element_blank(),panel.border=element_blank())
  
  p2 
  
}

#0229 redox
if(F){
  #
  n <- length(names(Dat.fill[,c(5:ncol(Dat.fill))]))
  
  v <- c(0)
  set.seed(100)
  for (i in 1:(n-1)){
    mtry_fit <- randomForest(Dat.fill[,c(5:ncol(Dat.fill))],Dat.fill$redox,mtry = i)
    error <- mean(mtry_fit$mse)
    v[i] = error
    
    print(error)
  }
  
  #n = 29 for redox
  #
  set.seed(100)
  rf1 = randomForest(Dat.fill[,c(5:ncol(Dat.fill))],Dat.fill$redox,importance=TRUE, nodesize=5,proximity= TRUE,ntree=299, mtry= 29)
  plot(rf1) #ntree
  print(rf1)
  #
  round(importance(rf1), 2)
  
  varImpPlot(rf1,main = "Feature importance measured by Random Forest",n.var = 14, bg = par("bg"),
             color = par("fg"), gcolor = par("fg"), lcolor = "gray",cex.lab=1,cex=0.8)
  
  str(rf1$proximity)
  
  a <- rf1$importance
  
  
  
  #0223
  rf_CRs <- rfPermute(
    Dat.fill[,c(5:ncol(Dat.fill))],Dat.fill$redox,ntree = 100, 
    na.action = na.omit, mtry = 29, nrep = 50, num.cores = 1)
  rp.importance(rf_CRs, scale = TRUE)
  plot(rp.importance(rf_CRs, scale = TRUE))
  
  #
  rf_otu <- as.data.frame(rp.importance(rf_CRs,scale = TRUE))
  
  names(rf_otu)
  
  #write.table(file="randomForest_redox_enrich_otu_in_Ms.txt", rf_otu, sep="\t")
  #
  rf_otu$pval1 <- ifelse(rf_otu$`%IncMSE.pval`<0.05, "A", "B")
  rf_otu$pval2 <- ifelse(rf_otu$`IncNodePurity.pval`<0.05, "A", "B")
  
  #
  rf_otu <- subset(rf_otu,rf_otu$`%IncMSE`>0 & rf_otu$`IncNodePurity`>0 )
  
  
  
  rf_otu$name <- row.names(rf_otu)
  #
  library(stringr)
  name=str_split(rf_otu$name,';',simplify = T)[,c(3,4,5,6)] #取行名，以';'号分割，提取。
  name <- as.data.frame(name)
  name$names <- str_c(name[,1], name[,4],seq(1:nrow(rf_otu)),sep = ";")
  rf_otu$name <- name$names
  #
  rf_otu <- rf_otu[order(rf_otu$`%IncMSE`,rf_otu$IncNodePurity,decreasing = TRUE),]
  names(rf_otu)[1] <- "IncMSE" 
  
  # importance1 <- subset(rf_otu, rf_otu$`%IncMSE.pval`<0.05)
  # importance2 <- subset(rf_otu, rf_otu$`IncNodePurity.pval`<0.05)
  #
  p1 <- ggplot(rf_otu, aes(x=reorder(name,IncMSE), y= IncMSE,fill = pval1)) + 
    geom_bar(stat = 'identity',color= "black") + 
    coord_flip() + labs(x="genus", y="IncMSE") + 
    ggtitle("IncMSE") + 
    theme(plot.title = element_text(size=18,face="bold"))
  p1
  ggsave(filename = "randomForest_redox_IncMSE.pdf",p1,  width = 6, height = 5)
  
  
  #IncNodePurity
  p1 <- ggplot(rf_otu, aes(x=reorder(name,IncNodePurity), y= IncNodePurity,fill = pval2)) + 
    geom_bar(stat = 'identity',color= "black") + 
    coord_flip() + labs(x="genus", y="IncNodePurity") + 
    ggtitle("IncNodePurity") + 
    theme(plot.title = element_text(size=18,face="bold"))
  p1
  ggsave(filename = "randomForest_redox_IncNodePurity.pdf",p1,  width = 6, height = 5)
  
  #
  p2 <- ggplot(data=rf_otu, aes(x=name, y=IncMSE, fill=pval1))+
    geom_bar(stat = "identity",position = "identity",width=0.9)+
    coord_flip()+
    scale_fill_manual(values = c("steelblue",'grey', "darkgreen"), guide=FALSE)+
    xlab("genus")+
    geom_abline(linetype="dashed",intercept = 50, slope = 0,size=1,colour='gray')+
    geom_abline(linetype="dashed",intercept = -50, slope = 0,size=1,colour='gray')+
    theme(panel.grid=element_blank(),panel.border=element_blank())
  
  p2 
  
}

#
obs <- aov(value ~ bgroup_soil, data = predicted2.melt)
#
Tukey_HSD_obs <- TukeyHSD(obs, ordered = FALSE, conf.level = 0.95)
# 
obs_table <- as.data.frame(Tukey_HSD_obs$bgroup_soil)
# 
obs_table #
