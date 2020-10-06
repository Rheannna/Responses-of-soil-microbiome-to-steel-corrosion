tropical=  c('darkorange', 'dodgerblue', 'red2', 'limegreen', 'yellow')
palette(tropical)


# #众数
# getmode <- function(v) {
#   uniqv <- unique(v)
#   uniqv[which.max(tabulate(match(v, uniqv)))]
# }
# 
# 
# library(RADanalysis)
# a <- matrix(nrow=3795,ncol=3)
# for(i in 1:3795) {a[i]<-length(which(otu_al20[i,]>0))} 
# row.names(a) <- row.names(otu_al20)
# 
# b <- t(otu_al20)
# for(i in 1:3795) {a[i,2] <- sum(b[,i])/a[i,1]} 
# #for(i in 1:3795) {a[i,3] <- mean(as.numeric(otu_al20[i,]),trim=0)} 
# colnames(a) <- c("number","means_reads","mean_local_abundance")
# h <- t(otu_norm)
# for(i in 1:3795) {a[i,3] <- sum(h[,i])/a[i,1]} 
# a <- as.data.frame(a)
# 
# a <- a[order(a$means_reads,decreasing=TRUE),]
# a[,4] <- row.names(a)
# 
# #c <- a$means_reads <800 for draw
# #d <- a[c,]
# m <- l$means_reads >5
# d <- l[m,]
# e <- a$number >15   #e <- a$number >100
# f <- a[e,]
# 
# a$number = as.factor(ifelse(a$number<60, "special",ifelse(a$means_reads>50, "general","nosig")))
# a$means_reads = as.factor(ifelse(a$means_reads>50, "high",ifelse(a$means_reads<40, "low","medium")))
# e <- a$number=="special"
# e <- a[e,]
# f <- e$means_reads=="high"
# f <- e[f,]
# special_abundant <- intersect(row.names(otu_al20),row.names(f))
# special_abundant <- otu_norm[special_abundant, ]
# special_abundant <- special_abundant*10000
# write.table(special_abundant,file="special_abundant.txt",sep="\t",col.names = NA)
#
###############
###############
###############
# #draw OTU plot
# pdf(file="otu1220.pdf", onefile=FALSE, paper="special", width=3, height=3, pointsize=6)
# plot(d$number, d$means_reads,xlab="sample numbers",ylab="OTU mean reads")
# points(d[rownames(core),]$number,d[rownames(core),]$means_reads,col=1,pch=16,cex = 1.2)
# points(d[rownames(abd),]$number,d[rownames(abd),]$means_reads,col=2,pch=16,cex = 1.2)
# points(d[setdiff(row.names(abd),row.names(core)),]$number,d[setdiff(row.names(abd),row.names(core)),]$means_reads,col=3,pch=16,cex = 1.2)
# points(d[intersect(row.names(core),row.names(ra)),]$number,d[intersect(row.names(core),row.names(ra)),]$means_reads,col=4,pch=16,cex = 1.2)
# text(d[setdiff(row.names(abd),row.names(core)),]$number,d[setdiff(row.names(abd),row.names(core)),]$means_reads,labels=d[setdiff(row.names(abd),row.names(core)),]$V4,cex=0.6,offset=-0.7,font=2,pos=1,col="gray50")
# legend("topleft",c("Core","Abundant","Speicialists"),col=1:3,pch=16)
# dev.off()
# 
# #draw CRs OTU plot
# CRs <- union(rownames(enriched_CRs2),rownames(enriched_CRs3))
# CRs_otu <- a[CRs,]   #OTU abundance and sample informations
# CRs_otu [,4] <- row.names(CRs_otu)
# pdf(file="otu_CRs_weight_loss.pdf", onefile=FALSE, paper="special", width=3, height=2, pointsize=6)
# plot(CRs_otu$number, CRs_otu$means_reads,xlab="Sample numbers",ylab="Mean reads of CRs OTUs ")
# points(CRs_otu[(setdiff(rownames(enriched_CRs2),CRs2_enriched_always)),]$number,CRs_otu[(setdiff(rownames(enriched_CRs2),CRs2_enriched_always)),]$means_reads,col=2,pch=16,cex = 1.2)
# points(CRs_otu[(setdiff(rownames(enriched_CRs3),CRs2_enriched_always)),]$number,CRs_otu[(setdiff(rownames(enriched_CRs3),CRs2_enriched_always)),]$means_reads,col=3,pch=16,cex = 1.2)
# points(CRs_otu[CRs2_enriched_always,]$number,CRs_otu[CRs2_enriched_always,]$means_reads,col=1,pch=16,cex = 1.2)
# CRs1 <- c("OTU_432","OTU_7","OTU_19","OTU_60","OTU_1050","OTU_79","OTU_285","OTU_389","OTU_151","OTU_4077","OTU_271","OTU_363","OTU_4884","OTU_2754","OTU_702","OTU_1379","OTU_765","OTU_482","OTU_1582","OTU_145","OTU_2137","OTU_258","OTU_201")
# points(CRs_otu[CRs1,]$number,CRs_otu[CRs1,]$means_reads,col=4,pch=16,cex = 1.2)
# text(CRs_otu[CRs1,]$number,CRs_otu[CRs1,]$means_reads,labels=CRs_otu[CRs1,]$V4,cex=0.3,offset=0.3,font=2,pos=3,col="gray30")
# legend("topleft",c("CRs_always","CRs2_only","CRs3_only","weight_loss assosiated"),col=1:4,pch=16,cex=0.6)
# dev.off()



design <- read.table("./20180308-201804 R result/beta/beta10000/design31.txt", header=T, row.names= 1, sep="\t")
alpha = read.table("./20180308-201804 R result/alpha/alpha_al20.txt", header=T, row.names= 1, sep="\t")
se <-read.table("./MENA/20181126/Fasting_map_all_beta_1126.txt",sep="\t",header=T,row.names=1)
env.data.log <- log1p(abs(se))
index = cbind(alpha, env.data.log[match(rownames(alpha), rownames(env.data.log)), ])

####environmental factors######
cc <- se[,c(1:8,12:15)]


# e <- a$means_rebundance >1
# f <- a[e,]
# 
# #core/rare/abundant/intermediate
# aa <- otu_norm[rownames(otu_al20),]
# aa <- t(aa)
# e <- a$number >80
# f <- a[e,]
# whole <- otu_norm[row.names(otu_al20),]
# core <- otu_norm[row.names(f), ]
# ra<-aa[,colSums(aa)/nrow(aa)<0.00005]
# ra <- as.data.frame(t(ra))
# abd<-aa[,colSums(aa)/nrow(aa)>0.0005]
# abd <- as.data.frame(t(abd))
# b<-aa[,colSums(aa)/nrow(aa)>=0.00005]
# int<-b[,colSums(b)/nrow(b)<=0.0005]

# count <- t(otu_al20) #OTU数 3795
# sum(colSums(count)) #2098376
# 
# 
# ra_count <- otu_al20[rownames(ra),]
# ra_count <- t(ra_count)
# sum(colSums(ra_count)) #103753 reads OTU 2226/3795
# 
# 
# abd_count <- otu_al20[rownames(abd),]
# abd_count <- t(abd_count)
# sum(colSums(abd_count)) #1614858 OTU 305/3795
# # >1% 589389 OTU 10/3795
# 
# int_count <- otu_al20[rownames(int),]
# int_count <- t(int_count)
# sum(colSums(int_count)) #379765 OTU 1264/3795
# 
# core_count <- otu_al20[rownames(core),]
# core_count <- t(core_count)
# sum(colSums(core_count)) #1776323 OTU850/3795
# 
# m <- a$mean_local_abundance >0.005
# n <- a[m,]
# local_abundant <- otu_al20[rownames(n),] #295
# 
# 
# s <- a$mean_local_abundance <0.0001
# u <- a[s,]
# local_rare <- otu_al20[rownames(u),] #1905
#########################

