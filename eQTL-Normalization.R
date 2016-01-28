##Eqtl-normalization steps

##Run the quantile normalization
library(preprocessCore)
require(bit64)
library(data.table)
library(peer) ##Only works in Linux for PEER
setwd("/home/rcf-proj/dn/eQTL/Expression")

#cere_r<-fread("cerebellum/gene_expression_r_intensity.txt",header=T)
#cere_g<-fread("cerebellum/gene_expression_g_intensity.txt",header=T)
#cere<-log2(cere_r[,-1,with=F]/cere_g[,-1,with=F])
cere<-fread("cerebellum/gene_expression.txt",header=T)


q.cere<-normalize.quantiles(as.matrix(cere[,-1,with=F]),copy=TRUE)


q.cere.rankN<-t(apply(q.cere,2,function(x)scale(rank(x))))
q.cere.rankN<-round(q.cere.rankN,5)
colnames(q.cere.rankN)<-colnames(cere)[-1]


q.cere.rankN<-q.cere.rankN[apply(q.cere,1,function(x)1-(sum(is.na(x))/length(x)))>.99,]



##Covars are PCAs and Age
model = PEER()
PEER_setPhenoMean(model,as.matrix(t(q.cere.rankN)))
PEER_setNk(model,15)
PEER_getNk(model)
PEER_update(model)
cere.peer.factors = PEER_getX(model)
cere.resid<-t(PEER_getResiduals(model))
#precision = PEER_getAlpha(model)

#png("cere.png")
#plot(1.0 / precision,xlab="Factors", ylab="Factor relevance", main="")
#dev.off()


##Extract values for use in analysis
ids<-as.character(colnames(cere)[-1])
peer.factors<-data.frame(ids,cere.peer.factors)


###Remove Duplicates
q.cere.rankN<-q.cere.rankN[,!duplicated(ids)]
q.cere.rankN<-cbind(cere[apply(q.cere,1,function(x)1-(sum(is.na(x))/length(x)))>.99,1,with=F],q.cere.rankN)
peer.factors<-peer.factors[!duplicated(ids)| !ids=="NA",]

cere.resid<-cere.resid[,!duplicated(ids)]
cere.resid<-cbind(cere[apply(q.cere,1,function(x)1-(sum(is.na(x))/length(x)))>.99,1,with=F],cere.resid)
colnames(cere.resid)[-1]<-ids[!duplicated(ids)]
cere.resid<-cbind(cere.resid[,1,with=F],round(cere.resid[,-1,with=F],6))
###Write Data
write.csv(peer.factors,"cere.peer.factors.csv",quote=F,row.names=F)
write.table(q.cere.rankN,"cere.rank.txt",quote=F,row.names=F)
write.table(cere.resid,"cere.resid.txt",quote=F,row.names=F)





#######PREFRONTAL
#pref_r<-fread("prefrontalcortex/gene_expression_r_intensity.txt",header=T)
#pref_g<-fread("prefrontalcortex/gene_expression_g_intensity.txt",header=T)
#pref<-log2(pref_r[,-1,with=F]/pref_g[,-1,with=F])
rm(list=ls())
pref<-fread("prefrontalcortex/gene_expression.txt",header=T)

q.pref<-normalize.quantiles(as.matrix(pref[,-1,with=F]),copy=TRUE)
q.pref.rankN<-t(apply(q.pref,2,function(x)scale(rank(x))))
q.pref.rankN<-round(q.pref.rankN,5)
colnames(q.pref.rankN)<-colnames(pref)[-1]


q.pref.rankN<-q.pref.rankN[apply(q.pref,1,function(x)1-(sum(is.na(x))/length(x)))>.99,]


model = PEER()
PEER_setPhenoMean(model,as.matrix(t(q.pref.rankN)))
PEER_setNk(model,15)
PEER_getNk(model)
PEER_update(model)
pref.peer.factors = PEER_getX(model)
pref.resid<-t(PEER_getResiduals(model))

ids<-as.character(colnames(pref)[-1])
peer.factors<-data.frame(ids,pref.peer.factors)

###Remove Duplicates
q.pref.rankN<-q.pref.rankN[,!duplicated(ids)]
q.pref.rankN<-cbind(pref[apply(q.pref,1,function(x)1-(sum(is.na(x))/length(x)))>.99,1,with=F],q.pref.rankN)
peer.factors<-peer.factors[!duplicated(ids) | !ids=="NA",]


pref.resid<-pref.resid[,!duplicated(ids)]
pref.resid<-cbind(pref[apply(q.pref,1,function(x)1-(sum(is.na(x))/length(x)))>.99,1,with=F],pref.resid)
colnames(pref.resid)[-1]<-ids[!duplicated(ids)]
pref.resid<-cbind(pref.resid[,1,with=F],round(pref.resid[,-1,with=F],6))


###Write Data
write.csv(peer.factors,"pref.peer.factors.csv",quote=F,row.names=F)
write.table(q.pref.rankN,"pref.rank.txt",quote=F,row.names=F)
write.table(pref.resid,"pref.resid.txt",quote=F,row.names=F)




#######vis
#vis_r<-fread("visualcortex/gene_expression_r_intensity.txt",header=T)
#vis_g<-fread("visualcortex/gene_expression_g_intensity.txt",header=T)
#vis<-log2(vis_r[,-1,with=F]/vis_g[,-1,with=F])
vis<-fread("visualcortex/gene_expression.txt",header=T)
q.vis<-normalize.quantiles(as.matrix(vis[,-1,with=F]),copy=TRUE)
q.vis.rankN<-t(apply(q.vis,1,function(x)scale(rank(x))))
q.vis.rankN<-round(q.vis.rankN,5)
colnames(q.vis.rankN)<-colnames(vis)[-1]

q.vis.rankN<-q.vis.rankN[apply(q.vis,1,function(x)1-(sum(is.na(x))/length(x)))>.99,]



model = PEER()
PEER_setPhenoMean(model,as.matrix(t(q.vis.rankN)))
PEER_setNk(model,15)
PEER_getNk(model)
PEER_update(model)
vis.peer.factors = PEER_getX(model)
vis.resid<-t(PEER_getResiduals(model))

ids<-as.character(colnames(vis)[-1])
peer.factors<-data.frame(ids,vis.peer.factors)


###Remove Duplicates
q.vis.rankN<-q.vis.rankN[,!duplicated(ids)]
q.vis.rankN<-cbind(vis[apply(q.vis,1,function(x)1-(sum(is.na(x))/length(x)))>.99,1,with=F],q.vis.rankN)
peer.factors<-peer.factors[!duplicated(ids)| !ids=="NA",]

vis.resid<-vis.resid[,!duplicated(ids)]
vis.resid<-cbind(vis[apply(q.vis,1,function(x)1-(sum(is.na(x))/length(x)))>.99,1,with=F],vis.resid)
colnames(vis.resid)[-1]<-ids[!duplicated(ids)]
vis.resid<-cbind(vis.resid[,1,with=F],round(vis.resid[,-1,with=F],6))


###Write Data
write.csv(peer.factors,"vis.peer.factors.csv",quote=F,row.names=F)
write.table(q.vis.rankN,"vis.rank.txt",quote=F,row.names=F)
write.table(vis.resid,"vis.resid.txt",quote=F,row.names=F)




