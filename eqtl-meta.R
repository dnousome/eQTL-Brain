###Harvard Rep MetaTissue
setwd("~/Desktop/HarvardReplication/")
vis<-fread("Expression/vis.resid.txt",header=T)
pref<-fread("Expression/pref.resid.txt",header=T)
cere<-fread("Expression/cere.resid.txt",header=T)


meta.annot<-rbindlist(final.annot)
meta.annot<-meta.annot[!duplicated(meta.annot$reporter_id),]



###Vis
covs = read.table('Pheno/finalcovsimputed.txt',header=T)
vis.trans.exp<-vis[vis$reporter_id %in% meta.annot$reporter_id,]
vis.trans.id<-vis.trans.exp[,1,with=F]

#peer<-read.csv("Expression/vis.peer.factors.csv")
#covs<-merge(covs,peer,by.x="IID",by.y="ids")
covs<-merge(covs,eigen,by.x="IID",by.y="V1")
covs$gender<-ifelse(covs$gender=="M",0,1)

temp.cov<-covs[covs$IID %in% colnames(vis.trans.exp),]
temp.cov<-temp.cov[temp.cov$IID %in% snp.final[,1],]
vis.trans.exp<-vis.trans.exp[,colnames(vis.trans.exp) %in% temp.cov$IID,with=F]
vis.trans.exp<-vis.trans.exp[,order(match(colnames(vis.trans.exp),temp.cov$IID)),with=F]

resid.vis<-data.frame(t(apply(vis.trans.exp,1,function(x)
  resid(lm(as.numeric(x)~temp.cov$gender+temp.cov$rin+temp.cov$V3+temp.cov$V4+temp.cov$V5))
)))
colnames(resid.vis)<-temp.cov$IID


##Cere
covs = read.table('Pheno/finalcovsimputed.txt',header=T)
cere.trans.exp<-cere[cere$reporter_id %in% meta.annot$reporter_id,]
cere.trans.id<-cere.trans.exp[,1,with=F]

#peer<-read.csv("Expression/cere.peer.factors.csv")
#covs<-merge(covs,peer,by.x="IID",by.y="ids")
covs<-merge(covs,eigen,by.x="IID",by.y="V1")
covs$gender<-ifelse(covs$gender=="M",0,1)

temp.cov<-covs[covs$IID %in% colnames(cere.trans.exp),]
temp.cov<-temp.cov[temp.cov$IID %in% snp.final[,1],]
cere.trans.exp<-cere.trans.exp[,colnames(cere.trans.exp) %in% temp.cov$IID,with=F]
cere.trans.exp<-cere.trans.exp[,order(match(colnames(cere.trans.exp),temp.cov$IID)),with=F]

resid.cere<-data.frame(t(apply(cere.trans.exp,1,function(x)
  resid(lm(as.numeric(x)~temp.cov$gender+temp.cov$rin+temp.cov$V3+temp.cov$V4+temp.cov$V5))
)))
colnames(resid.cere)<-temp.cov$IID


##Pref
covs = read.table('Pheno/finalcovsimputed.txt',header=T)
pref.trans.exp<-pref[pref$reporter_id %in% meta.annot$reporter_id,]
pref.trans.id<-pref.trans.exp[,1,with=F]

#peer<-read.csv("Expression/pref.peer.factors.csv")
#covs<-merge(covs,peer,by.x="IID",by.y="ids")
covs<-merge(covs,eigen,by.x="IID",by.y="V1")
covs$gender<-ifelse(covs$gender=="M",0,1)

temp.cov<-covs[covs$IID %in% colnames(pref.trans.exp),]
temp.cov<-temp.cov[temp.cov$IID %in% snp.final[,1],]
pref.trans.exp<-pref.trans.exp[,colnames(pref.trans.exp) %in% temp.cov$IID,with=F]
pref.trans.exp<-pref.trans.exp[,order(match(colnames(pref.trans.exp),temp.cov$IID)),with=F]

resid.pref<-data.frame(t(apply(pref.trans.exp,1,function(x)
  resid(lm(as.numeric(x)~temp.cov$gender+temp.cov$rin+temp.cov$V3+temp.cov$V4+temp.cov$V5))
)))
colnames(resid.pref)<-temp.cov$IID



write.table(resid.pref,"pref.txt",quote=F,sep="\t",row.names=F)
write.table(resid.cere,"cere.txt",quote=F,sep="\t",row.names=F)
write.table(resid.vis,"vis.txt",quote=F,sep="\t",row.names=F)

###IDs
covs = read.table('Pheno/finalcovsimputed.txt',header=T)
vis.ids<-data.frame(id=colnames(resid.vis),"vis")
pref.ids<-data.frame(id=colnames(resid.pref),"pref")
cere.ids<-data.frame(id=colnames(resid.cere),"cere")

id1<-merge(vis.ids,pref.ids,by="id",all=T)
id2<-merge(id1,cere.ids,by='id',all=T)
id2<-id2[id2$id %in% covs$IID,]
id2<-id2[order(match(id2$id,covs$IID)),]
id.fin<-data.frame(id2[,1],ifelse(is.na(id2[,-1]),0,1))
names(id.fin)<-c("#TISSUE","vis",'pref','cere')

write.table(id.fin,"tissue_info.txt",sep="\t",quote=F,row.names=F)

##SNP
meta.snp<-snp.final[snp.final$FID %in% id.fin[,1],]
meta.snp1<-t(meta.snp[,-1])
write.table(meta.snp1,"snp.meta.txt",sep='\t',row.names=F,col.names=F,quote=F)

ind<-data.frame(meta.snp[,1],"U","case")
write.table(ind,"ind.txt",quote=F,row.names=F,col.names=F)



##Probe
meta.annot.out<-data.frame(pref.trans.id$reporter_id,"chr1",1:length(pref.trans.id$reporter_id))
write.table(meta.annot.out,"probe.info.txt",row.names=F,quote=F,col.names=F)


##SNP
snps<-"rs1920116	G	A
rs2736100	G	T
rs2252586	G	A
rs11979158	A	G
rs4295627	A	G
rs147958197	T	C
rs55705857	A	G
rs72714236	G	A
rs72714295	C	A
rs72714302	G	C
rs72716319	A	G
rs72716328	C	T
rs4977756	A	G
rs1412829	T	C
rs498872	C	T
rs4809324	T	C
rs6010620	G	A"

temp.snp<-data.frame(matrix(unlist(strsplit(snps,"[\n\t]+")),ncol=3,byrow=T))
temp.snp<-temp.snp[order(match(temp.snp$X1,names(meta.snp)[-1])),]
temp.snp<-data.frame(temp.snp$X1,1,0,1:nrow(temp.snp),temp.snp$X2,temp.snp$X3)

write.table(temp.snp,"snp.txt",row.names=F,quote=F,col.names=F,sep='\t')




###Meta
source /usr/usc/java/1.7.0_65/setup.sh
java -jar /auto/rcf-proj/dn/eQTL/Meta-Tissue.v.0.4/MetaTissueInputGenerator.jar \
-i /auto/rcf-proj/dn/eQTL/Meta/tissue_info.txt \
-l /auto/rcf-proj/dn/eQTL/Meta/gene_list.txt \
-m /auto/rcf-proj/dn/eQTL/Meta/probe.info.txt \
-b /auto/rcf-proj/dn/eQTL/Meta/ind.txt \
-a /auto/rcf-proj/dn/eQTL/Meta/snp.meta.txt \
-c /auto/rcf-proj/dn/eQTL/Meta/snp.txt \
-p /auto/rcf-proj/dn/eQTL/Meta/2_MetaTissue_input/output_gene.txt \
-q /auto/rcf-proj/dn/eQTL/Meta/2_MetaTissue_input/output_snp.txt \
-r /auto/rcf-proj/dn/eQTL/Meta/2_MetaTissue_input/matrix.txt \
-v -x


/auto/rcf-proj/dn/eQTL/Meta-Tissue.v.0.4/MetaTissueMM \
--expr /auto/rcf-proj/dn/eQTL/Meta/2_MetaTissue_input/output_gene.txt \
--geno /auto/rcf-proj/dn/eQTL/Meta/2_MetaTissue_input/output_snp.txt \
--matrix /auto/rcf-proj/dn/eQTL/Meta/2_MetaTissue_input/matrix.txt \
--output /auto/rcf-proj/dn/eQTL/Meta/3_MetaTissue_output/MetaTissue \
--metatissue_bin_path /auto/rcf-proj/dn/eQTL/Meta-Tissue.v.0.4/Metasoft \
--dosage




##Once done use the final annot to paste to targetsnp and reporter ids
reporter_id<-lapply(final.annot,function(x)x[,1,with=F])
annot<-rbindlist(final.annot)[,-5:-6,with=F]
annot<-annot[!duplicated(annot$reporter_id),]
annot$reporter_id<-as.character(annot$reporter_id)

out<-list()
for (i in 1:length(targetSNP)){
  out[[i]]<-sapply(reporter_id[[i]],function(x)paste(targetSNP[i],x,sep=':'))
}
out<-unlist(out)


output<-read.table("Meta/MetaTissue.SNP.0.metasoft.output.txt",header=T)

new<-output[output$RSID %in% out,]
new$rs<-sapply(strsplit(as.character(new$RSID),":"),'[',1)
new$reporter_id<-sapply(strsplit(as.character(new$RSID),":"),'[',2)

final.new<-merge(new,annot,by="reporter_id")
final.new<-final.new[,c(1:15,17,19,26)]
final.new<-final.new[order(match(final.new$rs,targetSNP)),]
names(final.new)[11:13]<-c("VisM","PrefM","CereM")

new.list<-split(final.new,final.new$rs)
final.new1<-lapply(new.list,function(x){
  x$FinalP<-ifelse(x$PVALUE_Q<0.05,x$PVALUE_RE2,x$PVALUE_FE)
  x$FDR<-p.adjust(x$FinalP,method="fdr")
  x[order(x$FDR),]
})

setwd("~/Desktop/HarvardReplication/Meta/Results")
sapply(names(final.new1),function(x)write.table(final.new1[[x]],paste(x,"Meta.txt",sep="."),row.names=F,quote=F,sep="\t"))
