##Reanalysis for HarvardData
##ssh drrayl@68.181.127.129 
##Grep genotypes for all targetSNPs
setwd("~/Desktop/HarvardReplication/")
library(data.table)
require(bit64)

##Create all file for plink###
geno<-fread("genotype.txt")
geno1<-t(geno)
write.table(geno1,"Geno1",quote=F,row.names=F,col.names=F)

sed 's/NA/0\/0/g' Geno1 >Geno2
sed 's/\// /g' Geno2 >Geno3
tail -n +2 Geno3 > Geno4

mv Geno4 harvard.ped

map<-read.table("features.txt",header=T)
finmap<-data.frame(chr=map[,2],rsid=map[,4],0,pos=map[,3])
write.table(finmap,"harvard.map",quote=F,row.names=F,col.names=F)

plink --file harvard --no-fid --no-parents --no-sex  --no-pheno --allow-extra-chr --out harvard_bin

plink --bfile harvard_bin --pca --allow-extra-chr --out harvard_pca --filter Pheno.final.txt --out harvard_pca_ctrl
eigen<-read.table("~/Desktop/HarvardReplication/PC/harvard_pca.eigenvec")

plot(eigen$V5,eigen$V6)
eigenval<-read.table("~/Desktop/HarvardReplication/PC/harvard_pca.eigenval")
sum(eigenval$V1[1:5])/sum(eigenval$V1)
###First 3 explain 22% of the variance


###Use the strand to convert
./update_build.sh harvard_bin HumanHap650Yv3_A-b37.strand harvard_bin_impute

##Chr8 128363940-132365228

plink --bfile harvard_bin --chr 8 --from-bp 128363940 --to-bp 132365228 --geno .99 -mind .99 --recode --out chr8_qc --allow-extra-chr
plink --bfile harvard_bin_impute --chr 8 --from-bp 128363940 --to-bp 132365228 --geno .99 -mind .99 --recode --out chr8_toimpute

### chr8
shapeit --input-ped chr8_toimpute.ped chr8_toimpute.map -T 8 \
-M /auto/rcf-proj/dn/nousome/annotations/1000G/genetic_map_chr8_combined_b37.txt \
-O chr8_phase



impute2 \
-h $ANNO/1000G/chr8.haps.gz \
-l $ANNO/1000G/chr8.legend.gz \
-known_haps_g chr8_phase.haps \
-int 128363940 132365228 \
-buffer 500 -Ne 20000 \
-phase \
-m $ANNO/1000G/genetic_map_chr8_combined_b37.txt \
-o results_chr8.gen \
-iter 30 \
-k 80 \
-allow_large_regions \
-k_hap 1038 \
-i results_chr8.info \
-burnin 10 -o_gz -r results_chr8.summary


zcat results_chr8.gen.gz| grep "rs55705857\|rs72714236\|rs72714295\|rs72714302\|rs72716319\|rs72716328\|rs147958197" >chr8_impute.gen
grep "rs55705857\|rs72714236\|rs72714295\|rs72714302\|rs72716319\|rs72716328\|rs147958197" results_chr8.info >chr8_impute.info



###Read.pheno
pheno<-data.frame(fread("~/Desktop/HarvardReplication/Pheno/phenotype.txt",header=T))
pheno.d<-data.frame(fread("~/Desktop/HarvardReplication/Pheno/phenotype_description.txt"))
pheno.f<-merge(pheno.d,pheno,by.x="trait_id")

pheno.f1<-pheno.f[pheno.f$trait_name %in% c('age_death','age_onset','ph','gender','pmi','batch','rin','alz_status','huntingtons_status'),]
traitnames<-pheno.f1[,2]
pheno.f2<-as.data.frame(t(pheno.f1[,-1:-7]))
names(pheno.f2)<-traitnames
pheno.f2$ID_1<-rownames(pheno.f2)
pheno.f2$IID<-substring(rownames(pheno.f2),2)

table(pheno.f2$alz_status,pheno.f2$huntingtons_status)

pheno.final<-pheno.f2[,c("IID","IID",'ID_1','age_death','age_onset','ph','gender','pmi','batch','rin','alz_status','huntingtons_status')]
pheno.ctrl.final<-pheno.final[pheno.final$alz_status=="control" & pheno.final$huntingtons_status=="control",]
write.table(pheno.ctrl.final,"pheno.ctrl.txt",quote=F,row.names=F)


##Extract SNP
targetSNP<-("rs1920116
            rs2736100
            rs2252586
            rs11979158
            rs4295627
            rs147958197
            rs55705857
            rs72714236
            rs72714295
            rs72714302
            rs72716319
            rs72716328
            rs4977756
            rs1412829
            rs498872
            rs4809324
            rs6010620")
targetSNP<-unlist(strsplit(targetSNP,"\n"))


write.table(targetSNP,"targetSNPs.txt",row.names=F,quote=F,col.names=F)

###Run in plink
plink -bfile harvard_bin --allow-extra-chr --extract targetSNPs.txt --recodeA --out harvard-targetSNPs



##
pheno.ctrl.final<-merge(pheno.ctrl.final,eigen[,-2],by.x="IID",by.y="V1")
pheno.final<-merge(pheno.final,eigen[,-2],by.x="IID",by.y="V1")
write.table(pheno.ctrl.final,"Pheno.final.CTRL.txt",quote=F,row.names=F)
write.table(pheno.final,"Pheno.final.all.txt",quote=F,row.names=F)





###Run back in R in MAC
library(RMySQL)
library(data.table)
require(bit64)


###Read in SNP Data
snp<-read.table("Genotype/harvard-targetSNPs.raw",header=T)
snp<-snp[snp$FID %in% covs$IID,]
snp<-snp[order(match(snp$FID,covs$IID)),]
snp<-snp[,-2:-6]
names(snp)<-sapply(strsplit(names(snp),"_"),'[',1)


##Chr8
chr8.dosage<-read.table("Genotype/chr8_impute.gen",header=F)

chr8.dosage1<-data.frame(apply(chr8.dosage[,-1:-5],1,function(x) 
  (0*x[seq(1,length(x),by=3)])+(1*x[seq(2,length(x),by=3)])+(2*x[seq(3,length(x),by=3)])
))
names(chr8.dosage1)<-chr8.dosage[,2]
chr8.sample<-read.table("genotype/chr8_phase.sample",header=T)
chr8.sample<-chr8.sample[-1,]
chr8<-data.frame(id=chr8.sample[,1],chr8.dosage1)


snp.final<-merge(snp,chr8,by.x="FID",by.y="id",all.x=T)
snp.final<-data.frame(FID=snp.final[,1],snp.final[,-1][,order(match(names(snp.final)[-1],targetSNP))])






setwd("~/Desktop/HarvardReplication/")
covs <- read.table('Pheno/pheno.ctrl.txt',header=T)


##Run Imputetation

vis.ids<-read.table("~/Desktop/HarvardReplication/Expression/vis.resid.txt",nrows = 1)
pref.ids<-read.table("~/Desktop/HarvardReplication/Expression/pref.resid.txt",nrows = 1)
cere.ids<-read.table("~/Desktop/HarvardReplication/Expression/cere.resid.txt",nrows = 1)

vis.ids<-data.frame(t(vis.ids[1,-1]),"vis")
pref.ids<-data.frame(t(pref.ids[1,-1]),"pref")
pref.ids<-pref.ids[complete.cases(pref.ids),]
cere.ids<-data.frame(t(cere.ids[1,-1]),"cere")

id1<-merge(vis.ids,pref.ids,by="X1",all=T)
id2<-merge(id1,cere.ids,by='X1',all=T)
id2<-id2[id2$X1 %in% covs$IID,]
id2<-id2[order(match(id2$X1,covs$IID)),]
id.fin<-data.frame(id2[,1],ifelse(is.na(id2[,-1]),0,1))
names(id.fin)<-c("#TISSUE","vis",'pref','cere')

id.fin<-id.fin[id.fin[,1] %in% snp.final$FID,]

fin.covs<-covs[covs$IID %in% id.fin[,1],]

summary(fin.covs$age_death)
sd(fin.covs$age_death)

summary(fin.covs$rin)
sd(fin.covs$rin)



prop.table(table(fin.covs$gender))

library(mice)
library(VIM)
imp.mi <- mice(fin.covs[,c("age_death","ph","gender","pmi","batch","rin")],seed=1000)

final.covs<-cbind(IID=fin.covs[,1],complete(imp.mi,1))
write.table(final.covs,"Pheno/finalcovsimputed.txt",row.names=F,quote=F)

##MAF
mafs<-round(apply(snp.final[,-1],2,function(x)sum(x)/(length(x)*2)),4)
clip <- pipe("pbcopy", "w")                       
write.table(mafs, file=clip,quote=F,row.names=F,col.names=T)                               
close(clip)



##MYSQL QUERY

mychannel <- dbConnect(MySQL(), user="genome", host="genome-mysql.cse.ucsc.edu")
query <- function(...) dbGetQuery(mychannel, ...)

snploc<-rbindlist(lapply(targetSNP,function(x)query(paste0("SELECT chrom, chromStart, chromEnd FROM hg19.snp142 WHERE name='",x,"';"))))
snploc$start<-ifelse(snploc$chrom=="chr8",snploc$chromStart-1900000,snploc$chromStart-1000000)
snploc$end<-snploc$chromEnd+1000000



gencode<-apply(snploc,1,function(x)query(paste0("SELECT name2 FROM hg19.wgEncodeGencodeBasicV19 WHERE chrom='",x[1],"' AND txStart >",x[4]," AND txEnd<",x[5],";")))
refseq<-apply(snploc,1,function(x)query(paste0("SELECT name FROM hg19.refGene WHERE chrom='",x[1],"' AND txStart >",x[4]," AND txEnd<",x[5],";")))


annot<-fread('Expression/reporter.txt',header=T)
annot.new<-fread("Expression/GPL4372.annot",skip=27)
annot.gen<-lapply(gencode,function(x)annot[unique(unlist(sapply(x$name2,function(x)grep(paste0("\\b",x,"\\b"),annot$gene_name)))),])
annot.ref<-lapply(refseq,function(x)annot[unique(unlist(sapply(x$name,function(x)grep(paste0("\\b",x,"\\b"),annot$substance_id)))),])


final.annot<-mapply(function(x,y){
  temp<-rbind(x,y)
  temp<-temp[!duplicated(temp$reporter_id),]
  temp<-merge(temp,annot.new,by.x="reporter_id",by.y="ID")
  loc<-sapply(strsplit(temp$`Chromosome location`,"[pqcen-]+"),'[',1)
  tabloc<-table(loc)
  if (length(tabloc)>1){
    temp[loc %in% c(names(tabloc[order(-tabloc)][1]),NA),]
  }else{
    temp
  }
}, annot.gen, annot.ref, SIMPLIFY=FALSE)





##eQTL for Cere

library(ggplot2)
setwd("~/Desktop/HarvardReplication/")
eigen<-read.table("~/Desktop/HarvardReplication/PC/harvard_pca.eigenvec")
covs = read.table('Pheno/finalcovsimputed.txt',header=T)


#peer<-read.csv("Expression/cere.peer.factors.csv")

#covs<-merge(covs,peer,by.x="IID",by.y="ids")
covs<-merge(covs,eigen,by.x="IID",by.y="V1")

covs$gender<-ifelse(covs$gender=="M",0,1)

covs.analysis<-covs[covs$IID %in% snp.final$FID,]
covs.analysis<-covs.analysis[order(match(covs.analysis$IID,snp.final$FID)),]

cere<-fread("Expression/cere.resid.txt",header=T)


results<-list()
results.trans<-list()

for (i in 1:length(final.annot)){
  temp.exp<-cere[cere$reporter_id %in% final.annot[[i]]$reporter_id,]
  temp.exp1<-temp.exp[,-1,with=F]
  temp.exp1<-temp.exp1[,colnames(temp.exp1) %in% covs.analysis$IID,with=F]
  temp.exp2<-temp.exp1[,order(match(colnames(temp.exp1),covs.analysis$IID)),with=F]
  
  temp.cov<-covs.analysis[covs.analysis$IID %in% colnames(temp.exp2),]
  temp.cov<-temp.cov[order(match(temp.cov$IID,colnames(temp.exp2))),]
  
  temp.snp<-snp.final[snp.final$FID %in% temp.cov$IID,]
  temp.snp<-temp.snp[order(match(temp.snp$FID,temp.cov$IID)),]
  
  #temp.peer<-peer[peer$ids %in% temp.cov$IID,]
  #temp.peer<-temp.peer[order(match(temp.peer$ids,temp.cov$IID)),]
  
  
  ###For Cis
  reg<-apply(temp.exp2,1,function(x){
    summary(lm(as.numeric(x)~temp.snp[,i+1]+temp.cov$gender+temp.cov$rin
               +temp.cov$V3+temp.cov$V4+temp.cov$V5))
  })
  
  
  temp.res<-data.frame(t(sapply(reg,function(x)x$coefficients[2,c(1,4)])),reporter_id=temp.exp$reporter_id)
  temp.res$fdr<-p.adjust(temp.res$Pr...t..,method="fdr")
  temp.res<-merge(temp.res,final.annot[[i]],by="reporter_id")
  temp.res<-temp.res[,c(1,2,3,4,5,7)]
  temp.res1<-temp.res[order(temp.res$fdr),]
  
  names(temp.res1)<-c("id",'Beta',"p","fdr","refseq","gene")
  temp.res1<-temp.res1[,c('id','refseq','gene','Beta','p','fdr')]
  temp.res1[,4:6]<-round(temp.res1[,4:6],3)
  results[[i]]<-temp.res1
  
  
  
  
  ####Trans Analysis
  temp.trans<-cere[!cere$reporter_id %in% final.annot[[i]]$reporter_id,]
  temp.trans1<-temp.trans[,-1,with=F]
  temp.trans1<-temp.trans1[,colnames(temp.trans1) %in% covs.analysis$IID,with=F]
  temp.trans2<-temp.trans1[,order(match(colnames(temp.trans1),covs.analysis$IID)),with=F]
  
  reg.trans<-sapply(apply(temp.trans2,1,function(x){
    summary(lm(as.numeric(x)~temp.snp[,i+1]+temp.cov$gender+temp.cov$rin+temp.cov$V3+temp.cov$V4+temp.cov$V5+.,
               data=temp.cov[,c(13:27)]))}),
    function(x)x$coefficients[2,c(1,4)])
  
  temp.res<-data.frame(t(reg.trans),reporter_id=temp.trans$reporter_id)
  temp.res$fdr<-p.adjust(temp.res$Pr...t..,method="fdr")
  temp.res<-merge(temp.res,annot,by="reporter_id")
  temp.res<-temp.res[,c(1,2,3,4,5,7)]
  temp.res1<-temp.res[order(temp.res$fdr),]
  names(temp.res1)<-c("id",'Beta',"p","fdr","refseq","gene")
  temp.res1<-temp.res1[,c('id','refseq','gene','Beta','p','fdr')]
  results.trans[[i]]<-temp.res1
  
  
  
  ###Plots for cis
  setwd("~/Desktop/HarvardReplication/Cere")
  dir.create(colnames(temp.snp)[-1][i])
  setwd(colnames(temp.snp)[-1][i])
  for (j in 1:nrow(temp.exp2)){
    plot.data<-data.frame(snp=round(temp.snp[,i+1]),rankZ=as.numeric(temp.exp2[j,]))
    png(paste0(sapply(strsplit(temp.res[j,6],","),'[',1),".png"))
    print(ggplot(data=plot.data,aes(x=factor(snp),y=rankZ))+ geom_boxplot())
    dev.off()
    
  }
}
names(results)<-targetSNP
setwd("~/Desktop/HarvardReplication/Cere")
sapply(names(results),function(x)write.table(results[[x]],paste(x,"txt",sep="."),row.names=F,quote=F))




###VISUAL
setwd("~/Desktop/HarvardReplication/")
vis<-fread("Expression/vis.rank.txt",header=T)
covs = read.table('Pheno/pheno.ctrl.txt',header=T)
peer<-read.csv("Expression/vis.peer.factors.csv")
covs<-merge(covs,peer,by.x="IID",by.y="ids")
covs<-merge(covs,eigen,by.x="IID",by.y="V1")
covs$gender<-ifelse(covs$gender=="M",0,1)
covs.analysis<-covs[covs$IID %in% snp.final$FID,]
covs.analysis<-covs.analysis[order(match(covs.analysis$IID,snp.final$FID)),]




results<-list()
results.trans<-list()
for (i in 1:length(final.annot)){
  temp.exp<-vis[vis$reporter_id %in% final.annot[[i]]$reporter_id,]
  temp.exp1<-temp.exp[,-1,with=F]
  temp.exp1<-temp.exp1[,colnames(temp.exp1) %in% covs.analysis$IID,with=F]
  temp.exp2<-temp.exp1[,order(match(colnames(temp.exp1),covs.analysis$IID)),with=F]
  
  temp.cov<-covs.analysis[covs.analysis$IID %in% colnames(temp.exp2),]
  temp.cov<-temp.cov[order(match(temp.cov$IID,colnames(temp.exp2))),]
  
  temp.snp<-snp.final[snp.final$FID %in% temp.cov$IID,]
  temp.snp<-temp.snp[order(match(temp.snp$FID,temp.cov$IID)),]
  
  temp.peer<-peer[peer$ids %in% temp.cov$IID,]
  temp.peer<-temp.peer[order(match(temp.peer$ids,temp.cov$IID)),]
  
  reg<-apply(temp.exp2,1,function(x){
    summary(lm(as.numeric(x)~temp.snp[,i+1]+temp.cov$gender+temp.cov$rin+temp.cov$V3+temp.cov$V4+temp.cov$V5+.,
               data=temp.cov[,c(13:27)]))
  })
  temp.res<-data.frame(t(sapply(reg,function(x)x$coefficients[2,c(1,4)])),reporter_id=temp.exp$reporter_id)
  temp.res$fdr<-p.adjust(temp.res$Pr...t..,method="fdr")
  temp.res<-merge(temp.res,final.annot[[i]],by="reporter_id")
  temp.res<-temp.res[,c(1,2,3,4,5,7)]
  temp.res1<-temp.res[order(temp.res$fdr),]
  names(temp.res1)<-c("id",'Beta',"p","fdr","refseq","gene")
  temp.res1<-temp.res1[,c('id','refseq','gene','Beta','p','fdr')]
  temp.res1[,4:6]<-round(temp.res1[,4:6],3)
  results[[i]]<-temp.res1
  
  ########TRANS
  temp.trans<-vis[!vis$reporter_id %in% final.annot[[i]]$reporter_id,]
  temp.trans1<-temp.trans[,-1,with=F]
  temp.trans1<-temp.trans1[,colnames(temp.trans1) %in% covs.analysis$IID,with=F]
  temp.trans2<-temp.trans1[,order(match(colnames(temp.trans1),covs.analysis$IID)),with=F]
  
  reg.trans<-sapply(apply(temp.trans2,1,function(x){
    summary(lm(as.numeric(x)~temp.snp[,i+1]+temp.cov$gender+temp.cov$rin+temp.cov$V3+temp.cov$V4+temp.cov$V5+.,
               data=temp.cov[,c(13:27)]))}),
    function(x)x$coefficients[2,c(1,4)])
  
  temp.res<-data.frame(t(reg.trans),reporter_id=temp.trans$reporter_id)
  temp.res$fdr<-p.adjust(temp.res$Pr...t..,method="fdr")
  temp.res<-merge(temp.res,annot,by="reporter_id")
  temp.res<-temp.res[,c(1,2,3,4,5,7)]
  temp.res1<-temp.res[order(temp.res$fdr),]
  names(temp.res1)<-c("id",'Beta',"p","fdr","refseq","gene")
  temp.res1<-temp.res1[,c('id','refseq','gene','Beta','p','fdr')]
  results.trans[[i]]<-temp.res1
  
  
  
  setwd("~/Desktop/HarvardReplication/Vis")
  dir.create(colnames(temp.snp)[-1][i])
  setwd(colnames(temp.snp)[-1][i])
  for (j in 1:nrow(temp.exp2)){
    plot.data<-data.frame(snp=round(temp.snp[,i+1]),rankZ=as.numeric(temp.exp2[j,]))
    png(paste0(sapply(strsplit(temp.res[j,6],","),'[',1),".png"))
    print(ggplot(data=plot.data,aes(x=factor(snp),y=rankZ))+ geom_boxplot())
    dev.off()
  }
}

names(results)<-targetSNP
setwd("~/Desktop/HarvardReplication/Vis")
sapply(names(results),function(x)write.table(results[[x]],paste(x,"txt",sep="."),row.names=F,quote=F))







###PREFRONTAL
setwd("~/Desktop/HarvardReplication/")
pref<-fread("Expression/pref.rank.txt",header=T)
covs = read.table('Pheno/pheno.ctrl.txt',header=T)
peer<-read.csv("Expression/pref.peer.factors.csv")
covs<-merge(covs,peer,by.x="IID",by.y="ids")
covs<-merge(covs,eigen,by.x="IID",by.y="V1")
covs$gender<-ifelse(covs$gender=="M",0,1)
covs.analysis<-covs[covs$IID %in% snp.final$FID,]
covs.analysis<-covs.analysis[order(match(covs.analysis$IID,snp.final$FID)),]


results<-list()
results.trans<-list()
for (i in 1:length(final.annot)){
  temp.exp<-pref[pref$reporter_id %in% final.annot[[i]]$reporter_id,]
  temp.exp1<-temp.exp[,-1,with=F]
  temp.exp1<-temp.exp1[,colnames(temp.exp1) %in% covs.analysis$IID,with=F]
  temp.exp2<-temp.exp1[,order(match(colnames(temp.exp1),covs.analysis$IID)),with=F]
  
  temp.cov<-covs.analysis[covs.analysis$IID %in% colnames(temp.exp2),]
  temp.cov<-temp.cov[order(match(temp.cov$IID,colnames(temp.exp2))),]
  
  temp.snp<-snp.final[snp.final$FID %in% temp.cov$IID,]
  temp.snp<-temp.snp[order(match(temp.snp$FID,temp.cov$IID)),]
  
  
  reg<-apply(temp.exp2,1,function(x){
    summary(lm(as.numeric(x)~temp.snp[,i+1]+temp.cov$gender+temp.cov$rin+temp.cov$V3+temp.cov$V4+temp.cov$V5+.,
               data=temp.cov[,c(13:27)]))
  })
  temp.res<-data.frame(t(sapply(reg,function(x)x$coefficients[2,c(1,4)])),reporter_id=temp.exp$reporter_id)
  temp.res$fdr<-p.adjust(temp.res$Pr...t..,method="fdr")
  temp.res<-merge(temp.res,final.annot[[i]],by="reporter_id")
  temp.res<-temp.res[,c(1,2,3,4,5,7)]
  temp.res1<-temp.res[order(temp.res$fdr),]
  names(temp.res1)<-c("id",'Beta',"p","fdr","refseq","gene")
  temp.res1<-temp.res1[,c('id','refseq','gene','Beta','p','fdr')]
  temp.res1[,4:6]<-round(temp.res1[,4:6],3)
  results[[i]]<-temp.res1
  
  ###Trans
  reg.trans<-sapply(apply(temp.trans2,1,function(x){
    summary(lm(as.numeric(x)~temp.snp[,i+1]+temp.cov$gender+temp.cov$rin+temp.cov$V3+temp.cov$V4+temp.cov$V5+.,
               data=temp.cov[,c(13:27)]))}),
    function(x)x$coefficients[2,c(1,4)])
  
  temp.res<-data.frame(t(reg.trans),reporter_id=temp.trans$reporter_id)
  temp.res$fdr<-p.adjust(temp.res$Pr...t..,method="fdr")
  temp.res<-merge(temp.res,annot,by="reporter_id")
  temp.res<-temp.res[,c(1,2,3,4,5,7)]
  temp.res1<-temp.res[order(temp.res$fdr),]
  names(temp.res1)<-c("id",'Beta',"p","fdr","refseq","gene")
  temp.res1<-temp.res1[,c('id','refseq','gene','Beta','p','fdr')]
  results.trans[[i]]<-temp.res1
  
  
  
  setwd("~/Desktop/HarvardReplication/Pref")
  dir.create(colnames(temp.snp)[-1][i])
  setwd(colnames(temp.snp)[-1][i])
  for (j in 1:nrow(temp.exp2)){
    plot.data<-data.frame(snp=round(temp.snp[,i+1]),rankZ=as.numeric(temp.exp2[j,]))
    png(paste0(sapply(strsplit(temp.res[j,6],","),'[',1),".png"))
    print(ggplot(data=plot.data,aes(x=factor(snp),y=rankZ))+ geom_boxplot())
    dev.off()
  }
}

names(results)<-targetSNP
setwd("~/Desktop/HarvardReplication/Pref")
sapply(names(results),function(x)write.table(results[[x]],paste(x,"txt",sep="."),row.names=F,quote=F))



