## name: geneset_functional_enrichment.r
## date: 03/20/2018

path1="../sub2/prediction/breast/"
path2="../sub2/prediction/ova/"

fold=5
bg=bi=bf=NULL
og=oi=of=NULL
for(i in 1:fold){
    tmp=read.delim(paste0(path1,"rna/cor_nrmse_b",i,".txt"),header=F,row.names=1)
    bg=c(bg,tmp[,1])
    tmp=read.delim(paste0(path1,"global_individual/cor_nrmse_b",i,".txt"),header=F,row.names=1)
    bi=c(bi,tmp[,1])
    tmp=read.delim(paste0(path1,"final/cor_nrmse_b",i,".txt"),header=F,row.names=1)
    bf=c(bf,tmp[,1])

    tmp=read.delim(paste0(path2,"rna/cor_nrmse_o",i,".txt"),header=F,row.names=1)
    og=c(og,tmp[,1])
    tmp=read.delim(paste0(path2,"global_individual/cor_nrmse_o",i,".txt"),header=F,row.names=1)
    oi=c(oi,tmp[,1])
    tmp=read.delim(paste0(path2,"final/cor_nrmse_o",i,".txt"),header=F,row.names=1)
    of=c(of,tmp[,1])
}

bg=apply(matrix(bg,ncol=fold),1,mean,na.rm=T)
bi=apply(matrix(bi,ncol=fold),1,mean,na.rm=T)
bf=apply(matrix(bf,ncol=fold),1,mean,na.rm=T)
dbf=bf-bg

og=apply(matrix(og,ncol=fold),1,mean,na.rm=T)
oi=apply(matrix(oi,ncol=fold),1,mean,na.rm=T)
of=apply(matrix(of,ncol=fold),1,mean,na.rm=T)
dof=of-og

tmp=read.delim(paste0(path1,"rna/cor_nrmse_b1.txt"),header=F,row.names=1)
bname=rownames(tmp)
tmp=read.delim(paste0(path2,"rna/cor_nrmse_o1.txt"),header=F,row.names=1)
oname=rownames(tmp)

names(bg)=names(bi)=names(bf)=names(dbf)=bname
names(og)=names(oi)=names(of)=names(dof)=oname

# 1. final
bg=bg[!is.na(bg)];bi=bi[!is.na(bi)];bf=bf[!is.na(bf)]
og=og[!is.na(og)];oi=oi[!is.na(oi)];of=of[!is.na(of)]

ind1=round(quantile(1:length(bf)),0)
ind2=round(quantile(1:length(of)),0)

for(i in 1:4){
    tmp=sort(bf)
    writeLines(names(tmp)[ind1[i]:ind1[i+1]],con=paste0("geneset/bf",i,".txt"))
    tmp=sort(of)
    writeLines(names(tmp)[ind2[i]:ind2[i+1]],con=paste0("geneset/of",i,".txt"))
}

# 2. final - global
dbf=dbf[!is.na(dbf)&dbf>0]
dof=dof[!is.na(dof)&dof>0]

ind1=round(quantile(1:length(dbf)),0)
ind2=round(quantile(1:length(dof)),0)

for(i in 1:4){
    tmp=sort(dbf)
    writeLines(names(tmp)[ind1[i]:ind1[i+1]],con=paste0("geneset/dbf",i,".txt"))
    tmp=sort(dof)
    writeLines(names(tmp)[ind2[i]:ind2[i+1]],con=paste0("geneset/dof",i,".txt"))
}


# 3. feature importance
mat=as.matrix(read.delim(paste0(path1,"feature_importance/fi.txt"),row.names=1,check.names=F))
bfi=apply(mat,2,sum)
mat=as.matrix(read.delim(paste0(path2,"feature_importance/fi.txt"),row.names=1,check.names=F))
ofi=apply(mat,2,sum)

# top genes above random selections
bfi=bfi[bfi>1]
ofi=ofi[ofi>1]

ind1=round(quantile(1:length(bfi)),0)
ind2=round(quantile(1:length(ofi)),0)

for(i in 1:4){
    tmp=sort(bfi)
    writeLines(names(tmp)[ind1[i]:ind1[i+1]],con=paste0("geneset/bfi",i,".txt"))
    tmp=sort(ofi)
    writeLines(names(tmp)[ind2[i]:ind2[i+1]],con=paste0("geneset/ofi",i,".txt"))
}

