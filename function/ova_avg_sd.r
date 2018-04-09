## name: entry_avg_sd.r
## date: 10/17/2017

## Here I caculate ova proteome and phospho avg&sd of each gene/entry

path1="/state3/hyangl/CPTAC_round2/sub3/data/raw/"
path2="/state3/hyangl/CPTAC_round2/function/"

# proteome
x1=as.matrix(read.delim(paste0(path1,"retrospective_ova_PNNL_proteome_sort_common_gene_7061.txt"),check.names=F,row.names=1))
x2=as.matrix(read.delim(paste0(path1,"retrospective_ova_JHU_proteome_sort_common_gene_7061.txt"),check.names=F,row.names=1))
x=cbind(x1,x2)
mat=matrix(NA,nrow=dim(x)[1],ncol=2)
rownames(mat)=rownames(x)
mat[,1]=apply(x,1,mean,na.rm=T)
mat[,2]=apply(x,1,sd,na.rm=T)
write.table(mat,file=paste0(path2,"ova_proteome_avg_sd.txt"),quote=F,sep="\t",col.names=F)

x=as.matrix(read.delim(paste0(path1,"retrospective_breast_proteome_filtered.txt"),check.names=F,row.names=1))
mat=matrix(NA,nrow=dim(x)[1],ncol=2)
rownames(mat)=rownames(x)
mat[,1]=apply(x,1,mean,na.rm=T)
mat[,2]=apply(x,1,sd,na.rm=T)
write.table(mat,file=paste0(path2,"breast_proteome_avg_sd.txt"),quote=F,sep="\t",col.names=F)


# phospho
x=as.matrix(read.delim(paste0(path1,"retrospective_ova_phospho_sort_common_gene_10057.txt"),check.names=F,row.names=1))
mat=matrix(NA,nrow=dim(x)[1],ncol=2)
rownames(mat)=rownames(x)
mat[,1]=apply(x,1,mean,na.rm=T)
mat[,2]=apply(x,1,sd,na.rm=T)
write.table(mat,file=paste0(path2,"ova_phospho_avg_sd.txt"),quote=F,sep="\t",col.names=F)

x=as.matrix(read.delim(paste0(path1,"retrospective_breast_phospho_filtered.txt"),check.names=F,row.names=1))
mat=matrix(NA,nrow=dim(x)[1],ncol=2)
rownames(mat)=rownames(x)
mat[,1]=apply(x,1,mean,na.rm=T)
mat[,2]=apply(x,1,sd,na.rm=T)
write.table(mat,file=paste0(path2,"breast_phospho_avg_sd.txt"),quote=F,sep="\t",col.names=F)




