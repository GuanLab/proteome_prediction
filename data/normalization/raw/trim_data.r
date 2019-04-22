## name: trim_data.r
## date: 07/07/2017

## Here I trim the raw dataset into trimmed_set which only contain samples with proteome, rna-seq and CNA results
## The sample with rna-seq results also have CNA results (copy number alteration)

#path1="/state3/hyangl/CPTAC_final/sub2/data/raw/"
#path2="/state3/hyangl/CPTAC_final/sub2/data/trimmed_set/"
path1="./";
path2="./";


############ ova ###############
proteome=as.matrix(read.delim(paste0(path1,"raw_ova_proteome.txt"),check.names=F,row.names=1))
rna=as.matrix(read.delim(paste0(path1,"raw_ova_rna.txt"),check.names=F,row.names=1))
cnv=as.matrix(read.delim(paste0(path1,"raw_ova_cnv.txt"),check.names=F,row.names=1))
dim(proteome);dim(rna);dim(cnv)
#[1] 7061   105
#[1] 15121   294
#[1] 11859   559


## rna&cnv
rna=rna[,colnames(rna) %in%  colnames(proteome)]
cnv=cnv[,colnames(cnv) %in%  colnames(proteome)]
dim(proteome);dim(rna);dim(cnv)
#[1] 7061  105
#[1] 15121   105
#[1] 11859   105

rna_new=proteome[,colnames(rna)];rna_new[]=NA
cnv_new=proteome[,colnames(cnv)];cnv_new[]=NA
ind=intersect(rownames(proteome),rownames(rna))
rna_new[ind,colnames(rna_new)]=rna[ind,colnames(rna_new)] # Here be carefull! Keep safe
ind=intersect(rownames(proteome),rownames(cnv))
cnv_new[ind,colnames(cnv_new)]=cnv[ind,colnames(cnv_new)]
dim(proteome);dim(rna_new);dim(cnv_new)
#[1] 7061  105
#[1] 7061  105
#[1] 7061  105

tmp=rna_new;tmp=cbind(rownames(tmp),tmp);colnames(tmp)[1]="Gene_ID"
write.table(tmp,file=paste0(path2,"trimmed_ova_rna.txt"),quote=F,sep="\t",row.names=F)
tmp=cnv_new;tmp=cbind(rownames(tmp),tmp);colnames(tmp)[1]="Gene_ID"
write.table(tmp,file=paste0(path2,"trimmed_ova_cnv.txt"),quote=F,sep="\t",row.names=F)
#######################################

############ breast ##################
proteome=as.matrix(read.delim(paste0(path1,"raw_breast_proteome.txt"),check.names=F,row.names=1))
rna=as.matrix(read.delim(paste0(path1,"raw_breast_rna.txt"),check.names=F,row.names=1))
cnv=as.matrix(read.delim(paste0(path1,"raw_breast_cnv.txt"),check.names=F,row.names=1))
dim(proteome);dim(rna);dim(cnv)
#[1] 10006   105
#[1] 15107    77
#[1] 16884    77

## rna&cnv
rna=rna[,colnames(rna) %in%  colnames(proteome)]
cnv=cnv[,colnames(cnv) %in%  colnames(proteome)]
dim(proteome);dim(rna);dim(cnv)
#[1] 10006   105
#[1] 15107    77
#[1] 16884    77

rna_new=proteome[,colnames(rna)];rna_new[]=NA
cnv_new=proteome[,colnames(cnv)];cnv_new[]=NA
ind=intersect(rownames(proteome),rownames(rna))
rna_new[ind,colnames(rna_new)]=rna[ind,colnames(rna_new)] # Here be carefull! Keep safe
ind=intersect(rownames(proteome),rownames(cnv))
cnv_new[ind,colnames(cnv_new)]=cnv[ind,colnames(cnv_new)]
dim(proteome);dim(rna_new);dim(cnv_new)
#[1] 10006   105
#[1] 10006    77
#[1] 10006    77

tmp=rna_new;tmp=cbind(rownames(tmp),tmp);colnames(tmp)[1]="Gene_ID"
write.table(tmp,file=paste0(path2,"trimmed_breast_rna.txt"),quote=F,sep="\t",row.names=F)
tmp=cnv_new;tmp=cbind(rownames(tmp),tmp);colnames(tmp)[1]="Gene_ID"
write.table(tmp,file=paste0(path2,"trimmed_breast_cnv.txt"),quote=F,sep="\t",row.names=F)
#######################################














