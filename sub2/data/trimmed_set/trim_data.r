## name: trim_data.r
## date: 12/05/2017

## Here I trim the raw dataset into trimmed_set which only contain samples with proteome, rna-seq and CNA results
## The sample with rna-seq results also have CNA results (copy number alteration)

#path1="/state3/hyangl/CPTAC_exp/sub2/data/raw/"
#path2="/state3/hyangl/CPTAC_exp/sub2/data/trimmed_set/"
path1="../raw/"
path2="./"

############ ova ###############
pnnl=as.matrix(read.delim(paste0(path1,"retrospective_ova_PNNL_proteome_sort_common_gene_7061.txt"),check.names=F,row.names=1))
jhu=as.matrix(read.delim(paste0(path1,"retrospective_ova_JHU_proteome_sort_common_gene_7061.txt"),check.names=F,row.names=1))
rna=as.matrix(read.delim(paste0(path1,"retrospective_ova_rna_seq_sort_common_gene_15121.txt"),check.names=F,row.names=1))
cnv=as.matrix(read.delim(paste0(path1,"retrospective_ova_CNA_sort_common_gene_11859.txt"),check.names=F,row.names=1))
dim(pnnl);dim(jhu);dim(rna);dim(cnv)
#[1] 7061   84
#[1] 7061  122
#[1] 15121   294
#[1] 11859   559

# trim sample with rna/cnv results
pnnl=pnnl[,colnames(pnnl) %in%  colnames(rna)]
jhu_overlap=jhu[,colnames(jhu) %in% colnames(pnnl)] # 24 overlapped samples 
jhu=jhu[,colnames(jhu) %in%  colnames(rna)]
jhu=jhu[,!(colnames(jhu) %in% colnames(pnnl))] # exclude overlapped examples in jhu
proteome=cbind(pnnl,jhu) ## Here the order of samples changes when you combine two matrices!!
proteome=proteome[,sort(colnames(proteome))] # sort sample ids!!
dim(pnnl);dim(jhu);dim(jhu_overlap);dim(proteome)
#[1] 7061   54
#[1] 7061   51
#[1] 7061   24
#[1] 7061  105

# write trimmed dataset without imputation; 7061*54/51
tmp=pnnl;tmp=cbind(rownames(tmp),tmp);colnames(tmp)[1]="Gene_ID"
write.table(tmp,file=paste0(path2,"ova_proteome_pnnl.txt"),quote=F,sep="\t",row.names=F)
tmp=jhu;tmp=cbind(rownames(tmp),tmp);colnames(tmp)[1]="Gene_ID"
write.table(tmp,file=paste0(path2,"ova_proteome_jhu.txt"),quote=F,sep="\t",row.names=F)
tmp=jhu_overlap;tmp=cbind(rownames(tmp),tmp);colnames(tmp)[1]="Gene_ID"
write.table(tmp,file=paste0(path2,"ova_proteome_jhu_overlap.txt"),quote=F,sep="\t",row.names=F)
tmp=proteome;tmp=cbind(rownames(tmp),tmp);colnames(tmp)[1]="Gene_ID"
write.table(tmp,file=paste0(path2,"ova_proteome.txt"),quote=F,sep="\t",row.names=F)

## rna&cnv
proteome=cbind(pnnl,jhu) ## Here the order of samples changes when you combine two matrices!!
proteome=proteome[,sort(colnames(proteome))] # sort sample ids!!
dim(proteome) #[1] 7061  105

rna=rna[,colnames(rna) %in%  colnames(proteome)]
cnv=cnv[,colnames(cnv) %in%  colnames(proteome)]
dim(proteome);dim(rna);dim(cnv)
#[1] 7061  105
#[1] 15121   105
#[1] 11859   105

rna_new=proteome;rna_new[]=NA
cnv_new=proteome;cnv_new[]=NA
ind=intersect(rownames(proteome),rownames(rna))
rna_new[ind,colnames(rna)]=rna[ind,colnames(rna)] # Here be carefull! Keep safe
ind=intersect(rownames(proteome),rownames(cnv))
cnv_new[ind,colnames(cnv)]=cnv[ind,colnames(cnv)]
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
proteome=as.matrix(read.delim(paste0(path1,"retrospective_breast_proteome_sort_common_gene_10005.txt"),check.names=F,row.names=1))
rna=as.matrix(read.delim(paste0(path1,"retrospective_breast_RNA_sort_common_gene_15107.txt"),check.names=F,row.names=1))
cnv=as.matrix(read.delim(paste0(path1,"retrospective_breast_CNA_sort_common_gene_16884.txt"),check.names=F,row.names=1))
dim(proteome);dim(rna);dim(cnv)
#[1] 10006   105
#[1] 15107    77
#[1] 16884    77

# trim sample with rna/cnv results
proteome=proteome[,colnames(proteome) %in%  colnames(rna)]
dim(proteome) #[1] 10006    77

# write trimmed dataset without imputation; 10623*77
tmp=proteome;tmp=cbind(rownames(tmp),tmp);colnames(tmp)[1]="Gene_ID"
write.table(tmp,file=paste0(path2,"breast_proteome.txt"),quote=F,sep="\t",row.names=F)

# rna&cnv
rna_new=proteome;rna_new[]=NA
cnv_new=proteome;cnv_new[]=NA
ind=intersect(rownames(proteome),rownames(rna))
rna_new[ind,colnames(rna)]=rna[ind,colnames(rna)] # Here be carefull! Keep safe
ind=intersect(rownames(proteome),rownames(cnv))
cnv_new[ind,colnames(cnv)]=cnv[ind,colnames(cnv)]
dim(proteome);dim(rna_new);dim(cnv_new)
#[1] 10006    77
#[1] 10006    77
#[1] 10006    77

tmp=rna_new;tmp=cbind(rownames(tmp),tmp);colnames(tmp)[1]="Gene_ID"
write.table(tmp,file=paste0(path2,"trimmed_breast_rna.txt"),quote=F,sep="\t",row.names=F)
tmp=cnv_new;tmp=cbind(rownames(tmp),tmp);colnames(tmp)[1]="Gene_ID"
write.table(tmp,file=paste0(path2,"trimmed_breast_cnv.txt"),quote=F,sep="\t",row.names=F)
#######################################














