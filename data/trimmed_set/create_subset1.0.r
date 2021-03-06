## name: create_subset1.0_list.r
## date: 11/05/2017

#path1="/state3/hyangl/CPTAC_final/sub3/data/raw/"
path1="../raw/"

##### 1. ova #################
proteome=as.matrix(read.delim(paste0(path1,"retrospective_ova_PNNL_proteome_sort_common_gene_7061.txt"),check.names=F,row.names=1))
rna=as.matrix(read.delim(paste0(path1,"retrospective_ova_rna_seq_sort_common_gene_15121.txt"),check.names=F,row.names=1))
cnv=as.matrix(read.delim(paste0(path1,"retrospective_ova_CNA_sort_common_gene_11859.txt"),check.names=F,row.names=1))

## rna
ind=rownames(rna) %in% rownames(proteome)
rna=rna[ind,]
#writeLines(rownames(rna),con="list_ova_rna_subset_proteome.txt") # 6021

num=apply(rna,1,function(x){
	sum(x==min(x,na.rm=T))
	})
ind=num==1
writeLines(rownames(rna)[ind],con="list_ova_rna_subset1.0.txt") # 5837

## cnv
ind=rownames(cnv) %in% rownames(proteome)
cnv=cnv[ind,]
#writeLines(rownames(cnv),con="list_ova_cnv_subset_proteome.txt") # 4835

num=apply(cnv,1,function(x){
        sum(x==min(x,na.rm=T))
        })
ind=num==1
writeLines(rownames(cnv)[ind],con="list_ova_cnv_subset1.0.txt") # 4834

##################################

##### 2. breast #################
proteome=as.matrix(read.delim(paste0(path1,"retrospective_breast_proteome_sort_common_gene_10005.txt"),check.names=F,row.names=1))
rna=as.matrix(read.delim(paste0(path1,"retrospective_breast_RNA_sort_common_gene_15107.txt"),check.names=F,row.names=1))
cnv=as.matrix(read.delim(paste0(path1,"retrospective_breast_CNA_sort_common_gene_16884.txt"),check.names=F,row.names=1))

## rna
ind=rownames(rna) %in% rownames(proteome)
rna=rna[ind,] # 9032
sum(is.na(rna)) #[1] 1400 breast rna has NA!

cutoff=dim(rna)[2]*1.0
ind=apply(rna,1,function(x){sum(!is.na(x))>=cutoff})

#num=apply(rna,1,function(x){ # breast rna didn't fill NA with min!!
#        sum(x==min(x,na.rm=T))
#        })
#ind=num==1
writeLines(rownames(rna)[ind],con="list_breast_rna_subset1.0.txt") # 8738

## cnv
ind=rownames(cnv) %in% rownames(proteome)
cnv=cnv[ind,] # 8853
sum(is.na(cnv)) #[1] no NA

cutoff=dim(cnv)[2]*1.0
ind=apply(cnv,1,function(x){sum(!is.na(x))>=cutoff})
cnv=cnv[ind,] # 8853

num=apply(cnv,1,function(x){ # breast cnv didn't fill NA with min!!
        sum(x==min(x,na.rm=T))
        })
ind=num==1
writeLines(rownames(cnv)[ind],con="list_breast_cnv_subset1.0.txt") # 8839










