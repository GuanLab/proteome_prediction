## name: create_gene_list.r
## date: 02/05/2018

mat=as.matrix(read.delim("fi.txt",row.names=1,check.names=F))

fi=apply(mat,2,sum)
fi_freq=apply(mat,2,function(x){sum(x!=0)})

writeLines(paste0(names(fi),"\t",fi),con="fi_score_breast.txt")


