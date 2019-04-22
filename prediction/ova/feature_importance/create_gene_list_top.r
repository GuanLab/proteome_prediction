## name: create_gene_list.r
## date: 02/05/2018

mat=as.matrix(read.delim("fi.txt",row.names=1,check.names=F))

fi=apply(mat,2,sum)
fi_freq=apply(mat,2,function(x){sum(x!=0)})

fi1000=sort(fi,decreasing=T)[1:1000]

writeLines(names(fi1000),con="top1000.txt")

writeLines(paste0(names(fi1000),"\t",fi1000),con="top1000_score.txt")


