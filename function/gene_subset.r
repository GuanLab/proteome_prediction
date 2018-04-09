## name: gene_subset.r
## date: 10/09/2017

## Here I trim input matrices into a sub-matrices according to the rowname list
## The 1st argument is the appendix name, the 2nd is the rowname list, others are input matrices

args <- commandArgs(trailingOnly = TRUE)

name=args[1]
gene=readLines(args[2])

for(i in args[-c(1:2)]){
	input=as.matrix(read.delim(i,check.names=F,row.names=1))
	output=input[rownames(input) %in% gene,]
	
	output=cbind(rownames(output),output);colnames(output)[1]="Gene_ID"
	write.table(output,file=paste0(i,".",name),quote=F,sep="\t",row.names=F)
}


