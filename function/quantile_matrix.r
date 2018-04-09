## name: quantile_matrix.r
## date: 09/25/2017

## This function quantile map query matrix to reference matrix, ignoring NA values

args <- commandArgs(trailingOnly = TRUE)

quantile_mat=function(mat_query,mat_ref){
	ref=as.vector(data.matrix(mat_ref))
	ref=ref[!is.na(ref)]
	query=as.vector(data.matrix(mat_query))
	ind=!is.na(query)

	quantiled_ref=quantile(ref,probs=seq(0,1,1/(length(query[ind])-1))) # quantile
	names(quantiled_ref)=NULL

	query_new=query
	query_new[ind]=quantiled_ref[rank(query[ind])]
	mat_out=matrix(query_new,nrow=dim(mat_query)[1])
	rownames(mat_out)=rownames(mat_query)
	colnames(mat_out)=colnames(mat_query)
	return(mat_out)
}

ref=as.matrix(read.delim(args[1],stringsAsFactors=F,check.names=F,row.names=1))

for (j in args[-1]){
	query=as.matrix(read.delim(j,stringsAsFactors=F,check.names=F,row.names=1))
	output=quantile_mat(query,ref)
	output=cbind(rownames(output),output);colnames(output)[1]="Gene_ID"
	write.table(output,file=paste0(j,".quantile_mat"),quote=F,sep="\t",row.names=F)
}



