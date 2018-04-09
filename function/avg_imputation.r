## name: avg_imputation.r
## date: 07/13/2017

## Here I imputate input missing data
## 1. For partially missing rows, simply fill with the average of non-missing values
## 2. For completely missing rows, fill with the sample average


args <- commandArgs(trailingOnly = TRUE)

for (j in args){
	input=as.matrix(read.delim(j,stringsAsFactors=F,check.names=F,row.names=1))

	output=input
	ind = which(!apply(output,1,function(x){all(is.na(x))}))
	for(i in ind){
	    tmp=is.na(output[i,])
	    output[i,tmp]=mean(output[i,!tmp],na.rm=T)
	}
	sample_avg=apply(output,2,mean,na.rm=T)
	ind = apply(output,1,function(x){all(is.na(x))}) # completely missing rows
	output[ind,]=rep(sample_avg,each=sum(ind))

	output=cbind(rownames(output),output);colnames(output)[1]="Gene_ID"
	write.table(output,file=paste0(j,".avg_imputation"),quote=F,sep="\t",row.names=F)
}


