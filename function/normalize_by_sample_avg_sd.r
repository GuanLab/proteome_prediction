## name: normalized_by_sample_avg_sd.r
## date: 10/09/2017

## Here I normalize the data by sample (subtracting) mean & (divided by) sd

args <- commandArgs(trailingOnly = TRUE)

for(i in args){
	input=as.matrix(read.delim(i,check.names=F,row.names=1))

	mat_mean=matrix(rep(apply(input,2,mean,na.rm=T),each=dim(input)[1]),nrow=dim(input)[1])
	mat_sd=matrix(rep(apply(input,2,sd,na.rm=T),each=dim(input)[1]),nrow=dim(input)[1])

	# in case abnormal input
	mat_mean[is.na(mat_mean)]=0
	mat_sd[is.na(mat_sd)]=1

	output=(input-mat_mean)/mat_sd

	output=cbind(rownames(output),output);colnames(output)[1]="Gene_ID"
	write.table(output,file=paste0(i,".scaled"),quote=F,sep="\t",row.names=F)	
}






