## name: recenter_prediction.r
## date: 10/17/2017

## Here I recenter each entry/row by target avg&sd
## Rscript recenter_prediction.r ref_avg_sd.txt pred1.txt pred2.txt ..

args <- commandArgs(trailingOnly = TRUE)

ref=as.matrix(read.delim(args[1],check.names=F,row.names=1,header=F))

for(i in args[-1]){
	input=as.matrix(read.delim(i,check.names=F,row.names=1))
	d1=dim(input)[1];d2=dim(input)[2]

	mat_avg1=matrix(rep(ref[,1],times=d2),nrow=d1)
	mat_sd1=matrix(rep(ref[,2],times=d2),nrow=d1)
	mat_avg2=matrix(rep(apply(input,1,mean,na.rm=T),times=d2),nrow=d1)
	mat_sd2=matrix(rep(apply(input,1,sd,na.rm=T),times=d2),nrow=d1)

	output=(input-mat_avg2)/mat_sd2*mat_sd1+mat_avg1
	# UPGRADE what if mat_sd2=0?	

        output=cbind(rownames(output),output);colnames(output)[1]="Gene_ID"
        write.table(output,file=paste0(i,".recentered"),quote=F,sep="\t",row.names=F)
}




