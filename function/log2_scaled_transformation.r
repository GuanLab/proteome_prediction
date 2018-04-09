## name: log2_scaled_transformation.r
## date: 09/25/2017

## Here I calculate log2 transformation of scaled input matrix log2(x/mean(x))
## This function can NOT be used generally - only for leaderboard RNA
## lb RNA are all positive numbers and NA are labeled as 0!!

args <- commandArgs(trailingOnly = TRUE)

options(warn=-1)

for(i in args){
	input=as.matrix(read.delim(i,check.names=F,row.names=1))

	output=input
	input[is.na(input)]=0 # what if there exist NA and 0? This version treat 0 as NA
	for(j in 1:dim(input)[1]){
		ind=input[j,]!=0
		output[j,ind]=log(input[j,ind]/mean(input[j,ind]),base=2) #UPGRADE: assign min first
		output[j,!ind]=min(output[j,ind]) # assign minimum values to NA positions
	}
	output[is.na(output)]=0 # set abnormal values to 0
	output[is.infinite(output)]=0

	tmp=output;tmp=cbind(rownames(tmp),tmp);colnames(tmp)[1]="Gene_ID"
	write.table(tmp,file=paste0(i,".log2"),quote=F,sep="\t",row.names=F)
}


