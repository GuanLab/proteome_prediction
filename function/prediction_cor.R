## name: prediction_cor.R
## date: 10/01/2017

## Here for each row/entry, I calculate the correlation between two matrices and report cor quantile
## This is used to compare two matrix (e.g. before and after quantile)
## Since the number of rows/genes could be 10,000+, I only report quantile

cor_quantile=function(mat1,mat2,probs=seq(0,1,0.01)){
	vec_cor=NULL
	for(i in 1:dim(mat1)[1]){
		vec_cor=c(vec_cor,cor(mat1[i,],mat2[i,],use="pairwise.complete.obs"))
	}
	output=quantile(vec_cor,probs=probs,na.rm=T)
	return(output)
}

options(warn=-1)

args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=2){print("Wrong arguments!\nRscript prediction_cor.R prediciton1.txt prediciton2.txt\n")}

mat1=as.matrix(read.delim(args[1],check.names=F,row.names=1))
mat2=as.matrix(read.delim(args[2],check.names=F,row.names=1))

print(cor_quantile(mat1,mat2,probs=seq(0,1,0.1)))






