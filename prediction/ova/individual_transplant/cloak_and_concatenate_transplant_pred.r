## name: cloak_and_concatenate_transplant_pred.r
## date: 10/10/2017

args <- commandArgs(trailingOnly = TRUE)
i=as.numeric(args[1])
num_parallel=as.numeric(args[2])

path2="../../../data/cv_set/";

test=as.matrix(read.delim(paste0(path2,"test_ova_proteome_",i,".txt"),check.names=F,row.names=1))
common_target=readLines("common_target.txt")

pred_transplant=NULL
for(j in 1:num_parallel){
	mat=as.matrix(read.delim(paste0("o",i,j,"/prediction.dat"),header=F))
	pred_transplant=rbind(pred_transplant,mat)
}

pred=test
pred[]=NA
pred[common_target,]=pred_transplant

tmp=pred;tmp=cbind(rownames(tmp),tmp);colnames(tmp)[1]="proteinID"
write.table(tmp,file=paste0("prediction_bo_proteome_",i,".txt"),quote=F,sep="\t",row.names=F)

