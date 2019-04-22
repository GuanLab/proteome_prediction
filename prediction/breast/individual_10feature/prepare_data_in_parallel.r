## name: prepare_data_in_parallel.r
## date: 10/10/2017

args <- commandArgs(trailingOnly = TRUE)
i=as.numeric(args[1])
num_parallel=as.numeric(args[2])

path0="../../../data/trimmed_set/";
path1="../../../data/normalization/raw/";
path2="../../../data/cv_set/";

feature=as.matrix(read.delim(paste0(path1,"raw_breast_rna.txt.avg_imputation.top10"),check.names=F,row.names=1))
train=as.matrix(read.delim(paste0(path2,"train_breast_proteome_",i,".txt"),check.names=F,row.names=1))
test=as.matrix(read.delim(paste0(path2,"test_breast_proteome_",i,".txt"),check.names=F,row.names=1))

feature_train=feature[,colnames(train)]
feature_test=feature[,colnames(test)]
# dim(train);dim(test);dim(feature_train);dim(feature_test)

write.table(feature_train,file="feature_train.dat",quote=F,sep="\t",na="NA",col.names=F,row.names=F)
write.table(feature_test,file="feature_test.dat",quote=F,sep="\t",na="NA",col.names=F,row.names=F)

system(paste0("mkdir o",i))
system(paste0("cp fast_prediction.py o",i))
system(paste0("mv feature_*.dat o",i))

size_subset=ceiling(dim(train)[1]/num_parallel)
for(j in 1:num_parallel){
	first=(j-1)*size_subset+1
	last=j*size_subset
	if(last>dim(train)[1]) {last=dim(train)[1]}
	train_subset=train[first:last,]

	system(paste0("cp -r o",i," o",i,j))
	write.table(train_subset,file=paste0("o",i,j,"/train.dat"),quote=F,sep="\t",na="NA",col.names=F,row.names=F)
}


