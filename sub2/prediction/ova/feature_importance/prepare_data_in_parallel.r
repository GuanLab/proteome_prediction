## name: prepare_data_in_parallel.r
## date: 10/10/2017

## this is adapted from CPTAC_round2/sub2/prediction/m183
## I didn't rerun it

args <- commandArgs(trailingOnly = TRUE)
i=as.numeric(args[1])
num_parallel=as.numeric(args[2])

path0="../../../data/trimmed_set/";
path1="../../../data/normalization/raw/";
path2="../../../data/cv_set/";

feature=as.matrix(read.delim(paste0(path1,"raw_ova_rna.txt.avg_imputation.subset1.0"),check.names=F,row.names=1))
pnnl=as.matrix(read.delim(paste0(path0,"ova_proteome_pnnl.txt"),check.names=F,row.names=1))
jhu=as.matrix(read.delim(paste0(path0,"ova_proteome_jhu.txt"),check.names=F,row.names=1))

##### fill pnnl NA with jhu_overlap ###
jhu_overlap=as.matrix(read.delim(paste0(path0,"ova_proteome_jhu_overlap.txt"),check.names=F,row.names=1))
name_overlap=colnames(jhu_overlap)[colnames(jhu_overlap) %in% colnames(pnnl)]
pnnl[,name_overlap][is.na(pnnl[,name_overlap])]=jhu_overlap[,name_overlap][is.na(pnnl[,name_overlap])]
#######################################

train=cbind(pnnl,jhu)
test=as.matrix(read.delim(paste0(path2,"test_ova_proteome_",i,".txt"),check.names=F,row.names=1))

feature_train=feature[,colnames(train)]
feature_test=feature[,colnames(test)]
# dim(train);dim(test);dim(feature_train);dim(feature_test)

write.table(feature_train,file="feature_train.dat",quote=F,sep="\t",na="NA",col.names=F,row.names=F)
write.table(feature_test,file="feature_test.dat",quote=F,sep="\t",na="NA",col.names=F,row.names=F)

system("mkdir o")
#system("cp fast_prediction.py o")
system("cp rf_save.py o")
system("mv feature_*.dat o")

size_subset=ceiling(dim(train)[1]/num_parallel)
for(j in 1:num_parallel){
	first=(j-1)*size_subset+1
	last=j*size_subset
	if(last>dim(train)[1]) {last=dim(train)[1]}
	train_subset=train[first:last,]

	system(paste0("cp -r o o",j))
	write.table(train_subset,file=paste0("o",j,"/train.dat"),quote=F,sep="\t",na="NA",col.names=F,row.names=F)
	writeLines(rownames(train_subset),con=paste0("o",j,"/genename.txt"))
}


