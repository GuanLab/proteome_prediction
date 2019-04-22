## name: cv_partition.r
## date: 12/07/2017

path1="../trimmed_set/"
path2="./"

# set.seed(449) I didn't set seed but experiments are done so I just keep the partition...

fold=5

# breast
proteome=as.matrix(read.delim(paste0(path1,"breast_proteome.txt"),row.names=1,check.names=F))
num_sample=dim(proteome)[2]
size_subset=round(num_sample/fold,0)
ind=sample(num_sample,num_sample)

for(j in 1:fold){
        first=(j-1)*size_subset+1
        last=j*size_subset
        if(j==fold & last!=num_sample) {last=num_sample}
        train_subset=proteome[,-sort(ind[first:last])]
        test_subset=proteome[,sort(ind[first:last])]
	tmp=train_subset;tmp=cbind(rownames(tmp),tmp);colnames(tmp)[1]="Gene_ID"
        write.table(tmp,file=paste0(path2,"train_breast_proteome_",j,".txt"),quote=F,sep="\t",na="NA",col.names=T,row.names=F)
	tmp=test_subset;tmp=cbind(rownames(tmp),tmp);colnames(tmp)[1]="Gene_ID"
        write.table(tmp,file=paste0(path2,"test_breast_proteome_",j,".txt"),quote=F,sep="\t",na="NA",col.names=T,row.names=F)
}

# ova
proteome=as.matrix(read.delim(paste0(path1,"ova_proteome.txt"),row.names=1,check.names=F))
num_sample=dim(proteome)[2]
size_subset=round(num_sample/fold,0)
ind=sample(num_sample,num_sample)

for(j in 1:fold){
        first=(j-1)*size_subset+1
        last=j*size_subset
        if(j==fold & last!=num_sample) {last=num_sample}
        train_subset=proteome[,-sort(ind[first:last])]
        test_subset=proteome[,sort(ind[first:last])]
	tmp=train_subset;tmp=cbind(rownames(tmp),tmp);colnames(tmp)[1]="Gene_ID"
        write.table(tmp,file=paste0(path2,"train_ova_proteome_",j,".txt"),quote=F,sep="\t",na="NA",col.names=T,row.names=F)
	tmp=test_subset;tmp=cbind(rownames(tmp),tmp);colnames(tmp)[1]="Gene_ID"
        write.table(tmp,file=paste0(path2,"test_ova_proteome_",j,".txt"),quote=F,sep="\t",na="NA",col.names=T,row.names=F)
}


