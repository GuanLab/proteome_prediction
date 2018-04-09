## name: stack_model.r
## date: 08/15/2017

path0="./"
path1="../rna/"
path2="../individual/" 

for (i in 1:5){

    w1=1
    w2=3
    w=w1+w2

    x1=as.matrix(read.delim(paste0(path1,"prediction_breast_proteome_",i,".txt"),check.names=F,row.names=1))
    x2=as.matrix(read.delim(paste0(path2,"prediction_breast_proteome_",i,".txt"),check.names=F,row.names=1))

    #x2[is.na(x2)]=x1[is.na(x2)]
    y=w1/w*x1 + w2/w*x2
#y=x2

    tmp=cbind(rownames(y),y);colnames(tmp)[1]="proteinID"
    write.table(tmp,file=paste0(path0, "prediction_breast_proteome_",i,".txt"),quote=F,sep="\t",row.names=F)
}
system("./sub2_score.sh")


