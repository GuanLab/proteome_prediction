## name: sub2_hongyang.r
## date: 07/07/2017

args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=2){print("Wrong arguments!\nRscript sub2_hongyang.r data_true.txt data_pred.txt\n")}


mat_true=as.matrix(read.delim(args[1],row.names=1,check.names=F))
mat_pred=as.matrix(read.delim(args[2],row.names=1,check.names=F))

d1=dim(mat_true)[1]
vec_cor=rep(NA,d1)
vec_nrmse=rep(NA,d1)
#cutoff=dim(mat_true)[2]*0.7
cutoff=5
for(i in 1:d1){
	ind=!is.na(mat_true[i,])
	y=mat_true[i,ind]
        x=mat_pred[i,ind]
	## 1.correlation
	if(sum(ind)>=cutoff){
	        if(sd(x)==0 | sd(y)==0 | sum(ind)<2){
        	        vec_cor[i]=0 
	        }else{
        	        vec_cor[i]=cor(x,y)
        	}
	}
	## 2.nrmse
	if(sum(ind)>4){ # for proteins with at least 4 observations
	        y_diff=range(y)[2]-range(y)[1]
	        vec_nrmse[i]=(sum((y-x)^2)/length(y))^0.5/y_diff
	}
}

avg_nrmse=mean(vec_nrmse,na.rm=T)
avg_cor=mean(vec_cor,na.rm=T)

#print(c(avg_cor,avg_nrmse))
#print(paste0(avg_cor,"_",avg_nrmse))
print(avg_cor)

mat=cbind(vec_cor,vec_nrmse)
rownames(mat)=rownames(mat_true)
write.table(mat, file="cor_nrmse.txt",quote=F,sep="\t",col.names=F)


