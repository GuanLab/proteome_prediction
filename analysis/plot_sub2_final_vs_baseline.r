
## name: plot_sub2_final_vs_baseline.r
## date: 02/22/2018

library(ggplot2)
library(reshape2)
library(plyr)
library(scales)

source("~/function/my_palette.r")
source("~/function/multiplot.R")

fold=5

## breast
path1="../sub2/prediction/breast/rna/"
path2="../sub2/prediction/breast/final/"
ymax=2

tmp=as.matrix(read.delim(paste0(path1,"cor_nrmse_b",1,".txt"),header=F,row.names=1))
mat_final=mat_baseline=matrix(NA,nrow=dim(tmp),ncol=fold)
rownames(mat_baseline)=rownames(mat_final)=rownames(tmp)
for(i in 1:fold){
    mat_baseline[,i]=as.matrix(read.delim(paste0(path1,"cor_nrmse_b",i,".txt"),header=F,row.names=1))[,1]
    mat_final[,i]=as.matrix(read.delim(paste0(path2,"cor_nrmse_b",i,".txt"),header=F,row.names=1))[,1]
}

bb=apply(mat_baseline,1,mean,na.rm=T)
bf=apply(mat_final,1,mean,na.rm=T)
ind=!is.na(bb)
bb=bb[ind]
bf=bf[ind]

tbl=rbind(data.frame(model=rep("baseline",length(bb)),correlation=bb),
        data.frame(model=rep("final",length(bf)),correlation=bf))
df=data.frame(model=c("baseline","final"), x1=c(mean(bb),mean(bf)), 
    x2=c(mean(bb),mean(bf)), y1=0, y2=ymax)
tmp_col=my_palette[c("bluegrey","red")]
names(tmp_col)=c("baseline","final")
pb=ggplot(tbl, aes(correlation, fill = model, colour = model)) +
    geom_density(alpha = 0.3) +
    scale_fill_manual(values=tmp_col)  +
    scale_colour_manual(values=tmp_col) +
    ylim(0,ymax) +
    geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2), linetype=2) +
    labs(title="breast") +
    theme(plot.title = element_text(hjust = 0.5))

## ova
path1="../sub2/prediction/ova/rna/"
path2="../sub2/prediction/ova/final/"
ymax=2

tmp=as.matrix(read.delim(paste0(path1,"cor_nrmse_o",1,".txt"),header=F,row.names=1))
mat_final=mat_baseline=matrix(NA,nrow=dim(tmp),ncol=fold)
rownames(mat_baseline)=rownames(mat_final)=rownames(tmp)
for(i in 1:fold){
    mat_baseline[,i]=as.matrix(read.delim(paste0(path1,"cor_nrmse_o",i,".txt"),header=F,row.names=1))[,1]
    mat_final[,i]=as.matrix(read.delim(paste0(path2,"cor_nrmse_o",i,".txt"),header=F,row.names=1))[,1]
}

ob=apply(mat_baseline,1,mean,na.rm=T)
of=apply(mat_final,1,mean,na.rm=T)
ind=!is.na(ob)
ob=ob[ind]
of=of[ind]

tbl=rbind(data.frame(model=rep("baseline",length(ob)),correlation=ob),
        data.frame(model=rep("final",length(of)),correlation=of))
df=data.frame(model=c("baseline","final"), x1=c(mean(ob),mean(of)), 
    x2=c(mean(ob),mean(of)), y1=0, y2=ymax)
tmp_col=my_palette[c("bluegrey","red")]
names(tmp_col)=c("baseline","final")
po=ggplot(tbl, aes(correlation, fill = model, colour = model)) +
    geom_density(alpha = 0.3) +
    scale_fill_manual(values=tmp_col)  +
    scale_colour_manual(values=tmp_col) +
    ylim(0,ymax) +
    geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2), linetype=2) +
    labs(title="ovary") +
    theme(plot.title = element_text(hjust = 0.5))

list_p=list()
list_p[[1]]=pb
list_p[[2]]=po

pdf(file="../figure/sub2_final_vs_baseline.pdf",width=8,height=6)
mat_layout=matrix(1:2,nrow=2)
multiplot(plotlist=list_p,layout = mat_layout)
dev.off()

