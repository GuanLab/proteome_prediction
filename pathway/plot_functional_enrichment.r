## name: plot_functional_enrichment.r
## date: 03/20/2018

library(ggplot2)
library(reshape2)
library(plyr)
library(scales)

source("~/function/my_palette.r")
source("~/function/multiplot.R")

path1="../sub2/prediction/breast/"
path2="../sub2/prediction/ova/"

fold=5
bg=bi=bf=NULL
og=oi=of=NULL
for(i in 1:fold){
    tmp=read.delim(paste0(path1,"rna/cor_nrmse_b",i,".txt"),header=F,row.names=1)
    bg=c(bg,tmp[,1])
    tmp=read.delim(paste0(path1,"global_individual/cor_nrmse_b",i,".txt"),header=F,row.names=1)
    bi=c(bi,tmp[,1])
    tmp=read.delim(paste0(path1,"final/cor_nrmse_b",i,".txt"),header=F,row.names=1)
    bf=c(bf,tmp[,1])

    tmp=read.delim(paste0(path2,"rna/cor_nrmse_o",i,".txt"),header=F,row.names=1)
    og=c(og,tmp[,1])
    tmp=read.delim(paste0(path2,"global_individual/cor_nrmse_o",i,".txt"),header=F,row.names=1)
    oi=c(oi,tmp[,1])
    tmp=read.delim(paste0(path2,"final/cor_nrmse_o",i,".txt"),header=F,row.names=1)
    of=c(of,tmp[,1])
}

bg=apply(matrix(bg,ncol=fold),1,mean,na.rm=T)
bi=apply(matrix(bi,ncol=fold),1,mean,na.rm=T)
bf=apply(matrix(bf,ncol=fold),1,mean,na.rm=T)
dbf=bf-bg

og=apply(matrix(og,ncol=fold),1,mean,na.rm=T)
oi=apply(matrix(oi,ncol=fold),1,mean,na.rm=T)
of=apply(matrix(of,ncol=fold),1,mean,na.rm=T)
dof=of-og

tmp=read.delim(paste0(path1,"rna/cor_nrmse_b1.txt"),header=F,row.names=1)
bname=rownames(tmp)
tmp=read.delim(paste0(path2,"rna/cor_nrmse_o1.txt"),header=F,row.names=1)
oname=rownames(tmp)

names(bg)=names(bi)=names(bf)=names(dbf)=bname
names(og)=names(oi)=names(of)=names(dof)=oname

mat=as.matrix(read.delim(paste0(path1,"feature_importance/fi.txt"),row.names=1,check.names=F))
bfi=apply(mat,2,sum);names(bfi)=colnames(mat)
mat=as.matrix(read.delim(paste0(path2,"feature_importance/fi.txt"),row.names=1,check.names=F))
ofi=apply(mat,2,sum);names(ofi)=colnames(mat)

# kegg
cutoff_fdr=0.01
map_kegg=readLines("kegg_category.txt")
cate=c("1. Metabolism","2. Genetic Information Processing","3. Environmental Information Processing",
    "4. Cellular Processes","5. Organismal Systems","6. Human Diseases","7. Drug Development")
ind_cate=NULL
for(i in 1:length(cate)){
    ind_cate=c(ind_cate,which(map_kegg==cate[i]))
}
ind_cate=c(ind_cate,length(map_kegg))

tmp_col=my_palette[c("bluegrey","red","red","blue","teal","yellow","cyan","purple","brown")]
names(tmp_col)=c("baseline","final",sub("[0-9]\\. ","",cate))

#scaleFUN <- function(x) sprintf("%.1f", x)

## 1. final
ymax=2

#tbl=rbind(data.frame(model=rep("baseline",length(bg)),correlation=bg),
#        data.frame(model=rep("final",length(bf)),correlation=bf))
tbl=data.frame(model=rep("final",length(bf)),correlation=bf)
df=data.frame(model=c("baseline","final"), x1=c(mean(bg,na.rm=T),mean(bf,na.rm=T)),
    x2=c(mean(bg,na.rm=T),mean(bf,na.rm=T)), y1=0, y2=ymax)
b1=ggplot(tbl, aes(correlation, fill = model, colour = model)) +
    geom_density(alpha = 0.3) +
    scale_fill_manual(values=tmp_col)  +
    scale_colour_manual(values=tmp_col) +
    ylim(0,ymax) +
    xlim(-1,1) +
    #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2), linetype=2) +
    labs(title="breast") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),legend.position="none")
    #annotate("text", x = c(mean(bg,na.rm=T),mean(bf,na.rm=T)), y = c(ymax,ymax), 
    #   label=format(c(mean(bg,na.rm=T),mean(bf,na.rm=T)),digits=3))

tbl=data.frame(model=rep("final",length(of)),correlation=of)
df=data.frame(model=c("baseline","final"), x1=c(mean(og,na.rm=T),mean(of,na.rm=T)),
    x2=c(mean(og,na.rm=T),mean(of,na.rm=T)), y1=0, y2=ymax)
o1=ggplot(tbl, aes(correlation, fill = model, colour = model)) +
    geom_density(alpha = 0.3) +
    scale_fill_manual(values=tmp_col)  +
    scale_colour_manual(values=tmp_col) +
    ylim(0,ymax) +
    xlim(-1,1) +
    #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2), linetype=2) +
    labs(title="ovary") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),legend.position="none")
    

tbl=NULL
for(i in 1:4){
    tbl=rbind(tbl,read.delim(paste0("david/bf",i,".txt"),sep="\t",stringsAsFactors=F,check.names=F))
}
ind=grep("KEGG_PATHWAY",tbl[,"Category"])
ind=ind[tbl[ind,"FDR"]<cutoff_fdr]
id=sub("hsa","",sub("\\:.*","",tbl[ind,"Term"]))
kegg=sub(".*\\:","",tbl[ind,"Term"])
for(i in 1:length(id)){ # find major categories
    tmp_num=which(map_kegg==id[i])
    for(j in 1:length(cate)){
        if(tmp_num>ind_cate[j] && tmp_num<=ind_cate[j+1]){
            names(kegg)[i]=sub("[0-9]\\. ","",cate[j])
        }
    }
}
tbl1=tbl;ind1=ind;bkegg=kegg

tbl=NULL
for(i in 1:4){
    tbl=rbind(tbl,read.delim(paste0("david/of",i,".txt"),sep="\t",stringsAsFactors=F,check.names=F))
}
ind=grep("KEGG_PATHWAY",tbl[,"Category"])
ind=ind[tbl[ind,"FDR"]<cutoff_fdr]
id=sub("hsa","",sub("\\:.*","",tbl[ind,"Term"]))
kegg=sub(".*\\:","",tbl[ind,"Term"])
for(i in 1:length(id)){ # find major categories
    tmp_num=which(map_kegg==id[i])
    for(j in 1:length(cate)){
        if(tmp_num>ind_cate[j] && tmp_num<=ind_cate[j+1]){
            names(kegg)[i]=sub("[0-9]\\. ","",cate[j])
        }
    }
}
tbl2=tbl;ind2=ind;okegg=kegg

df=NULL
ymax=4.9
pos=0
for(i in 1:length(ind1)){
    gene=unlist(strsplit(tbl1[ind1[i],"Genes"],split=", "))
    df=rbind(df,data.frame(kegg=names(bkegg)[i],x1=bf[gene],x2=bf[gene],y1=pos,y2=pos+0.1))
    pos=pos+0.15
}
df_kegg=data.frame(kegg=bkegg,x=rep(0,length(bkegg)),y=seq(0,pos-0.1,0.15))
bk1=ggplot(tbl, aes(correlation, fill = kegg, colour = kegg)) +
    #geom_density(alpha = 0) +
    scale_fill_manual(values=tmp_col)  +
    scale_colour_manual(values=tmp_col) +
    ylim(0,ymax) +
    xlim(-1,1) +
    geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2), linetype=1) +
    #scale_y_continuous(labels=scaleFUN) +
    labs(title="breast") +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5),legend.position="none") +
    annotate("text", x = rep(-1,length(bkegg)), y = seq(0,pos-0.1,0.15)+0.05, label=bkegg, hjust=0,size=3)


kegg_legend=ggplot(tbl, aes(correlation, fill = kegg, colour = kegg)) +
    #geom_density(alpha = 0) +
    scale_fill_manual(values=tmp_col)  +
    scale_colour_manual(values=tmp_col) +
    ylim(0,ymax) +
    xlim(-1,1) +
    geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2), linetype=1) +
    #scale_y_continuous(labels=scaleFUN) +
    labs(title="breast") +
    theme(plot.title = element_text(hjust = 0.5)) +
    annotate("text", x = rep(-1,length(bkegg)), y = seq(0,pos-0.1,0.15)+0.05, label=bkegg, hjust=0)

df=NULL
pos=0
for(i in 1:length(ind2)){
    gene=unlist(strsplit(tbl2[ind2[i],"Genes"],split=", "))
    df=rbind(df,data.frame(kegg=names(okegg)[i],x1=of[gene],x2=of[gene],y1=pos,y2=pos+0.1))
    pos=pos+0.15
}
df_kegg=data.frame(kegg=okegg,x=rep(0,length(okegg)),y=seq(0,pos-0.1,0.15))
ok1=ggplot(tbl, aes(correlation, fill = kegg, colour = kegg)) +
    #geom_density(alpha = 0) +
    scale_fill_manual(values=tmp_col)  +
    scale_colour_manual(values=tmp_col) +
    ylim(0,ymax) +
    xlim(-1,1) +
    geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2), linetype=1) +
    #scale_y_continuous(labels=scaleFUN) +
    labs(title="ovary") +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5),legend.position="none") +
    annotate("text", x = rep(-1,length(okegg)), y = seq(0,pos-0.1,0.15)+0.05, label=okegg, hjust=0,size=3)


## 2. final - global
ymax=0.9

dbf=dbf[!is.na(dbf)&dbf>0]
dbf=log(dbf,base=10)
tbl=rbind(data.frame(model=rep("final",length(dbf)),correlation=dbf))
df=data.frame(model=c("final"), x1=c(mean(dbf,na.rm=T)),
    x2=c(mean(dbf,na.rm=T)), y1=0, y2=ymax)
b2=ggplot(tbl, aes(correlation, fill = model, colour = model)) +
    geom_density(alpha = 0.3) +
    scale_fill_manual(values=tmp_col)  +
    scale_colour_manual(values=tmp_col) +
    ylim(0,ymax) +
    xlim(-6,0.2) +
    #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2), linetype=2) +
    labs(title="breast",x="increased correlation(log10)") +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5),legend.position="none") 
    #annotate("text", x = c(mean(dbf,na.rm=T)), y = c(ymax), 
    #   label=format(c(mean(dbf,na.rm=T)),digits=3))

dof=dof[!is.na(dof)&dof>0]
dof=log(dof,base=10)
tbl=rbind(data.frame(model=rep("final",length(dof)),correlation=dof))
df=data.frame(model=c("final"), x1=c(mean(dof,na.rm=T)),
    x2=c(mean(dof,na.rm=T)), y1=0, y2=ymax)
o2=ggplot(tbl, aes(correlation, fill = model, colour = model)) +
    geom_density(alpha = 0.3) +
    scale_fill_manual(values=tmp_col)  +
    scale_colour_manual(values=tmp_col) +
    ylim(0,ymax) +
    xlim(-6,0.2) +
    labs(title="ovary",x="increased correlation(log10)") +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5),legend.position="none") 

cutoff_fdr=1
tbl=NULL
for(i in 1:4){
    tbl=rbind(tbl,read.delim(paste0("david/dbf",i,".txt"),sep="\t",stringsAsFactors=F,check.names=F))
}
ind=grep("KEGG_PATHWAY",tbl[,"Category"])
ind=ind[tbl[ind,"FDR"]<cutoff_fdr]
id=sub("hsa","",sub("\\:.*","",tbl[ind,"Term"]))
kegg=sub(".*\\:","",tbl[ind,"Term"])
for(i in 1:length(id)){ # find major categories
    tmp_num=which(map_kegg==id[i])
    for(j in 1:length(cate)){
        if(tmp_num>ind_cate[j] && tmp_num<=ind_cate[j+1]){
            names(kegg)[i]=sub("[0-9]\\. ","",cate[j])
        }
    }
}
tbl1=tbl;ind1=ind;bkegg=kegg

cutoff_fdr=0.01
tbl=NULL
for(i in 1:4){
    tbl=rbind(tbl,read.delim(paste0("david/dof",i,".txt"),sep="\t",stringsAsFactors=F,check.names=F))
}
ind=grep("KEGG_PATHWAY",tbl[,"Category"])
ind=ind[tbl[ind,"FDR"]<cutoff_fdr]
id=sub("hsa","",sub("\\:.*","",tbl[ind,"Term"]))
kegg=sub(".*\\:","",tbl[ind,"Term"])
for(i in 1:length(id)){ # find major categories
    tmp_num=which(map_kegg==id[i])
    for(j in 1:length(cate)){
        if(tmp_num>ind_cate[j] && tmp_num<=ind_cate[j+1]){
            names(kegg)[i]=sub("[0-9]\\. ","",cate[j])
        }
    }
}
tbl2=tbl;ind2=ind;okegg=kegg

ymax=2.25
df=NULL
pos=0
for(i in 1:length(ind1)){
    gene=unlist(strsplit(tbl1[ind1[i],"Genes"],split=", "))
    df=rbind(df,data.frame(kegg=names(bkegg)[i],x1=dbf[gene],x2=dbf[gene],y1=pos,y2=pos+0.1))
    pos=pos+0.15
}
df_kegg=data.frame(kegg=bkegg,x=rep(0,length(bkegg)),y=seq(0,pos-0.1,0.15))
bk2=ggplot(tbl, aes(correlation, fill = kegg, colour = kegg)) +
    #geom_density(alpha = 0) +
    scale_fill_manual(values=tmp_col)  +
    scale_colour_manual(values=tmp_col) +
    ylim(0,ymax) +
    xlim(-6,0.2) +
    geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2), linetype=1) +
    labs(title="breast",x="increased correlation(log10)") +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5),legend.position="none") +
    annotate("text", x = rep(-6,length(bkegg)), y = seq(0,pos-0.1,0.15)+0.05, label=bkegg, hjust=0,size=3)

df=NULL
pos=0
for(i in 1:length(ind2)){
    gene=unlist(strsplit(tbl2[ind2[i],"Genes"],split=", "))
    df=rbind(df,data.frame(kegg=names(okegg)[i],x1=dof[gene],x2=dof[gene],y1=pos,y2=pos+0.1))
    pos=pos+0.15
}
df_kegg=data.frame(kegg=okegg,x=rep(0,length(okegg)),y=seq(0,pos-0.1,0.15))
ok2=ggplot(tbl, aes(correlation, fill = kegg, colour = kegg)) +
    #geom_density(alpha = 0) +
    scale_fill_manual(values=tmp_col)  +
    scale_colour_manual(values=tmp_col) +
    ylim(0,ymax) +
    xlim(-6,0.2) +
    geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2), linetype=1) +
    labs(title="ovary",x="increased correlation(log10)") +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5),legend.position="none") +
    annotate("text", x = rep(-6,length(okegg)), y = seq(0,pos-0.1,0.15)+0.05, label=okegg, hjust=0,size=3)

## feature importance
ymax=4.5

bfi=bfi[bfi>1]
bfi=log(bfi,base=10)
tbl=rbind(data.frame(model=rep("final",length(bfi)),correlation=bfi))
df=data.frame(model=c("final"), x1=c(mean(bfi,na.rm=T)),
    x2=c(mean(bfi,na.rm=T)), y1=0, y2=ymax)
b3=ggplot(tbl, aes(correlation, fill = model, colour = model)) +
    geom_density(alpha = 0.3) +
    scale_fill_manual(values=tmp_col)  +
    scale_colour_manual(values=tmp_col) +
    ylim(0,ymax) +
    xlim(0,0.8) +
    labs(title="breast",x="feature importance") +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5),legend.position="none") 

ofi=ofi[ofi>1]
ofi=log(ofi,base=10)
tbl=rbind(data.frame(model=rep("final",length(ofi)),correlation=ofi))
df=data.frame(model=c("final"), x1=c(mean(ofi,na.rm=T)),
    x2=c(mean(ofi,na.rm=T)), y1=0, y2=ymax)
o3=ggplot(tbl, aes(correlation, fill = model, colour = model)) +
    geom_density(alpha = 0.3) +
    scale_fill_manual(values=tmp_col)  +
    scale_colour_manual(values=tmp_col) +
    ylim(0,ymax) +
    xlim(0,0.8) +
    labs(title="ovary",x="feature importance") +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5),legend.position="none") 

cutoff_fdr=1
tbl=NULL
for(i in 1:4){
    tbl=rbind(tbl,read.delim(paste0("david/bfi",i,".txt"),sep="\t",stringsAsFactors=F,check.names=F))
}
ind=grep("KEGG_PATHWAY",tbl[,"Category"])
ind=ind[tbl[ind,"FDR"]<cutoff_fdr]
id=sub("hsa","",sub("\\:.*","",tbl[ind,"Term"]))
kegg=sub(".*\\:","",tbl[ind,"Term"])
for(i in 1:length(id)){ # find major categories
    tmp_num=which(map_kegg==id[i])
    for(j in 1:length(cate)){
        if(tmp_num>ind_cate[j] && tmp_num<=ind_cate[j+1]){
            names(kegg)[i]=sub("[0-9]\\. ","",cate[j])
        }
    }
}
tbl1=tbl;ind1=ind;bkegg=kegg

cutoff_fdr=1
tbl=NULL
for(i in 1:4){
    tbl=rbind(tbl,read.delim(paste0("david/ofi",i,".txt"),sep="\t",stringsAsFactors=F,check.names=F))
}
ind=grep("KEGG_PATHWAY",tbl[,"Category"])
ind=ind[tbl[ind,"FDR"]<cutoff_fdr]
id=sub("hsa","",sub("\\:.*","",tbl[ind,"Term"]))
kegg=sub(".*\\:","",tbl[ind,"Term"])
for(i in 1:length(id)){ # find major categories
    tmp_num=which(map_kegg==id[i])
    for(j in 1:length(cate)){
        if(tmp_num>ind_cate[j] && tmp_num<=ind_cate[j+1]){
            names(kegg)[i]=sub("[0-9]\\. ","",cate[j])
        }
    }
}
tbl2=tbl;ind2=ind;okegg=kegg

ymax=2.7
df=NULL
pos=0
for(i in 1:length(ind1)){
    gene=unlist(strsplit(tbl1[ind1[i],"Genes"],split=", "))
    df=rbind(df,data.frame(kegg=names(bkegg)[i],x1=bfi[gene],x2=bfi[gene],y1=pos,y2=pos+0.1))
    pos=pos+0.15
}
df_kegg=data.frame(kegg=bkegg,x=rep(0,length(bkegg)),y=seq(0,pos-0.1,0.15))
bk3=ggplot(tbl, aes(correlation, fill = kegg, colour = kegg)) +
    #geom_density(alpha = 0) +
    scale_fill_manual(values=tmp_col)  +
    scale_colour_manual(values=tmp_col) +
    ylim(0,ymax) +
    xlim(0,0.8) +
    geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2), linetype=1) +
    labs(title="breast",x="feature importance") +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5),legend.position="none") +
    annotate("text", x = rep(0.4,length(bkegg)), y = seq(0,pos-0.1,0.15)+0.05, label=bkegg, hjust=0,size=3)

df=NULL
pos=0
for(i in 1:length(ind2)){
    gene=unlist(strsplit(tbl2[ind2[i],"Genes"],split=", "))
    df=rbind(df,data.frame(kegg=names(okegg)[i],x1=ofi[gene],x2=ofi[gene],y1=pos,y2=pos+0.1))
    pos=pos+0.15
}
df_kegg=data.frame(kegg=okegg,x=rep(0,length(okegg)),y=seq(0,pos-0.1,0.15))
ok3=ggplot(tbl, aes(correlation, fill = kegg, colour = kegg)) +
    #geom_density(alpha = 0) +
    scale_fill_manual(values=tmp_col)  +
    scale_colour_manual(values=tmp_col) +
    ylim(0,ymax) +
    xlim(0,0.8) +
    geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2), linetype=1) +
    labs(title="ovary",x="feature importance") +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5),legend.position="none") +
    annotate("text", x = rep(0.4,length(okegg)), y = seq(0,pos-0.1,0.15)+0.05, label=okegg, hjust=0,size=3)

list_p=list()
list_p[[1]]=b1
list_p[[2]]=o1
list_p[[3]]=bk1
list_p[[4]]=ok1
pdf(file="figure/sub2_final_kegg.pdf",width=12,height=10)
mat_layout=matrix(c(rep(c(rep(1,6),rep(2,6)),2),rep(c(rep(3,6),rep(4,6)),8)),nrow=10,byrow=T)
multiplot(plotlist=list_p,layout = mat_layout)
dev.off()

list_p=list()
list_p[[1]]=b2
list_p[[2]]=o2
list_p[[3]]=bk2
list_p[[4]]=ok2
pdf(file="figure/sub2_delta_kegg.pdf",width=12,height=10)
mat_layout=matrix(c(rep(c(rep(1,6),rep(2,6)),2),rep(c(rep(3,6),rep(4,6)),8)),nrow=10,byrow=T)
multiplot(plotlist=list_p,layout = mat_layout)
dev.off()

list_p=list()
list_p[[1]]=b3
list_p[[2]]=o3
list_p[[3]]=bk3
list_p[[4]]=ok3
pdf(file="figure/sub2_fi_kegg.pdf",width=12,height=6)
mat_layout=matrix(c(rep(c(rep(1,6),rep(2,6)),2),rep(c(rep(3,6),rep(4,6)),4)),nrow=6,byrow=T)
multiplot(plotlist=list_p,layout = mat_layout)
dev.off()

pdf(file="figure/sub2_kegg_legend.pdf",width=4,height=8)
kegg_legend
dev.off()

