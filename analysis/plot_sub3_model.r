## name: plot_sub3_model.r
## date: 01/31/2018


library(ggplot2)
library(reshape2)
library(plyr)
library(scales)

source("~/function/my_palette.r")
source("~/function/multiplot.R")

## 1. number of features #########
b10=readLines("../sub3/prediction/breast/individual_10feature/cor_avg_breast.txt")
b10=as.numeric(sub("\\[1\\] ","",b10))
b100=readLines("../sub3/prediction/breast/individual_100feature/cor_avg_breast.txt")
b100=as.numeric(sub("\\[1\\] ","",b100))
b1000=readLines("../sub3/prediction/breast/individual_1000feature/cor_avg_breast.txt")
b1000=as.numeric(sub("\\[1\\] ","",b1000))
b=readLines("../sub3/prediction/breast/individual/cor_avg_breast.txt")
b=as.numeric(sub("\\[1\\] ","",b))
bgo=readLines("../sub3/prediction/breast/individual_go_phosphorylation/cor_avg_breast.txt")
bgo=as.numeric(sub("\\[1\\] ","",bgo))

o10=readLines("../sub3/prediction/ova/individual_10feature/cor_avg_ova.txt")
o10=as.numeric(sub("\\[1\\] ","",o10))
o100=readLines("../sub3/prediction/ova/individual_100feature/cor_avg_ova.txt")
o100=as.numeric(sub("\\[1\\] ","",o100))
o1000=readLines("../sub3/prediction/ova/individual_1000feature/cor_avg_ova.txt")
o1000=as.numeric(sub("\\[1\\] ","",o1000))
o=readLines("../sub3/prediction/ova/individual/cor_avg_ova.txt")
o=as.numeric(sub("\\[1\\] ","",o))
ogo=readLines("../sub3/prediction/ova/individual_go_phosphorylation/cor_avg_ova.txt")
ogo=as.numeric(sub("\\[1\\] ","",ogo))

tbl=rbind(cbind(b10,b100,b1000,bgo,b),cbind(o10,o100,o1000,ogo,o))
colnames(tbl)=c(10,100,1000,"1476/1079(GO)","6956/3217(all)")
rownames(tbl)=rep(c("breast","ovary"),each=5)
tbl.m=melt(tbl)
colnames(tbl.m)=c("cancer","model","correlation")
tmp_col=my_palette[c("blue","red")]
names(tmp_col)=c("breast","ovary")

sub3_feature = ggplot(tbl.m, aes(model, correlation, fill=cancer)) +
    geom_violin(colour="grey50",draw_quantiles = 0.5) +
    scale_fill_manual(values=tmp_col) +
    labs(title="Sub3 Individual Model", x="Number of features", y="Pearson's correlation") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylim(0,0.75) +
    geom_rect(mapping=aes(xmin=4.5,xmax=5.5,ymin=0.25,ymax=0.7,fill=NA),color=my_palette["yellow"])
#######################################

## 2. number of samples ###############
b0.4=readLines("../sub3/prediction/breast/individual_0.4sample/cor_avg_breast.txt")
b0.4=as.numeric(sub("\\[1\\] ","",b0.4))
b0.6=readLines("../sub3/prediction/breast/individual_0.6sample/cor_avg_breast.txt")
b0.6=as.numeric(sub("\\[1\\] ","",b0.6))
b0.8=readLines("../sub3/prediction/breast/individual_0.8sample/cor_avg_breast.txt")
b0.8=as.numeric(sub("\\[1\\] ","",b0.8))
b=readLines("../sub3/prediction/breast/individual/cor_avg_breast.txt")
b=as.numeric(sub("\\[1\\] ","",b))
bt=readLines("../sub3/prediction/breast/individual_transplant/cor_avg_breast.txt")
bt=as.numeric(sub("\\[1\\] ","",bt))

o0.4=readLines("../sub3/prediction/ova/individual_0.4sample/cor_avg_ova.txt")
o0.4=as.numeric(sub("\\[1\\] ","",o0.4))
o0.6=readLines("../sub3/prediction/ova/individual_0.6sample/cor_avg_ova.txt")
o0.6=as.numeric(sub("\\[1\\] ","",o0.6))
o0.8=readLines("../sub3/prediction/ova/individual_0.8sample/cor_avg_ova.txt")
o0.8=as.numeric(sub("\\[1\\] ","",o0.8))
o=readLines("../sub3/prediction/ova/individual/cor_avg_ova.txt")
o=as.numeric(sub("\\[1\\] ","",o))
ot=readLines("../sub3/prediction/ova/individual_transplant/cor_avg_ova.txt")
ot=as.numeric(sub("\\[1\\] ","",ot))

tbl=rbind(cbind(b0.4,b0.6,b0.8,b,bt),cbind(o0.4,o0.6,o0.8,o,ot))
colnames(tbl)=c("40%","60%","80%","100%(b=105;o=69)","100% b+o")
rownames(tbl)=rep(c("breast","ovary"),each=5)
tbl.m=melt(tbl)
colnames(tbl.m)=c("cancer","model","correlation")
tmp_col=my_palette[c("blue","red")]
names(tmp_col)=c("breast","ovary")

sub3_sample = ggplot(tbl.m, aes(model, correlation, fill=cancer)) +
    geom_violin(colour="grey50",draw_quantiles = 0.5) +
    scale_fill_manual(values=tmp_col) +
    labs(title="Sub3 Trans-cell Line Model", x="Number of samples", y="Pearson's correlation") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylim(0,0.75) +
    geom_rect(mapping=aes(xmin=4.5,xmax=5.5,ymin=0.25,ymax=0.7,fill=NA),color=my_palette["yellow"])
######################################

## 3. stacking #######################
bg=readLines("../sub3/prediction/breast/proteome/cor_avg_breast.txt")
bg=as.numeric(sub("\\[1\\] ","",bg))
bi=readLines("../sub3/prediction/breast/global_individual/cor_avg_breast.txt")
bi=as.numeric(sub("\\[1\\] ","",bi))
bt=readLines("../sub3/prediction/breast/global_individual_transplant/cor_avg_breast.txt")
bt=as.numeric(sub("\\[1\\] ","",bt))
bf=readLines("../sub3/prediction/breast/final/cor_avg_breast.txt")
bf=as.numeric(sub("\\[1\\] ","",bf))

og=readLines("../sub3/prediction/ova/proteome/cor_avg_ova.txt")
og=as.numeric(sub("\\[1\\] ","",og))
oi=readLines("../sub3/prediction/ova/global_individual/cor_avg_ova.txt")
oi=as.numeric(sub("\\[1\\] ","",oi))
ot=readLines("../sub3/prediction/ova/global_individual_transplant/cor_avg_ova.txt")
ot=as.numeric(sub("\\[1\\] ","",ot))
of=readLines("../sub3/prediction/ova/final/cor_avg_ova.txt")
of=as.numeric(sub("\\[1\\] ","",of))

tbl=rbind(cbind(bg,bi,bt,bf),cbind(og,oi,ot,of))
colnames(tbl)=c("Generic","Gene-specific","Trans-cell line","Multisite")
rownames(tbl)=rep(c("breast","ovary"),each=5)
tbl.m=melt(tbl)
colnames(tbl.m)=c("cancer","model","correlation")
tmp_col=my_palette[c("blue","red")]
names(tmp_col)=c("breast","ovary")

sub3_stacking = ggplot(tbl.m, aes(model, correlation, fill=cancer)) +
    geom_violin(colour="grey50",draw_quantiles = 0.5) +
    scale_fill_manual(values=tmp_col) +
    labs(title="Sub3 Model Comparison", x="Model", y="Pearson's correlation") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylim(0,0.75) +
    geom_rect(mapping=aes(xmin=3.5,xmax=4.5,ymin=0.30,ymax=0.72,fill=NA),color=my_palette["yellow"])
###################################################


list_p=list()
list_p[[1]]=sub3_feature
list_p[[2]]=sub3_sample
list_p[[3]]=sub3_stacking

pdf(file="../figure/sub3_model_comparison.pdf",width=9.5,height=12)
mat_layout=matrix(1:3,nrow=3)
multiplot(plotlist=list_p,layout = mat_layout)
dev.off()



