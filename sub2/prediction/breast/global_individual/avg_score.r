## name: avg_score.r
## date: 03/06/2018

args <- commandArgs(trailingOnly = TRUE)

x=readLines(args[1])
y=as.numeric(sub("\\[1\\] ","",x))

mat=matrix(y,nrow=5,ncol=length(y)/5)
writeLines(as.character(apply(mat,2,mean)),con="cor_scan.txt")

