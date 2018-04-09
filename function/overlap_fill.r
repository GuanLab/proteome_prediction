## name: overlap_fill.r
## date: 10/13/2017

## Here I fill NAs in the 1st file with the 2nd file

args <- commandArgs(trailingOnly = TRUE)

input1=as.matrix(read.delim(args[1],stringsAsFactors=F,check.names=F,row.names=1))
input2=as.matrix(read.delim(args[2],stringsAsFactors=F,check.names=F,row.names=1))

common_row=intersect(rownames(input1),rownames(input2))
common_col=intersect(colnames(input1),colnames(input2))

input1[common_row,common_col][is.na(input1[common_row,common_col])]=input2[common_row,common_col][is.na(input1[common_row,common_col])]

output=input1;output=cbind(rownames(output),output);colnames(output)[1]="proteinID"
write.table(output,file=paste0(args[1],".overlap_fill"),quote=F,sep="\t",row.names=F)


