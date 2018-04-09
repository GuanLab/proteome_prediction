## name: create_final_site.r
## date: 11/13/2017

## Three arguements:
# 1.subset list (Note this version has a header line)
# 2.input file
# 3.name of output file

args <- commandArgs(trailingOnly = TRUE)

final_site=readLines(args[1]) # final site list
final_site=final_site[-1]

x=as.matrix(read.delim(args[2],check.names=F,row.names=1))
y=x[final_site,]

tmp=cbind(rownames(y),y);colnames(tmp)[1]="phosphoID"
write.table(tmp,file=args[3],quote=F,sep="\t",row.names=F)


