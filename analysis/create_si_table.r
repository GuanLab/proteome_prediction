## name: create_si_table.r
## date: 04/09/2018


breast=cbind(sub("\\[1\\] ","",readLines("../sub2/prediction/breast/individual_10feature/cor_avg_breast.txt")),
    sub("\\[1\\] ","",readLines("../sub2/prediction/breast/individual_100feature/cor_avg_breast.txt")),
    sub("\\[1\\] ","",readLines("../sub2/prediction/breast/individual_1000feature/cor_avg_breast.txt")),
    sub("\\[1\\] ","",readLines("../sub2/prediction/breast/individual_go_expression/cor_avg_breast.txt")),
    sub("\\[1\\] ","",readLines("../sub2/prediction/breast/individual/cor_avg_breast.txt")))

ova=cbind(sub("\\[1\\] ","",readLines("../sub2/prediction/ova/individual_10feature/cor_avg_ova.txt")),
    sub("\\[1\\] ","",readLines("../sub2/prediction/ova/individual_100feature/cor_avg_ova.txt")),
    sub("\\[1\\] ","",readLines("../sub2/prediction/ova/individual_1000feature/cor_avg_ova.txt")),
    sub("\\[1\\] ","",readLines("../sub2/prediction/ova/individual_go_expression/cor_avg_ova.txt")),
    sub("\\[1\\] ","",readLines("../sub2/prediction/ova/individual/cor_avg_ova.txt")))

tbl=rbind(breast,ova)
rownames(tbl)=paste0("cross-validation",c(1:5,1:5))
colnames(tbl)=c("top10","top100","top1000","GO","all")

write.csv(tbl,file="table/SI_table1.csv",quote=F)


breast=cbind(sub("\\[1\\] ","",readLines("../sub2/prediction/breast/individual_0.4sample/cor_avg_breast.txt")),
    sub("\\[1\\] ","",readLines("../sub2/prediction/breast/individual_0.6sample/cor_avg_breast.txt")),
    sub("\\[1\\] ","",readLines("../sub2/prediction/breast/individual_0.8sample/cor_avg_breast.txt")),
    sub("\\[1\\] ","",readLines("../sub2/prediction/breast/individual/cor_avg_breast.txt")),
    sub("\\[1\\] ","",readLines("../sub2/prediction/breast/individual_transplant/cor_avg_breast.txt")))

ova=cbind(sub("\\[1\\] ","",readLines("../sub2/prediction/ova/individual_0.4sample/cor_avg_ova.txt")),
    sub("\\[1\\] ","",readLines("../sub2/prediction/ova/individual_0.6sample/cor_avg_ova.txt")),
    sub("\\[1\\] ","",readLines("../sub2/prediction/ova/individual_0.8sample/cor_avg_ova.txt")),
    sub("\\[1\\] ","",readLines("../sub2/prediction/ova/individual/cor_avg_ova.txt")),
    sub("\\[1\\] ","",readLines("../sub2/prediction/ova/individual_transplant/cor_avg_ova.txt")))

tbl=rbind(breast,ova)
rownames(tbl)=paste0("cross-validation",c(1:5,1:5))
colnames(tbl)=c("40%sample","60%sample","80%sample","100%sample","trans-tissue")

write.csv(tbl,file="table/SI_table2.csv",quote=F)


breast=cbind(sub("\\[1\\] ","",readLines("../sub2/prediction/breast/rna/cor_avg_breast.txt")),
    sub("\\[1\\] ","",readLines("../sub2/prediction/breast/global_individual/cor_avg_breast.txt")),
    sub("\\[1\\] ","",readLines("../sub2/prediction/breast/final/cor_avg_breast.txt")))

ova=cbind(sub("\\[1\\] ","",readLines("../sub2/prediction/ova/rna/cor_avg_ova.txt")),
    sub("\\[1\\] ","",readLines("../sub2/prediction/ova/global_individual/cor_avg_ova.txt")),
    sub("\\[1\\] ","",readLines("../sub2/prediction/ova/final/cor_avg_ova.txt")))

tbl=rbind(breast,ova)
rownames(tbl)=paste0("cross-validation",c(1:5,1:5))
colnames(tbl)=c("generic","gene-specific","trans-tissue")

write.csv(tbl,file="table/SI_table3.csv",quote=F)


