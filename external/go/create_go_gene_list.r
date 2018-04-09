## name: create_go_gene_list.r
## date: 12/07/2017

## Here I download a tsv table from go website and extract genes accosciated with GO terms
## e.g.
## 1. http://amigo.geneontology.org/amigo/term/GO:0001932
## 2. go to link to all genes and gene products
## 3. filter Organism - Homo sapiens
## 4. Custom DL (1) Gene/product(bioentity) (2) Gene/product(bioentity_label) (3) Gene/product name(bioentity_name)
## 5. save them in go_0001932.txt
## note that a parent entry doesn't nessarily cover a children entry - so download all of related ones

## protein expression:
#0010467:gene expression
#0010468:regulation of gene expression

## protein phosphorylation:
#0001932:regulation of protein phosphorylation
#0006468:protein phosphorylation
#0042325:regulation of phosphorylation
#0016310:phosphorylation

## 1. gene expression
go=c("0010467","0010468")

tbl=NULL
for(i in go){
	tbl=rbind(tbl,as.matrix(read.delim(paste0("download/go_",i,".txt"),header=F)))
}

gene=sort(unique(tbl[,2]))
writeLines(gene,con="gene_go_gene_expression.txt")

## 2. protein phosphorylation
go=c("0001932","0006468","0042325","0016310")

tbl=NULL
for(i in go){
        tbl=rbind(tbl,as.matrix(read.delim(paste0("download/go_",i,".txt"),header=F)))
}

gene=sort(unique(tbl[,2]))
writeLines(gene,con="gene_go_protein_phosphorylation.txt")






