library(GSAR)
library(GSVAdata)
proteomes <- read.csv("C:/Users/koush/OneDrive/Desktop/Thesis/other datasets/final_proteomes.csv")
data("c2BroadSets")

library(org.Hs.eg.db)
library(GSEABase)

#proteomes2 <- proteomes[c(1:len), c(1,2:20,34:56)] (Basal & Luminal-A)

rownames(proteomes) <- proteomes[,1]
proteomes <- proteomes[,-1]
View(proteomes)
dim(proteomes)

proteomes <- data.matrix(proteomes)


C2 <- as.list(geneIds(c2BroadSets))
len <- length(C2)
genes.entrez <- unique(unlist(C2))
genes.symbol <- array("",c(length(genes.entrez),1))
x <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])
for (ind in 1:length(genes.entrez)){
  if (length(xx[[genes.entrez[ind]]])!=0)
    genes.symbol[ind] <- xx[[genes.entrez[ind]]]
}
## discard genes with no mapping to gene symbol identifiers
genes.no.mapping <- which(genes.symbol == "")
if(length(genes.no.mapping) > 0){
  genes.entrez <- genes.entrez[-genes.no.mapping]
  genes.symbol <- genes.symbol[-genes.no.mapping]
}
names(genes.symbol) <- genes.entrez

##discard genes in C2 pathways which do not exist in proteomes dataset
overlap <- rownames(proteomes)


remained <- array(0,c(1,len))
for (k in seq(1, len, by=1)) {
  remained[k] <- sum((genes.symbol[C2[[k]]] %in% overlap) & 
                       (C2[[k]] %in% genes.entrez))
}

## discard C2 pathways which have less than 10 or more than 500 genes
C2 <- C2[(remained>=10)&&(remained<=500)]
pathway.names <- names(C2)
c2.pathways <- list()
for (k in seq(1, length(C2), by=1)){
  selected.genes <- which(overlap %in% genes.symbol[C2[[k]]])
  c2.pathways[[length(c2.pathways)+1]] <- overlap[selected.genes]
}
names(c2.pathways) <- pathway.names
path.index <- which(names(c2.pathways) == "PYEON_HPV_POSITIVE_TUMORS_UP")


target.pathway <- proteomes[c2.pathways[["PYEON_HPV_POSITIVE_TUMORS_UP"]],]
target.pathway <- data.matrix(target.pathway)
View(target.pathway)
group.label <- c(rep(1,19), rep(2,61))
#dim(target.pathway)

WW_pvalue <- WWtest(target.pathway, group.label)
KS_pvalue <- KStest(target.pathway, group.label)
MD_pvalue <- MDtest(target.pathway, group.label)
RKS_pvalue <- RKStest(target.pathway, group.label)
RMD_pvalue <- RMDtest(target.pathway, group.label)
F_pvalue <- AggrFtest(target.pathway, group.label)
GSNCA_pvalue <- GSNCAtest(target.pathway, group.label)

WW_pvalue
KS_pvalue
MD_pvalue
RKS_pvalue
RMD_pvalue
F_pvalue
GSNCA_pvalue

plotMST2.pathway(object=proteomes[c2.pathways[[path.index]],],
                 group=c(rep(1,19), rep(2,61)), name="PYEON_HPV_POSITIVE_TUMORS_UP", 
                 legend.size=1.2, leg.x=-1.2, leg.y=2, 
                 label.size=1, label.dist=0.8, cor.method="pearson")

results <- TestGeneSets(object=proteomes, group=group.label, 
                        geneSets=c2.pathways[1:3], min.size=10, max.size=100, test="WWtest")
results
