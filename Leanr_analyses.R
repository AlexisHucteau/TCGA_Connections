library(LEANR)
library(stringr)
library(igraph)
library(dplyr)
library(clusterProfiler)
library(hugene20sttranscriptcluster.db)
library(ggplot2)


FIsInGene_020720_with_annotations <- read.delim("~/Genes_network/Genes_Expansion/Networks/FIsInGene_020720_with_annotations.tsv")
reactom_network <- graph_from_data_frame(FIsInGene_020720_with_annotations)
IDH <- read.csv("GitHub/TCGA_Connections/Results/IDH_WT_differencial_gene_expression.csv")

IDHm_pvalues <- unique(IDH[!duplicated(IDH$hgnc_symbol),c(6,9)])

gene_pvalues <- IDHm_pvalues$P.Value
names(gene_pvalues) <- IDHm_pvalues$hgnc_symbol

lean <- run.lean(gene_pvalues, reactom_network, ncores = 8)

lean[["restab"]] %>% as.data.frame(.) %>% dplyr::filter(., PLEAN < 0.1) %>% dplyr::select(., PLEAN) %>% unique(.) %>% rownames(.) %>% length(.)

res <- lean[["restab"]] %>% as.data.frame(.) %>% dplyr::filter(., PLEAN < 0.1) %>% dplyr::filter(., !duplicated(.$PLEAN))

filter_strong_networks <- lean$nhs %>% names(.) %>% as.character(.) %>% sapply(., function(x) {x %in% rownames(res)}) %>% .[.] %>% names(.)
filter_strong_networks

strong_network <- list()
for (i in filter_strong_networks) {
  strong_network[[i]] <- lean$nhs[[i]]
  if ("ACTB" %in% lean$nhs[[i]]) {
    print("oui")
    print(i)
  }
    
}
ego <- list()
for (i in names(strong_network)) {
  ego[[i]] <- enrichGO(
    strong_network[[i]],
    OrgDb = "org.Hs.eg.db",
    keyType = "SYMBOL",
    ont = "ALL"
  )
}

genes_list <- unlist(strong_network) %>% as.character(.) %>% unique(.)

dotplot(ego[[]], showCategory = 100 ) + ggtitle("GO")

ego2 <- enrichGO(
  genes_list,
  OrgDb = "org.Hs.eg.db",
  keyType = "SYMBOL",
  ont = "ALL"
)
