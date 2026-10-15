#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE, scipen=999)
projects <- if (.Platform$OS.type=="windows") "//wsl.localhost/Ubuntu/home/dodi/projects" else file.path(Sys.getenv("HOME"),"projects")
root <- file.path(projects,"qe_cross_species_discovery","05_host_focused_DESeq2_sensitivity")
plant <- file.path(projects,"qe_plant_denovo")
res <- read.delim(file.path(root,"results","plant_gene_DESeq2_GP_vs_UP_host_focused_all.tsv"),check.names=FALSE)
out <- file.path(root,"enrichment"); dir.create(out,recursive=TRUE,showWarnings=FALSE)
if(any(file.exists(file.path(out,c("host_focused_GO_GP_up.tsv","host_focused_KEGG_GP_up.tsv"))))) stop("Host-focused enrichment outputs already exist; no files changed.")
if(!all(c("gene_id","padj","DESeq2_status_host_focused")%in%names(res))) stop("Host result columns missing")
valid <- res$gene_id[!is.na(res$padj) & res$padj!="NA"]
gp <- res$gene_id[res$DESeq2_status_host_focused=="GP_up"]
up <- res$gene_id[res$DESeq2_status_host_focused=="UP_up"]
run_ora <- function(map, foreground) {
 map <- unique(map[map$gene_id %in% valid & nzchar(map$term_id) & map$term_id!="-",c("gene_id","term_id")])
 bg <- unique(map$gene_id); fg <- intersect(unique(foreground),bg)
 sets <- lapply(split(map$gene_id,map$term_id),unique); sets <- sets[lengths(sets)>=5]
 ans <- do.call(rbind,lapply(names(sets),function(t) { s<-sets[[t]]; a<-length(intersect(fg,s)); b<-length(fg)-a; c<-length(s)-a; d<-length(bg)-a-b-c; z<-fisher.test(matrix(c(a,b,c,d),2,byrow=TRUE),alternative="greater"); data.frame(term_id=t,foreground_genes_with_term=a,foreground_gene_set_size=length(fg),background_genes_with_term=length(s),background_gene_set_size=length(bg),fold_enrichment=if(a) (a/length(fg))/(length(s)/length(bg)) else 0,odds_ratio=unname(z$estimate),pvalue=z$p.value) }))
 ans$padj <- p.adjust(ans$pvalue,"BH"); ans[order(ans$padj,ans$pvalue),]
}
go <- read.delim(file.path(plant,"04_enrichment","primary_FDR001_LFC1","go","plant_tested_gene_to_GO.tsv"),check.names=FALSE)
kegg <- read.delim(file.path(plant,"04_enrichment","primary_FDR001_LFC1","kegg","plant_tested_gene_to_KEGG_pathway_canonical.tsv"),check.names=FALSE)
# Preserve original canonical-pathway correction.
kegg$term_id[kegg$term_id=="map00072"] <- "map00650"; kegg$term_id[kegg$term_id=="map01130"] <- "map01110"
lookup <- read.delim(file.path(plant,"04_enrichment","kegg","KEGG_reference_pathway_lookup.tsv"),check.names=FALSE)
for (x in c("GP_up","UP_up")) {
 fg <- if(x=="GP_up") gp else up
 g <- run_ora(go,fg); k <- run_ora(kegg,fg)
 g$significant_BH_0.05 <- g$padj<.05
 k$KEGG_pathway_name <- lookup$KEGG_pathway_name[match(k$term_id,lookup$KEGG_pathway_ID)]
 k$significant_BH_0.05 <- k$padj<.05
 write.table(g,file.path(out,paste0("host_focused_GO_",x,".tsv")),sep="\t",quote=FALSE,row.names=FALSE)
 write.table(k,file.path(out,paste0("host_focused_KEGG_",x,".tsv")),sep="\t",quote=FALSE,row.names=FALSE)
}
cat("HOST-FOCUSED GO/KEGG ENRICHMENT: PASS\n")
