writeLines("
title: Rafaela_alpha_beta
author: Chiara Gini
date: 2023-08-31
output: html_document")


## Rscript
## script to normalize the filtered OTU table starting from CSV file

#############################################################################
## This script is mainly meant to be run locally where R packages can more easily be installed/updated
#############################################################################

## SET UP
library("ape")
library("ggplot2")
library("phyloseq")
library("tidyverse")
library("data.table")
library("metagenomeSeq")
library("factoextra")
library("vegan")
library("broom")
library("DT")
library("scales")

## PARAMETERS
HOME <- Sys.getenv("HOME")
prj_folder = file.path(HOME, "Results")
analysis_folder = "RAFAELA_KOMBUCHA"
outdir = file.path(analysis_folder)

repo = file.path(HOME, "milkqua_microbiome")
source(file.path(repo, "r_scripts/support_functions/dist2list.R")) ## from: https://github.com/vmikk/metagMisc/
source(file.path(repo, "r_scripts/support_functions/phyloseq_transform.R")) ## from: https://github.com/vmikk/metagMisc/

## making results folder
if(!file.exists(file.path(prj_folder, analysis_folder, "results"))) dir.create(file.path(prj_folder, analysis_folder, "results"), showWarnings = FALSE)

## calling metadata
fname = file.path(prj_folder, analysis_folder, "mapping_rafaela_kombucha.csv")
metadata <- fread(fname)

##calling otu table 

otu.tab = read.csv("~/Results/RAFAELA_KOMBUCHA/Kombucha_bacteria-fungi_csv.csv")
otu.tab = subset (otu.tab, Kingdom == "Bacteria")
tax.tab <- otu.tab [, c(1,14)]
row.names(tax.tab) <- tax.tab[,1]
tax.tab[,1] <- NULL
tax.tab = as.matrix(tax.tab)
row.names(otu.tab) <- otu.tab[,1]
otu.tab[,1] <- NULL
otu.tab$Bacteria <- NULL
otu.tab$Kingdom <- NULL
otu.tab = as.matrix(otu.tab)
OTU = otu_table(otu.tab, taxa_are_rows = T)
TAX = tax_table(tax.tab)
physeq=phyloseq(OTU, TAX)
#plot_bar(physeq, fill = "Bacteria") #check for functioning of the pyseq object

## Alpha diversity (alpha diversity is calculated on the original count data, not normalised) 

writeLines(" - calculate alpha diversity indices")
otu_tax_sample <- physeq
alpha = estimate_richness(otu_tax_sample, split = TRUE)
alpha$"sample-id" = row.names(alpha)
alpha = relocate(alpha, `sample-id`)
fwrite(x = alpha, file = file.path(prj_folder, analysis_folder, "results", "alpha_bacteria.csv"))

## drawing alpha diversity plots

# fname = file.path(prj_folder, analysis_folder, "results", "alpha.csv")
# alpha <- fread(fname)
alpha = select(alpha, -c(se.chao1, se.ACE))
names(alpha)[1] <- "sample-id"
alpha$`sample-id` <- gsub("\\.", "-", alpha$`sample-id`)

mAlpha <- reshape2::melt(alpha, id.vars = "sample-id", variable.name = "metric", value.name = "value")
mAlpha$group <- metadata$group[match(mAlpha$`sample-id`,metadata$`sample-id`)]

C <- mAlpha %>%
  group_by(metric, group) %>%
  summarize(N=n(),avg=round(mean(value),3)) %>%
  spread(key = metric, value = avg)
write.csv(C, file="~/Results/RAFAELA_KOMBUCHA/results/alpha_bacteria_base.csv")

D <- mAlpha %>%
  group_by(metric) %>%
  do(tidy(lm(value ~ group, data = .))) %>%
  filter(term != "(Intercept)")
write.csv(D, file="~/Results/RAFAELA_KOMBUCHA/results/alpha_bacteria_signif.csv")

datatable(D, options = list(pageLength=100)) %>% 
  formatStyle('p.value', backgroundColor = styleInterval(0.05, c('yellow', 'white')))

mAlpha <- mAlpha %>% group_by(metric) %>%
  mutate(scaled_value = rescale(value, to = c(0,100)))

p <- ggplot(data = D, mapping= aes(x=term, y=p.value))
p <- p + geom_point(aes(color = metric, stroke = 1), position=position_jitter(h=0, w=0.27)) #,shape=metric
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=0.5)
p <- p + geom_hline(yintercept=0.10, linetype="dashed", color = "darkorange", size=0.5)
p <- p + coord_trans(y="log2")
p <- p + scale_y_continuous(breaks=pretty_breaks(n=20)) 
p <- p + scale_y_continuous(breaks = c(0, 0.05, 0.10, 1)) +  theme(axis.text.x = element_text(angle=90))
#p <- p + scale_color_manual(values = c("#f8766d", "#d39200", "#93aa00", "#00ba38", "#00c19f", "#00b9e3", "#619cff", "#db72fb", "#ff61c3"))
p
ggsave(filename = "~/Results/RAFAELA_KOMBUCHA/results/scatterplot_alpha_bacteria_base.png", width = 6, height = 4, plot = p, device = "png")


## Preprocessing: e.g. filtering
## making normalization folder
if(!file.exists(file.path(prj_folder, outdir))) dir.create(file.path(prj_folder, outdir), showWarnings = FALSE)
writeLines(" - CSS normalization")
otu_tax_sample_norm = phyloseq_transform_css(otu_tax_sample, norm = TRUE, log = FALSE)

## add taxonomy to normalised counts
otu_css_norm = base::as.data.frame(otu_table(otu_tax_sample_norm))
otu_css_norm$tax_id = row.names(otu_css_norm)
otu_css_norm <- relocate(otu_css_norm, tax_id)
taxonomy = as.data.frame(tax_table(otu_tax_sample_norm))
taxonomy$tax_id = row.names(taxonomy)
taxonomy <- relocate(taxonomy, tax_id)
otu_css_norm = otu_css_norm %>% inner_join(taxonomy, by = "tax_id")

writeLines(" - writing out the CSS normalized OTU table")
fwrite(x = otu_css_norm, file = file.path(prj_folder, outdir, "results", "otu_norm_CSS_bacteria.csv"))

###############
## distances ##
###############
writeLines(" - beta diversity: distance matrices")
writeLines(" - available distance metrics")
dist_methods <- unlist(distanceMethodList)
print(dist_methods)

## bray-curtis
writeLines(" - calculate Bray-Curtis distances")
distances = distance(otu_tax_sample_norm, method="bray", type = "samples")
iMDS  <- ordinate(otu_tax_sample_norm, "MDS", distance=distances)

writeLines(" - write out distance matrix")
dd = dist2list(distances, tri = FALSE)
dx = spread(dd, key = "col", value = "value")
fwrite(x = dx, file = file.path(prj_folder, analysis_folder, "results", "bray_curtis_distances_bacteria.csv"))

## euclidean
writeLines(" - calculate Euclidean distances")
distances = distance(otu_tax_sample_norm, method="euclidean", type = "samples")
iMDS  <- ordinate(otu_tax_sample_norm, "MDS", distance=distances)

writeLines(" - write out euclidean distance matrix")
dd = dist2list(distances, tri = FALSE)
dx = spread(dd, key = "col", value = "value")
fwrite(x = dx, file = file.path(prj_folder, analysis_folder, "results", "euclidean_distances_bacteria.csv"))

### add tree ###
random_tree = rtree(ntaxa((otu_tax_sample_norm)), rooted=TRUE, tip.label=taxa_names(otu_tax_sample_norm))
otu_norm_tree = merge_phyloseq(otu_tax_sample_norm, random_tree)

## unifrac
writeLines(" - calculate Unifrac distances")
distances = distance(otu_norm_tree, method="unifrac", type = "samples")
iMDS  <- ordinate(otu_norm_tree, "MDS", distance=distances)

writeLines(" - write out Unifrac distance matrix")
dd = dist2list(distances, tri = FALSE)
dx = spread(dd, key = "col", value = "value")
fwrite(x = dx, file = file.path(prj_folder, analysis_folder, "results", "unifrac_distances_bacteria.csv"))

## weighted unifrac
writeLines(" - calculate weighted Unifrac distances")
distances = distance(otu_norm_tree, method="wunifrac", type = "samples")
iMDS  <- ordinate(otu_norm_tree, "MDS", distance=distances)

writeLines(" - write out weighted Unifrac distance matrix")
dd = dist2list(distances, tri = FALSE)
dx = spread(dd, key = "col", value = "value")
fwrite(x = dx, file = file.path(prj_folder, analysis_folder, "results", "weighted_unifrac_distances_bacteria.csv"))

## DRAWING BETA DIVERSITY PLOT(S)
#take into consideration that the different alternatives we generally use are euclidean, unifrac, weighted unifrac & bray curtis

matrice= read.table("/home/mycelium/Results/RAFAELA_KOMBUCHA/results/bray_curtis_distances_bacteria.csv", row.names = 1, header=T, sep = ",")
vec <- filter(metadata, sample_type == "kombucha") %>% select("sample-id") %>% pull()
vex <- names(matrice) %in% vec
mat_kombucha = matrice[vex,vex]
mat_kombucha$group <- as.character(metadata$group[match(row.names(mat_kombucha),metadata$`sample-id`)])
matx= data.matrix(select(mat_kombucha, -c(group)))

## MDS
mds <- cmdscale(as.dist(matx))
mds <- as.data.frame(mds)
metadata$`sample-id` <- gsub("sample.","sample-", metadata$`sample-id`)
mds$group <- metadata$group[match(rownames(mds), metadata$`sample-id`)]
# mds$sample <- metadata$`sample-id`[match(rownames(mds), metadata$`sample-id`)]

p <- ggplot(mds, aes(V1,V2)) + geom_point(aes(colour = group), size = 3) + stat_ellipse(aes(x=V1, y=V2,color=group), type="t") + xlab("dim1") + ylab("dim2")

res.pca <-prcomp(mds[,-3],  scale = TRUE)
f <- fviz_pca_ind(res.pca, label='none',alpha.ind = 1, habillage=mds$group, repel = TRUE, invisible='quali') + 
  ggforce::geom_mark_ellipse(aes(fill = Groups, color = Groups)) +
  theme(legend.position = 'bottom') + ggtitle("bray curtis distances")
f
fname = file.path("~/Results/RAFAELA_KOMBUCHA/results/beta_bray-curtis_bacteria.png")
ggsave(filename = fname, plot = f, device = "png", dpi = 300, width = 10, height = 7)

## This is basically a PCA plot. The 50% means that the component of the PC(principal component) accounts for 50% of the total variation. The second PC accounts for 50% of the variation. So together they can explain 100% variation in the dataset.

## Significance of between-group distances: significance values based on permuted analysis of variance (999 permutations), repeated 100 times.
pv_group <- replicate(100, adonis2(matx ~ mat_kombucha$group, permutations = 999)$"Pr(>F)"[1], simplify = "vector")
mean(pv_group)
sd(pv_group)



#########
## FUNGI
#########

##calling otu table 

otu.tab = read.csv("~/Results/RAFAELA_KOMBUCHA/Kombucha_bacteria-fungi_csv.csv")
otu.tab = subset (otu.tab, Kingdom == "Fungi")
tax.tab <- otu.tab [, c(1,14)]
row.names(tax.tab) <- tax.tab[,1]
tax.tab[,1] <- NULL
tax.tab = as.matrix(tax.tab)
row.names(otu.tab) <- otu.tab[,1]
otu.tab[,1] <- NULL
otu.tab$Bacteria <- NULL
otu.tab$Kingdom <- NULL
otu.tab = as.matrix(otu.tab)
OTU = otu_table(otu.tab, taxa_are_rows = T)
TAX = tax_table(tax.tab)
physeq=phyloseq(OTU, TAX)
#plot_bar(physeq, fill = "Bacteria") #check for functioning of the pyseq object

## Alpha diversity (alpha diversity is calculated on the original count data, not normalised) 

writeLines(" - calculate alpha diversity indices")
otu_tax_sample <- physeq
alpha = estimate_richness(otu_tax_sample, measures=c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")) #excluded Fisher
alpha$"sample-id" = row.names(alpha)
alpha = relocate(alpha, `sample-id`)
fwrite(x = alpha, file = file.path(prj_folder, analysis_folder, "results", "alpha_fungi.csv"))

## drawing alpha diversity plots

# fname = file.path(prj_folder, analysis_folder, "results", "alpha.csv")
# alpha <- fread(fname)
alpha = select(alpha, -c(se.chao1, se.ACE))
names(alpha)[1] <- "sample-id"
alpha$`sample-id` <- gsub("\\.", "-", alpha$`sample-id`)

mAlpha <- reshape2::melt(alpha, id.vars = "sample-id", variable.name = "metric", value.name = "value")
mAlpha$group <- metadata$group[match(mAlpha$`sample-id`,metadata$`sample-id`)]

C <- mAlpha %>%
  group_by(metric, group) %>%
  summarize(N=n(),avg=round(mean(value),3)) %>%
  spread(key = metric, value = avg)
write.csv(C, file="~/Results/RAFAELA_KOMBUCHA/results/alpha_fungi_base.csv")

D <- mAlpha %>%
  group_by(metric) %>%
  do(tidy(lm(value ~ group, data = .))) %>%
  filter(term != "(Intercept)")
write.csv(D, file="~/Results/RAFAELA_KOMBUCHA/results/alpha_fungi_signif.csv")

datatable(D, options = list(pageLength=100)) %>% 
  formatStyle('p.value', backgroundColor = styleInterval(0.05, c('yellow', 'white')))

mAlpha <- mAlpha %>% group_by(metric) %>%
  mutate(scaled_value = rescale(value, to = c(0,100)))

p <- ggplot(data = D, mapping= aes(x=term, y=p.value))
p <- p + geom_point(aes(color = metric, stroke = 1), position=position_jitter(h=0, w=0.27)) #,shape=metric
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=0.5)
p <- p + geom_hline(yintercept=0.10, linetype="dashed", color = "darkorange", size=0.5)
p <- p + coord_trans(y="log2")
p <- p + scale_y_continuous(breaks=pretty_breaks(n=20)) 
p <- p + scale_y_continuous(breaks = c(0, 0.05, 0.10, 1)) +  theme(axis.text.x = element_text(angle=90))
#p <- p + scale_color_manual(values = c("#f8766d", "#d39200", "#93aa00", "#00ba38", "#00c19f", "#00b9e3", "#619cff", "#db72fb", "#ff61c3"))
p
ggsave(filename = "~/Results/RAFAELA_KOMBUCHA/results/scatterplot_alpha_fungi_base.png", width = 6, height = 4, plot = p, device = "png")


## Preprocessing: e.g. filtering
## making normalization folder
if(!file.exists(file.path(prj_folder, outdir))) dir.create(file.path(prj_folder, outdir), showWarnings = FALSE)
writeLines(" - CSS normalization")
otu_tax_sample_norm = phyloseq_transform_css(otu_tax_sample, norm = TRUE, log = FALSE)

## add taxonomy to normalised counts
otu_css_norm = base::as.data.frame(otu_table(otu_tax_sample_norm))
otu_css_norm$tax_id = row.names(otu_css_norm)
otu_css_norm <- relocate(otu_css_norm, tax_id)
taxonomy = as.data.frame(tax_table(otu_tax_sample_norm))
taxonomy$tax_id = row.names(taxonomy)
taxonomy <- relocate(taxonomy, tax_id)
otu_css_norm = otu_css_norm %>% inner_join(taxonomy, by = "tax_id")

writeLines(" - writing out the CSS normalized OTU table")
fwrite(x = otu_css_norm, file = file.path(prj_folder, outdir, "results", "otu_norm_fungi_CSS.csv"))

###############
## distances ##
###############
writeLines(" - beta diversity: distance matrices")
writeLines(" - available distance metrics")
dist_methods <- unlist(distanceMethodList)
print(dist_methods)

## bray-curtis
writeLines(" - calculate Bray-Curtis distances")
distances = distance(otu_tax_sample_norm, method="bray", type = "samples")
iMDS  <- ordinate(otu_tax_sample_norm, "MDS", distance=distances)

writeLines(" - write out distance matrix")
dd = dist2list(distances, tri = FALSE)
dx = spread(dd, key = "col", value = "value")
fwrite(x = dx, file = file.path(prj_folder, analysis_folder, "results", "bray_curtis_distances_fungi.csv"))

## euclidean
writeLines(" - calculate Euclidean distances")
distances = distance(otu_tax_sample_norm, method="euclidean", type = "samples")
iMDS  <- ordinate(otu_tax_sample_norm, "MDS", distance=distances)

writeLines(" - write out euclidean distance matrix")
dd = dist2list(distances, tri = FALSE)
dx = spread(dd, key = "col", value = "value")
fwrite(x = dx, file = file.path(prj_folder, analysis_folder, "results", "euclidean_distances_fungi.csv"))

### add tree ###
random_tree = rtree(ntaxa((otu_tax_sample_norm)), rooted=TRUE, tip.label=taxa_names(otu_tax_sample_norm))
otu_norm_tree = merge_phyloseq(otu_tax_sample_norm, random_tree)

## unifrac
writeLines(" - calculate Unifrac distances")
distances = distance(otu_norm_tree, method="unifrac", type = "samples")
iMDS  <- ordinate(otu_norm_tree, "MDS", distance=distances)

writeLines(" - write out Unifrac distance matrix")
dd = dist2list(distances, tri = FALSE)
dx = spread(dd, key = "col", value = "value")
fwrite(x = dx, file = file.path(prj_folder, analysis_folder, "results", "unifrac_distances_fungi.csv"))

## weighted unifrac
writeLines(" - calculate weighted Unifrac distances")
distances = distance(otu_norm_tree, method="wunifrac", type = "samples")
iMDS  <- ordinate(otu_norm_tree, "MDS", distance=distances)

writeLines(" - write out weighted Unifrac distance matrix")
dd = dist2list(distances, tri = FALSE)
dx = spread(dd, key = "col", value = "value")
fwrite(x = dx, file = file.path(prj_folder, analysis_folder, "results", "weighted_unifrac_distances_fungi.csv"))

## DRAWING BETA DIVERSITY PLOT(S)
#take into consideration that the different alternatives we generally use are euclidean, unifrac, weighted unifrac & bray curtis

matrice= read.table("/home/mycelium/Results/RAFAELA_KOMBUCHA/results/bray_curtis_distances_fungi.csv", row.names = 1, header=T, sep = ",")
vec <- filter(metadata, sample_type == "kombucha") %>% select("sample-id") %>% pull()
vex <- names(matrice) %in% vec
mat_kombucha = matrice[vex,vex]
mat_kombucha$group <- as.character(metadata$group[match(row.names(mat_kombucha),metadata$`sample-id`)])
matx= data.matrix(select(mat_kombucha, -c(group)))

## MDS
mds <- cmdscale(as.dist(matx))
mds <- as.data.frame(mds)
metadata$`sample-id` <- gsub("sample.","sample-", metadata$`sample-id`)
mds$group <- metadata$group[match(rownames(mds), metadata$`sample-id`)]
# mds$sample <- metadata$`sample-id`[match(rownames(mds), metadata$`sample-id`)]

p <- ggplot(mds, aes(V1,V2)) + geom_point(aes(colour = group), size = 3) + stat_ellipse(aes(x=V1, y=V2,color=group), type="t") + xlab("dim1") + ylab("dim2")

res.pca <-prcomp(mds[,-3],  scale = TRUE)
f <- fviz_pca_ind(res.pca, label='none',alpha.ind = 1, habillage=mds$group, repel = TRUE, invisible='quali') + 
  ggforce::geom_mark_ellipse(aes(fill = Groups, color = Groups)) +
  theme(legend.position = 'bottom') + ggtitle("bray curtis distances")
f
fname = file.path("~/Results/RAFAELA_KOMBUCHA/results/beta_fungi-bray-curtis.png")
ggsave(filename = fname, plot = f, device = "png", dpi = 300, width = 10, height = 7)

## This is basically a PCA plot. The 50% means that the component of the PC(principal component) accounts for 50% of the total variation. The second PC accounts for 50% of the variation. So together they can explain 100% variation in the dataset.

## Significance of between-group distances: significance values based on permuted analysis of variance (999 permutations), repeated 100 times.
pv_group <- replicate(100, adonis2(matx ~ mat_kombucha$group, permutations = 999)$"Pr(>F)"[1], simplify = "vector")
mean(pv_group)
sd(pv_group)

