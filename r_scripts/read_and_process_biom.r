## Rscript
## script to normalize the filtered OTU table

#############################################################################
## This script is mainly meant to be run locally
## where R packages can more easily be installed/updated
#############################################################################

## SET UP
library("ape")
library("ggplot2")
library("phyloseq")
library("tidyverse")
library("data.table")
library("metagenomeSeq")
source("dada2_etc/dist2list.R") ## from: https://github.com/vmikk/metagMisc/
source("dada2_etc/phyloseq_transform.R") ## from: https://github.com/vmikk/metagMisc/

## PARAMETERS
HOME <- Sys.getenv("HOME")
prj_folder = file.path(HOME, "Documents/MILKQUA")
fname = "dada2_etc/otu_table/otu_table_filtered.biom"
conf_file = "dada2_etc/Config/mapping_milkqua_skinswabs.csv"
min_tot_counts = 500 ## minimum number of total counts per sample to be included in the analysis
outdir = "results"

writeLines(" - reading the filtered (OTU-wise) biom file into phyloseq")
## both the OTU table and the taxonomic classification are available from the biom file (qiime 1.9)
biom_otu_tax <- phyloseq::import_biom(BIOMfilename = file.path(prj_folder,fname))

writeLines(" - removing samples with too few total counts")
biom_otu_tax = prune_samples(sample_sums(biom_otu_tax)>=min_tot_counts, biom_otu_tax)

otu = otu_table(biom_otu_tax, taxa_are_rows = TRUE)
taxa = tax_table(biom_otu_tax)

print(paste("N. of OTUs read from biom file is:", nrow(otu)))
print(paste("N .of samples retained after filtering is:", ncol(otu)))

colnames(otu) <- paste("sample-",colnames(otu),sep="")
print(head(otu))

writeLines(" - change the names of taxonomic levels to Kngdom, Class etc.")
colnames(taxa) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
print(head(taxa))

## metadata
writeLines(" - reading the metadata")
metadata = fread(file.path(prj_folder,conf_file))
names(metadata)[1] <- "sample-id"
if(is.numeric(metadata$`sample-id`)) metadata$`sample-id` = paste("sample",metadata$`sample-id`,sep="-")
metadata <- as.data.frame(metadata)
row.names(metadata) <- metadata$`sample-id`
metadata$`sample-id` <- NULL

## read into phyloseq
writeLines(" - add metadata to the phyloseq object")
samples = sample_data(metadata)
otu_tax_sample = phyloseq(otu,taxa,samples)
sample_data(otu_tax_sample)

# sample_names(otu_tax_sample)
# rank_names(otu_tax_sample)
# taxa_names(otu_tax_sample)

# plot_bar(otu_tax_sample, "Class", fill="Kingdom") + facet_grid(timepoint~treatment)

## Alpha diversity
## alpha diversity is calculated on the original count data, not normalised 
## (see https://www.bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis)
## (see https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531)
writeLines(" - calculate alpha diversity indices")
alpha = estimate_richness(otu_tax_sample, split = TRUE)
fwrite(x = alpha, file = file.path(prj_folder, outdir, "alpha.csv"))
p <- plot_richness(otu_tax_sample, x="treatment", color="timepoint")
ggsave(filename = file.path(prj_folder, outdir, "alpha_plot.png"), plot = p, device = "png", width = 11, height = 7)

## Preprocessing: e.g. filtering
writeLines(" - CSS normalization")
otu_tax_sample_norm = phyloseq_transform_css(otu_tax_sample, norm = TRUE, log = FALSE)
otu_css_norm = otu_table(otu_tax_sample_norm)
writeLines(" - writing out the CSS normalized OTU table")
fwrite(x = otu_css_norm, file = file.path(prj_folder, outdir, "otu_norm_CSS.csv"))

otu_relative = transform_sample_counts(otu_tax_sample, function(x) x/sum(x) )
otu_rel_filtered = filter_taxa(otu_relative, function(x) mean(x) > 5e-3, TRUE)
nrow(otu_table(otu_rel_filtered))

writeLines(" - additionas plots")
random_tree = rtree(ntaxa((otu_rel_filtered)), rooted=TRUE, tip.label=taxa_names(otu_rel_filtered))
plot(random_tree)

biom1 = merge_phyloseq(otu_rel_filtered, random_tree)
plot_tree(biom1, color="treatment", label.tips="taxa_names", ladderize="left", plot.margin=0.3)
plot_tree(biom1, color="Genus", shape="treatment", size="abundance")


## distances
writeLines(" - beta diversity: distance matrices")
writeLines(" - available distance metrics")
dist_methods <- unlist(distanceMethodList)
print(dist_methods)

writeLines(" - calculate Bray-Curtis distances")
distances = distance(otu_tax_sample_norm, method="bray", type = "samples")
iMDS  <- ordinate(otu_tax_sample_norm, "MDS", distance=distances)
p <- plot_ordination(otu_tax_sample_norm, iMDS, color="treatment", shape="timepoint")
ggsave(filename = file.path(prj_folder, outdir, "mds_plot_beta.png"), plot = p, device = "png")

writeLines(" - write out distance matrix")
dd = dist2list(distances, tri = FALSE)
dx = spread(dd, key = "col", value = "value")
fwrite(x = dx, file = file.path(prj_folder, outdir, "bray_curtis_distances.csv"))
