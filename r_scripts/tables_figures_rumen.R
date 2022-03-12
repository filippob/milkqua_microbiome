## script for tables

## libraries
library("ggplot2")
library("tidyverse")
library("data.table")

## path parameters
main_path = "~/Documents/MILKQUA/rumen"
path_to_results = "qiime_1.9/results_48"

## metadata
#project_folder = "~/Documents/MILKQUA"
metadata <- readxl::read_xlsx(file.path(main_path, "mapping_file_rumen.xlsx"), sheet = 1)
names(metadata)[1] <- "sample"

metadata$treatment[ which(metadata$treatment == "no treatment (ruminal liquid + diet)")] <- "Control"
metadata$treatment[ which(metadata$treatment == "AE1")] <- "NEO"
metadata$treatment[ which(metadata$treatment == "AE sintético 1")] <- "SEO"
metadata$treatment[ which(metadata$treatment == "γ-terpinene")] <- "g-terpinene"

meta_subset <- filter(metadata, treatment !="ruminal liquid")
# meta_subset %>%
#   group_by(treatment) %>%
#   dplyr::summarise(N=n())

## OTU - phylum
otu <- fread(file.path(main_path, path_to_results, "taxa_summary_abs/CSS_normalized_otu_table_L2.txt"), header = TRUE, skip = 1)
otu$`#OTU ID` <- gsub("^.*;","",otu$"#OTU ID")
otu <- gather(otu, key = "sample", value ="counts", -`#OTU ID`) %>% spread(key = `#OTU ID`, value = counts)
otu$treatment = metadata$treatment[match(otu$sample,metadata$sample)]
otu <- filter(otu, sample %in% meta_subset$sample)


## relative abundances
metadata_cols = names(meta_subset)[c(1,5)]
M <- dplyr::select(otu,-all_of(metadata_cols))
M <- M/rowSums(M)
M <- bind_cols(dplyr::select(otu, all_of(metadata_cols)),M)

M <- subset(M, treatment !="ruminal liquid")

## plot of phylum abundance
mm <- gather(M, key = "phylum", value = "abundance", -c(sample,treatment))
phyls = group_by(mm, phylum) %>% summarise(avg = mean(abundance)) %>% arrange(desc(avg))
oldc <- phyls$phylum[phyls$avg < 0.01]
newc <- rep("Other", length(oldc))
vec <- newc[match(mm$phylum,oldc)]
mm$phylum <- ifelse(mm$phylum %in% oldc, "Other", as.character(mm$phylum))

mm$phylum <- factor(mm$phylum, levels = c(phyls$phylum[1: (length(phyls$phylum) - length(oldc))],"Other"))

# p <- ggplot(mm, aes( x = phylum, y = abundance)) + geom_boxplot(aes(fill = phylum))
# p <- p + theme(text = element_text(size = 6),
#                axis.text.x = element_text(angle=90))
# p

library("ggpubr")
p <- ggboxplot(mm, "phylum", "abundance", color = "phylum", legend = "none")
p <- p + rotate_x_text(90) + font("xy.text", size=8)
p

## OTU - everything
otu <- fread(file.path(main_path, path_to_results, "taxa_summary_abs/CSS_normalized_otu_table_L6.txt"), header = TRUE, skip = 1)
otu$`#OTU ID` <- gsub("^.*;","",otu$"#OTU ID")
M <- otu[,-1] > 0
vec <- rowSums(M)/ncol(M) > 0.99
M <- otu[vec,]
vec <- !grepl("uncultured", M$`#OTU ID`)
M <- M[vec,]
M$avg <- rowMeans(M[,-1])
M <- arrange(M, desc(avg)) %>% rename(taxon = `#OTU ID`)

## write out table of the core microbiota
select(M, c(taxon,avg)) %>% rename(avg_normalised_counts = avg) %>% fwrite("../tables/core_microbiota.csv", sep = ",")

oldc <- M$taxon[M$avg < 25]
newc <- rep("Other", length(oldc))
vec <- newc[match(M$taxon,oldc)]
M$taxon <- ifelse(M$taxon %in% oldc, "Other", as.character(M$taxon))
M$taxon <- gsub("group","",M$taxon)

M <- group_by(M, taxon) %>% summarise(avg = mean(avg)) %>% arrange(desc(avg))
M$taxon -> M$short_name
M$short_name <- substr(M$short_name,start = 1, stop = 17)
M$short_name <- factor(M$short_name, levels = M$short_name)

require('RColorBrewer')
mycolors = c(brewer.pal(name="Set3", n = 11), brewer.pal(name="Paired", n = 12))

q <- ggpie(M, "avg", label="short_name", color = "white", fill = "short_name",  legend = "none",
      lab.pos = "out", palette = mycolors, ggtheme = theme_pubr()) 
q <- q + font("xy.text", size = 8, color = "gray20", face="bold")
q

g <- ggarrange(p, q, ncol = 2, labels = c("A","B"), heights = c(0.1,4))
# g
# ggsave(filename = "../figures/Figure1.png", plot = g, device = "png", dpi = 250)


library("cowplot")
png(filename = "../figures/Figure1.png", width = 12, height = 7, units = "in", res = 300)
ggdraw() +
  draw_plot(p, x = 0, y = 0.1, width = 0.4, height = 0.75) +
  draw_plot(q, x = 0.4, y = 0, width = 0.6, height = 1) +
  draw_plot_label(label = c("A", "B"), size = 14,
                  x = c(0, 0.5), y = c(1, 1)) 
dev.off()
