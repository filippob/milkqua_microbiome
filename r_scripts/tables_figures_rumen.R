## script for tables

## libraries
library("scales")
library("ggpubr")
library("ggplot2")
library("ggrepel")
library("reshape2")
library("tidyverse")
library("data.table")
library("gghighlight")


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
fname = file.path(main_path, "tables", "core_microbiota.csv")
select(M, c(taxon,avg)) %>% rename(avg_normalised_counts = avg) %>% fwrite(fname, sep = ",")

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
fname = file.path(main_path, "figures", "Figure1.png")
png(filename = fname, width = 12, height = 7, units = "in", res = 300)
ggdraw() +
  draw_plot(p, x = 0, y = 0.1, width = 0.4, height = 0.75) +
  draw_plot(q, x = 0.4, y = 0, width = 0.6, height = 1) +
  draw_plot_label(label = c("A", "B"), size = 14,
                  x = c(0, 0.5), y = c(1, 1)) 
dev.off()

########################################
## differentially abundant taxa - figure

load("taxonomy_ .RData")
D <- to_save[[1]]
DX <- to_save[[2]]
D0 <- to_save[[3]]

D0 <- mutate(D0, avg_counts = avg_counts+1) %>% spread(key = treatment, value = avg_counts)

D1 <- DX %>%
  inner_join(D0, by = c("level" = "level", "new_taxa" = "new_taxa")) %>%
  mutate(p.value = (p.value)) %>%
  gather(key = "treatment", value = "counts", -c(level,new_taxa, p.value)) ###why was it mutate(p.value = -log10(p.value))?

D1$level <- factor(D1$level, levels = c("phylum","class","order","family","genus"))
D1$treatment <- factor(D1$treatment, levels = c("Control", "NEO", "SEO","Carvacrol", "p-cymene", "g-terpinene"))

D1 <- D1 %>% group_by(level) %>% mutate(tot = sum(counts), relab = counts/tot)

p <- ggplot(D1, aes(x = treatment, y = new_taxa))
p <- p + geom_tile(aes(fill = relab), colour = "white")
p <- p + facet_grid(level~treatment, space = "free", scales = "free")
p <- p + scale_fill_gradient(low = "orange", high = "blue")
p <- p + theme(strip.text.y = element_text(size = 5), 
               strip.text.x = element_text(size = 6),
               # axis.text.y = element_text(size = 4),
               axis.text.x = element_text(size = 6),
               axis.title = element_text(size = 6))
#p <- p + guides(fill="none") 
p <- p + theme(axis.title.y = element_blank(),
                                    axis.text.y = element_blank(),
                                    axis.ticks.y = element_blank())
p <- p + xlab ("Treatments")
p

dd <- filter(D1, treatment == "NEO") %>% mutate(variable = "p-value")

q <- ggplot(dd, aes(x = factor(1), y = new_taxa, group=level))
q <- q + geom_tile(aes(fill = p.value), colour = "white")
q <- q + facet_grid(level~variable, space="free", scales = "free_y")
q <- q + scale_fill_gradient(low = "orange", high = "blue")
q <- q + theme(strip.text = element_text(size = 4), 
               strip.text.x = element_text(size = 6),
               axis.text.y = element_text(size = 5),
               axis.title = element_text(size = 6))
#q <- q + guides(fill=FALSE) 
q <- q + theme(
  # axis.title.x = element_blank(),
  # axis.text.x=element_blank(),
  # axis.ticks.x=element_blank(),
  strip.text.y = element_blank(),
  # axis.text.x = element_blank()
  axis.text.x = element_text(size = 6)
)
q <- q + theme(legend.title = element_text(size = 6)) 
q <- q + theme(legend.text = element_text(size = 6))
q <- q + xlab("") + ylab("Taxa")

q

figure_final <- ggarrange(q, p, widths=c(0.5, 1), labels = "AUTO", common.legend = TRUE, legend = "left",  ncol = 2, nrow = 1, hjust = c(-1, +0.4), vjust = c(1, 1))

print(figure_final)
ggsave(filename = "heatmap_rumen.png", plot = figure_final, device = "png", width = 8, height = 5)

## differentially abundant taxa - table

load("taxonomy_ .RData")
D <- to_save[[1]]
DX <- to_save[[2]]
D0 <- to_save[[3]]

dd <- spread(D0, key = treatment, value = avg_counts)
temp <- inner_join(DX,dd, by = c("level" = "level", "new_taxa" = "new_taxa"))
fwrite(temp, file = "rumen_significant_otus.csv", col.names = TRUE, sep = ",")
print (temp)

################################
### F:B RATIO               ####
################################

## figure
fb_ratio = fread(file.path(main_path, "intermediate_results/fb_ratio_stats.csv"))
load(file.path(main_path,"intermediate_results/fb.RData"))

mO = to_save[[1]] ## F:B ratio from the original data configuration
mR = to_save[[2]] ## F:B ratio from bootstrapping

mR$treatment <- factor(mR$treatment, levels = levels(mO$treatment))

temp <- mO %>%
  mutate(ratio=Firmicutes/Bacteroidetes) 

p <- ggplot(temp, aes(x=treatment,y=ratio)) + geom_boxplot(aes(fill=treatment)) 
# p

q <- ggplot(mR, aes(x=treatment, y=ratio))
q <- q + geom_boxplot(aes(fill=treatment))
# q

g <- ggarrange(p,q,ncol = 2, labels = c("A","B"), common.legend = TRUE)

fname = file.path(main_path, "figures", "fb_ratio.png")
ggsave(filename = fname, plot = g, device = "png", width = 8.5, height = 5, dpi = 300)

## table
dd <- fread(file.path(main_path, "intermediate_results/fb_ratio_stats.csv"))
temp = group_by(filter(mR, !is.na(ratio)), treatment) %>% summarise(boot_avg_FB = mean(ratio), boot_med_FB = median(ratio))

dd %>% left_join(temp, by = "treatment") %>% 
  select(-statistic) %>%
  rename(FB_avg = `F/B_avg`, FB_med = `F/B_med`, estimate_diff = estimate) %>%
  fwrite(file.path(main_path, "tables", "fb_ratio.csv"))


################################
#### ALPHA DIVERSITY          ##
################################
fname <- file.path(main_path, path_to_results,"alpha_diversity/alpha.txt")
alpha <- read.table(fname, header = TRUE)
alpha$sample <- row.names(alpha)
alpha$observed_species <- NULL

mAlpha <- reshape2::melt(alpha, id.vars = "sample", variable.name = "metric", value.name = "value")

mAlpha$treatment <- meta_subset$treatment[match(mAlpha$sample,meta_subset$sample)]

mAlpha$treatment <- factor(mAlpha$treatment, levels = c("Control","NEO","SEO","Carvacrol", "p-cymene", "g-terpinene"))
D <- mAlpha %>%
  group_by(metric,treatment) %>%
  summarize(N=n(),avg=round(mean(value),3)) %>%
  spread(key = metric, value = avg)

fname <- file.path(main_path, "tables","alpha_diversity.csv")
fwrite(D, file = fname, col.names = TRUE)

fname = file.path(main_path, "intermediate_results/alpha_significance.csv")
alpha_significance = fread(fname)

p <- ggplot(data = alpha_significance, mapping= aes(x=term, y=p.value))
p <- p + geom_point(aes(color = metric, stroke = 1.5), position=position_jitter(h=0, w=0.25))
p <- p + geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=0.5)
p <- p + geom_hline(yintercept=0.10, linetype="dashed", color = "darkorange", size=0.5)
p <- p + scale_y_continuous(breaks=pretty_breaks(n=20)) 
p <- p + gghighlight(p.value < 0.10, label_key = metric, use_direct_label = TRUE)
p <- p + coord_trans(y="log2")
p <- p + theme(axis.text.x = element_text(angle = 90))
p

fname = file.path(main_path, "figures", "alpha_significance.png")
ggsave(filename = fname, plot = p, device = "png", width = 7, height = 5, dpi = 300)


#######################################
#### BETA DIVERSITY                 ###
#######################################
# matrice= read.table(file.path(project_folder,"rumen/qiime_1.9/results/beta_diversity/weighted_unifrac_CSS_normalized_otu_table.txt"), row.names=1, header=T)
fname = file.path(main_path, path_to_results, "beta_diversity", "weighted_unifrac_CSS_normalized_otu_table.txt")
matrice= read.table(fname, row.names=1, header=T)

names(matrice) <- gsub("X","",names(matrice))

samples = filter(meta_subset, treatment != "ruminal liquid") %>% pull(sample)
vec <- rownames(matrice) %in% samples
matrice = matrice[vec,vec]

matrice$treatment <- as.character(meta_subset$treatment[match(row.names(matrice),meta_subset$sample)])
matx= data.matrix(select(matrice, -c(treatment)))

## MDS
mds <- cmdscale(as.dist(matx))
mds <- as.data.frame(mds)
mds$treatment <- meta_subset$treatment[match(rownames(mds), meta_subset$sample)]
mds$cow <- meta_subset$cow[match(rownames(mds), meta_subset$sample)]
mds <- mutate(mds, cow = as.factor(cow))

p <- ggplot(mds, aes(V1,V2)) + geom_point(aes(colour = treatment, shape = cow), size = 3)
p <- p + xlab("dim1") + ylab("dim2")
# p

fname = file.path(main_path, "figures", "beta_diversity.png")
ggsave(filename = fname, plot = p, device = "png", dpi = 300, width = 6, height = 5)
