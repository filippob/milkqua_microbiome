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
# main_path = "~/Documents/MILKQUA/STOOLS"
# path_to_results = "qiime_1.9/results_48"
basedir="~/Results"
prjdir = "STOOLS"
outdir = "results"

## metadata
#project_folder = "~/Documents/MILKQUA"
# metadata <- readxl::read_xlsx(file.path(main_path, "mapping_file_STOOLS.xlsx"), sheet = 1)
fname = file.path(basedir, prjdir, "mapping_milkqua_stools.csv")
metadata <- fread(fname)
names(metadata)[1] <- "sample"
names(metadata)[9] <- "cow"

# metadata %>%
#   group_by(treatment) %>%
#   dplyr::summarise(N=n())

## OTU - Genus
# otu <- fread(file.path(main_path, path_to_results, "taxa_summary_abs/CSS_normalized_otu_table_L2.txt"), header = TRUE, skip = 1)
otu <- fread(file.path(basedir, prjdir, outdir, "otu_norm_CSS.csv"), header = TRUE)
otu <- filter(otu, Family !="Mitochondria")
otu <- filter(otu, Class !="Chloroplast")
otu = select(otu, -1, -34:-38, -40)
otu <- otu %>% dplyr::select("Genus", everything())
otu <- otu %>% group_by(Genus) %>% summarise_all(funs(sum))

uncult <- slice(otu, 1, 230:238)
uncult$Genus <- "Uncultured or unknown"
uncult <- uncult %>%
  group_by(Genus) %>%
  summarise(across(everything(), sum))

otu <- otu[-c(1, 230:238), ]
otu <- rbind(otu, uncult)

otu <- gather(otu, key = "sample", value ="counts", -Genus) %>% spread(key = Genus, value = counts)

## relative abundances
metadata_cols = names(metadata)[1]
M <- dplyr::select(otu,-all_of(metadata_cols))
M <- M/rowSums(M)
M <- bind_cols(dplyr::select(otu, all_of(metadata_cols)),M)


## plot of genus abundance
mm <- gather(M, key = "genus", value = "abundance", -c(sample))
Genus = group_by(mm, genus) %>% summarise(avg = mean(abundance)) %>% arrange(desc(avg))
oldc <- Genus$genus[Genus$avg < 0.01]
newc <- rep("Lower than 1%", length(oldc))
vec <- newc[match(mm$genus,oldc)]
mm$genus <- ifelse(mm$genus %in% oldc, "Lower than 1%", as.character(mm$genus))

mm$genus <- factor(mm$genus, levels = c(Genus$genus[1: (length(Genus$genus) - length(oldc))],"Lower than 2%"))

fwrite(Genus, file = "~/Results/STOOLS/results/test_abund.csv", sep = ",")

# mm2 <- mm %>% 
#   group_by(genus) %>% 
#   summarise(across(abundance, sum))

# p <- ggplot(mm, aes( x = Genus, y = abundance)) + geom_boxplot(aes(fill = genus))
# p <- p + theme(text = element_text(size = 6),
#                axis.text.x = element_text(angle=90))
# p

library("ggpubr")

require('RColorBrewer')
mycolors = c(brewer.pal(name="Paired", n = 11), brewer.pal(name="Paired", n = 12))
mycolors2 = c(brewer.pal(name="Set2", n = 11), brewer.pal(name="Paired", n = 12))

p <- ggboxplot(mm, "genus", "abundance", color = "genus", legend = "none", palette = mycolors2)
p <- p + rotate_x_text(90) + font("xy.text", size=14)
p

## OTU - everything
# otu <- fread(file.path(main_path, path_to_results, "taxa_summary_abs/CSS_normalized_otu_table_L6.txt"), header = TRUE, skip = 1)
otu <- fread(file.path(basedir, prjdir, outdir, "otu_norm_CSS.csv"), header = TRUE)
otu <- filter(otu, Family !="Mitochondria")
otu <- filter(otu, Class !="Chloroplast")
otu = select(otu, -1, -34:-38, -40) #(otu, -1, -34, -36: -40)
otu <- otu %>%
  group_by(Genus) %>%
  summarise(across(everything(), sum))

uncult <- slice(otu, 1,230:238)
uncult$Genus <- "Uncultured or unknown"
uncult <- uncult %>%
  group_by(Genus) %>%
  summarise(across(everything(), sum))

otu <- otu[-c(1, 230:238), ]
otu <- rbind(otu, uncult)

A <- otu[,-1] > 0
vec <- rowSums(A)/ncol(A) > 0.99
A <- otu[vec,]
vec <- !grepl("uncultured", A$Genus)
A <- A[vec,]
A$avg <- rowMeans(A[,-1])
A <- arrange(A, desc(avg)) %>% rename(taxon = Genus)


## write out table of the core microbiota
ffname = file.path(basedir, prjdir, outdir, "core_microbiota.csv")
select(A, c(taxon,avg)) %>% rename(avg_normalised_counts = avg) %>% fwrite(ffname, sep = ",")


oldc <- A$taxon[A$avg < 160]
newc <- rep("Lower than 1%", length(oldc))
vec <- newc[match(A$taxon,oldc)]
A$taxon <- ifelse(A$taxon %in% oldc, "Lower than 1%", as.character(A$taxon))
A$taxon <- gsub("group","",A$taxon)

A <- group_by(A, taxon) %>% summarise(avg = mean(avg)) %>% arrange(desc(avg))
A$taxon -> A$short_name
A$short_name <- substr(A$short_name,start = 1, stop = 30)
A$short_name <- factor(A$short_name, levels = A$short_name)
names(A)[3] <- "Genera"

labs <- round(100*(A$avg/sum(A$avg)),2)
labs=paste0(labs, "%")

q <- ggpie(A, "avg", label=labs, color = "white", fill = "Genera",  legend = "right",
      lab.pos = "out",  palette = c("red", "green", "khaki", "steelblue1","maroon","lightblue","salmon","turquoise","violet","yellow","gray","tomato","navy","pink","springgreen4","peru","sienna1","plum4","mediumblue","darkorange","brown4","gold","bisque4"), font.legend=c(16, "black"), lab.font = "white", ggtheme = theme_pubr()) 
q <- q + font("xy.text", size = 16, color = "black")
q

A$pct <- A$avg/sum(A$avg)

A$pctper <- round(100*(A$pct),2)
A$pctper = paste0(A$pctper, "%")
A$Genera <- paste0(A$Genera, " ", "(",A$pctper,")")
A <- A %>%
  arrange(desc(pct))
A$Genera <- factor(A$Genera)

p <- ggplot(A, aes(x=factor(1), y=pct, fill=Genera)) + geom_bar(width=1,stat="identity") 
p <- p + coord_polar(theta='y') + guides(fill = guide_legend(title = "Genera")) 
p <- p + xlab("") + ylab("")
# p <- p + geom_label(aes(label=pctper), position = position_stack(0.5), show.legend = FALSE)
p <- p + scale_fill_manual(values = c("red", "green", "khaki", "steelblue1","maroon","lightblue","salmon","turquoise","violet","yellow","gray","tomato","navy","pink","springgreen4","peru","sienna1","plum4","mediumblue","darkorange","brown4","gold","bisque4"))
p <- p + theme(text = element_text(size=20),
               # axis.text.x = element_text(size=20),
               # strip.text = element_text(size = 20),
               axis.text = element_blank(),
               axis.ticks = element_blank(),
               legend.text=element_text(size=20),
               legend.title=element_text(size=20))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
p

# g <- ggarrange(p, q, ncol = 2, labels = c("A","B"), heights= c(1,1), widths = c(0.5,1)) #heights = c(0.1,4))
# g

ggsave(filename = file.path("~/Results/STOOLS/results/Figure2perc.png"), plot = p, device = "png", dpi = 250, width = 21, height = 13)


library("cowplot")
fname = file.path(basedir, prjdir, outdir,"Figure1.png")
png(filename = fname, width = 20, height = 8, units = "in", res = 300)
ggdraw() +
  draw_plot(p, x = 0, y = 0, width = 0.3, height = 0.9) +
  draw_plot(q, x = 0.14, y = 0, width = 1, height = 1) +
  draw_plot_label(label = c("A", "B"), size = 16,
                  x = c(0, 0.32), y = c(1, 1)) 
dev.off()

########################################
## differentially abundant taxa - figure

load("~/Results/STOOLS/results/taxonomy_ .RData")
D <- to_save[[1]]
DX <- to_save[[2]]
D0 <- to_save[[3]]

D0 <- mutate(D0, avg_counts = avg_counts+1) %>% spread(key = treatment, value = avg_counts)

D1 <- DX %>%
  inner_join(D0, by = c("level" = "level", "new_taxa" = "new_taxa")) %>%
  mutate(p.value = (p.value)) %>%
  gather(key = "treatment", value = "counts", -c(level,new_taxa, p.value)) ###why was it mutate(p.value = -log10(p.value))?

D1$level <- factor(D1$level, levels = c("Genus","class","order","family","Genus"))
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
ggsave(filename = "heatmap_STOOLS.png", plot = figure_final, device = "png", width = 8, height = 5)

## differentially abundant taxa - table

load("~/Results/STOOLS/results/taxonomy_ .RData")
D <- to_save[[1]]
DX <- to_save[[2]]
D0 <- to_save[[3]]

dd <- spread(D0, key = treatment, value = avg_counts)
temp <- inner_join(DX,dd, by = c("level" = "level", "new_taxa" = "new_taxa"))
fwrite(temp, file = "STOOLS_significant_otus.csv", col.names = TRUE, sep = ",")
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
# fname <- file.path(main_path, path_to_results,"alpha_diversity/alpha.txt")
fname <- file.path("~/Results/STOOLS/results/alpha.csv")
alpha <- read.table(fname, header = TRUE)
alpha$sample <- row.names(alpha)
alpha$observed_species <- NULL

mAlpha <- reshape2::melt(alpha, id.vars = "sample", variable.name = "metric", value.name = "value")

mAlpha$treatment <- metadata$treatment[match(mAlpha$sample,metadata$sample)]

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
# matrice= read.table(file.path(project_folder,"STOOLS/qiime_1.9/results/beta_diversity/weighted_unifrac_CSS_normalized_otu_table.txt"), row.names=1, header=T)
# fname = file.path(main_path, path_to_results, "beta_diversity", "weighted_unifrac_CSS_normalized_otu_table.txt")
fname = file.path("~/Results/STOOLS/results/bray_curtis_distances.csv")
matrice= fread(fname)
matrice <- matrice %>% remove_rownames %>% column_to_rownames(var = "row")
matrice <- matrice[,-1]
fname = file.path(basedir, prjdir, "mapping_milkqua_stools.csv")
metadata <- fread(fname)
names(metadata)[1] <- "sample"
names(metadata)[9] <- "cow"
metadata$samples <- "sample-"
metadata$sample <- paste0(metadata$samples,"",metadata$sample)

samples = filter(metadata) %>% pull(sample)
vec <- rownames(matrice) %in% samples
matrice = matrice[vec,vec]

matrice$type <- as.character(metadata$type[match(row.names(matrice),metadata$sample)])
matx= data.matrix(select(matrice, -c(type)))

## MDS
mds <- cmdscale(as.dist(matx))
mds <- as.data.frame(mds)
mds$Unit <- metadata$Unit[match(rownames(mds), metadata$sample)]
mds$type <- metadata$type[match(rownames(mds), metadata$sample)]
mds <- mutate(mds, type = as.factor(type))

p <- ggplot(mds, aes(V1,V2)) + geom_point(aes(colour = type, shape = Unit), size = 3) 
p <- p + xlab("dim1") + ylab("dim2")
p <- p + stat_ellipse(aes(V1, V2, color = type), type = "norm")
p

fname = file.path(main_path, "figures", "beta_diversity.png")
ggsave(filename = fname, plot = p, device = "png", dpi = 300, width = 6, height = 5)
