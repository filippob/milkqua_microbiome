library("plyr")
library("vegan")
library("dplyr")
library("tidyr")
library("broom")
library("ggpubr")
library("ggplot2")
library("reshape2")
library("data.table")

metadata <- fread("mapping_file.csv")

##########################
## FIGURE 1 (pie-chart)
#########################
## pie chart of overall data and by timepoint
otu_l2 <- fread("results/taxa_summary/otu_table_filtered_L2.txt", header=TRUE, skip=1)
names(otu_l2)[1] <- "OTU"

mO <- melt(otu_l2,id.vars = c("OTU"), value.name = "counts", variable.name = "sample")

mO$anatomic_portion <- metadata$anatomic_portion[match(mO$sample,metadata$`#SampleID`)]
names(mO)[1] <- "phylum"
mO$phylum <- gsub("^.*;","",mO$phylum)
mO$phylum <- gsub("^p__","",mO$phylum)

dd <- mO %>%
  dplyr::group_by(phylum) %>%
  dplyr::summarise(med_value = median(counts)) %>%
  arrange(desc(med_value))

mO$phylum <- factor(mO$phylum, levels = dd$phylum)
mO$anatomic_portion <- c("colon","caecum","stomach","ileum","duodenum","jejunum")[match(mO$anatomic_portion,unique(mO$anatomic_portion))]
mO$anatomic_portion <- factor(mO$anatomic_portion, levels = c("stomach","duodenum","jejunum","ileum","caecum","colon"))

p <- ggplot(mO, aes(x=phylum, y = counts/14)) + geom_bar(aes(fill=phylum), stat = "identity")
p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 1))
p <- p + facet_wrap(~anatomic_portion)
p <- p + ylab("relative abundance")
p

ggsave(filename = "figures/Figure2.png", plot = p, device = "png")

dd <- mO %>%
  mutate(rel_abund = counts) %>%
  group_by(anatomic_portion,phylum) %>%
  dplyr::summarise(avg_abund = mean(rel_abund)) %>%
  spread(key = anatomic_portion, value = avg_abund)

names(dd)[2:ncol(dd)] <- paste(names(dd)[-1],"avg",sep="_")

ds <- mO %>%
  mutate(rel_abund = counts) %>%
  group_by(anatomic_portion,phylum) %>%
  dplyr::summarise(std_abund = sd(rel_abund)) %>%
  spread(key = anatomic_portion, value = std_abund)

names(ds)[2:ncol(ds)] <- paste(names(ds)[-1],"std",sep="_")

D <- mO %>%
  group_by(phylum) %>%
  do(tidy(anova(lm(counts ~ anatomic_portion, data = .)))) %>%
  filter(term == "anatomic_portion")

dd <- dd %>% inner_join(ds, by = "phylum")

dd$p_value <- D$p.value[match(dd$phylum,D$phylum)]
fwrite(x = dd, file = "tables/phylum.csv", sep = ",", col.names = TRUE)

###############################################################
## FIGURE 3 - bubble chart
###############################################################
otu_l6 <- fread("results/taxa_summary_abs/mapping_file_L6.txt", header=TRUE)
names(otu_l6)[1] <- "sample"
otu_l6 <- select(otu_l6, -c(campione,rabbit,sex,n_sample))
mO <- melt(otu_l6,id.vars = c("sample","anatomic_portion"), value.name = "counts", variable.name = "OTU")

sapply(strsplit(x = as.character(mO$OTU), split = ";"), length)

mO <- mO %>%
  mutate(len = sapply(strsplit(x = as.character(OTU), split = ";"), length),
         taxon = c("phylum","class","order","family","genus")[match(len,seq(2,6))])

mO$anatomic_portion <- c("stomach","duodenum","jejunum","ileum","colon","caecum")[match(mO$anatomic_portion,unique(mO$anatomic_portion))]
mO$anatomic_portion <- factor(mO$anatomic_portion, levels = c("stomach","duodenum","jejunum","ileum","caecum","colon"))

m1 <- mO
m1$OTU <- gsub("^.*;","",m1$OTU)

DX <- m1 %>%
  filter(taxon %in% c("genus","family"),OTU!="") %>%
  group_by(taxon,anatomic_portion) %>%
  mutate(tot=sum(counts)) %>%
  group_by(taxon,anatomic_portion,OTU) %>%
  dplyr::summarise(rel_abund = sum(counts)/mean(tot)) %>%
  filter(rel_abund > 10^-2) %>%
  arrange(anatomic_portion, desc(rel_abund))

DX$taxon <- factor(DX$taxon, levels=rev(c("family","genus")))

taxa <- DX %>%
  group_by(taxon,OTU) %>%
  dplyr::summarize("s"=sum(rel_abund)) %>%
  arrange(desc(s)) %>%
  select(taxon,OTU,s)

DX$OTU <- factor(DX$OTU, levels = rev(taxa$OTU[order(taxa$s)]))
DX <- DX %>%
  arrange(OTU)

p <- ggplot(DX, aes(x = anatomic_portion, y = OTU))
p <- p + geom_point(aes(size = rel_abund, colour = anatomic_portion), alpha = 0.4)
p <- p + facet_wrap(~taxon, scales = "free")
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
               axis.text.y = element_text(size=6),
               axis.title.y = element_blank())
p <- p + labs(size = "abundancy") + guides(colour=FALSE)
p

ggsave(filename = "figures/Figure3.png", plot = p, device = "png")
ggsave(filename = "figures/bubble_chart.pdf", plot = p, device = "pdf")

xx <- DX %>%
  spread(key = anatomic_portion, value = rel_abund)

### figure 3 bis
mO$taxa <- gsub("^.*;","",mO$OTU)
xx <- mO %>%
  separate(col = OTU, into = c("domain","phylum","class","order","family","genus"), sep = ";")

library("tidytext")

DX <- xx %>%
  filter(taxon %in% c("genus","family"),taxa!="") %>%
  select(-c(domain,class,order,family,genus,len)) %>%
  group_by(taxon,anatomic_portion) %>%
  mutate(tot=sum(counts,na.rm = TRUE)) %>%
  group_by(taxon,phylum,anatomic_portion,taxa) %>%
  dplyr::summarise(rel_abund = sum(counts,na.rm = TRUE)/mean(tot, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(rel_abund > 10^-2) %>%
  mutate(phylum = as.factor(phylum),
         taxon_ord = ifelse(taxon == "family",2,1),
         taxa = reorder_within(x = taxa, by = taxon_ord, phylum))

p <- ggplot(DX, aes(x = anatomic_portion, y = taxa))
p <- p + geom_point(aes(size = rel_abund, colour = anatomic_portion), alpha = 0.4)
p <- p + facet_grid(phylum~taxon, scales = "free", space = "free")
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
               axis.text.y = element_text(size=6),
               axis.title.y = element_blank(),
               strip.text.y.right = element_text(angle = 0),
               strip.background.y = element_rect(colour = "black", fill = "white", size = 0, linetype = "blank"),
               panel.border = element_blank(),
               panel.spacing.y = unit(1, "lines"),
               legend.title = element_text(size = 7), 
               legend.text = element_text(size = 7))
p <- p + scale_y_discrete(labels = function(x) sub("_.*","",x))
p <- p + labs(size = "abundancy") + guides(colour=FALSE)
p

ggsave(filename = "figures/Figure3_bis.png", plot = p, device = "png", width = 5.5, height = 8)
ggsave(filename = "figures/Figure3_bis.jpg", plot = p, device = "jpg", width = 5.5, height = 8)
ggsave(filename = "figures/bubble_chart_3.pdf", plot = p, device = "pdf", width = 6.5, height = 10)

# dd <- DX %>%
#   group_by(taxon, anatomic_portion) %>%
#   top_n(5, rel_abund) %>%
#   arrange(anatomic_portion,taxon,desc(rel_abund)) %>%
#   filter(taxon!="family")

dd <- DX %>%
  spread(key = "anatomic_portion", value = "rel_abund") %>%
  arrange(taxon,phylum,taxa)

dd[is.na(dd)] <- 0

mm <- mO %>%
  select(-c(anatomic_portion,len)) %>%
  spread(key = "sample", value = "counts")

mm <- mm %>%
  gather(key = "sample", value = "counts", -c(OTU,taxon,taxa))

mm$anatomic_portion <- metadata$anatomic_portion[match(mm$sample,metadata$`#SampleID`)]
mm$anatomic_portion <- c("stomach","duodenum","jejunum","ileum","colon","caecum")[match(mm$anatomic_portion,unique(mm$anatomic_portion))]
mm$anatomic_portion <- factor(mm$anatomic_portion, levels = c("stomach","duodenum","jejunum","ileum","caecum","colon"))

D <- mm %>%
  group_by(taxon,taxa) %>%
  do(tidy(anova(lm(counts ~ anatomic_portion, data = .)))) %>%
  filter(term == "anatomic_portion")

dd$p_value <- D$p.value[match(dd$taxa,D$taxa)]
fwrite(x = dd, file = "tables/families_and_genera.csv", sep = ",", col.names = TRUE)

###########################################
## FIGURE 3 - CORE MICROBIOME
###########################################
corem_all <- fread("tables/Table_core.csv")

xx <- corem_all %>%
  group_by(anatomic_portion) %>%
  mutate(tot = sum(`avg counts`), relative_cnts = `avg counts`/tot,
         taxa = ifelse(genus == "uncultured", paste(family,genus,sep="_"), genus))

xx <- xx %>%
  group_by(anatomic_portion) %>%
  mutate(taxa = ifelse(relative_cnts > 0.01, taxa, "other"))

palette1 <- c("#E65644","#78ACDD","#7BCFDF","#7579D8","#D286DD","#D3E9DA","#74E175","#D847DA","#7DA786","#7B46E6","#8F6D90",
              "#CDED49","#D9B0E0","#D89D59","#D9CA53","#E2D1A3","#7CE1A2","#6CE3D1","#DA5B9E","#6FE740","#C8E48F","#D7CBD9",
              "#DD9290")

palette1 <- c("#E65644","dimgray","#7BCFDF","#7579D8","#D286DD","goldenrod","#74E175","lightpink","#7DA786","#7B46E6","#8F6D90",
              "#CDED49","#D9B0E0","#D89D59","#D9CA53","#E2D1A3","#7CE1A2","#6CE3D1","#DA5B9E","#6FE740","#C8E48F","#D7CBD9",
              "#DD9290")

xx <- xx %>% group_by(anatomic_portion,taxa) %>% summarise(rel_cnts = sum(relative_cnts))
xx <- xx %>%
  arrange(desc(rel_cnts))

xx$anatomic_portion <- c("jejunum","duodenum","whole_gut","stomach","ileum","caecum","colon")[match(xx$anatomic_portion,unique(xx$anatomic_portion))]
xx$anatomic_portion <- factor(xx$anatomic_portion, levels = c("stomach","duodenum","jejunum","ileum","caecum","colon","whole_gut"))

labs <- xx %>%
  group_by(taxa) %>%
  summarise(tot=sum(rel_cnts)) %>%
  arrange(desc(tot)) %>%
  filter(taxa != "other")

lls <- factor(labs$taxa, levels = labs$taxa)
xx$taxa <- factor(xx$taxa, levels = c(levels(lls), "other"))

p <- ggplot(xx, aes(x=factor(1), y = rel_cnts, fill=taxa))
p <- p + geom_bar(width=1,stat="identity")
p <- p + facet_wrap(~anatomic_portion, ncol = 2)
p <- p + coord_polar(theta='y')
p <- p + xlab("relative abundances") + ylab("")
p <- p + scale_fill_manual(values=palette1)
p <- p + scale_y_continuous(limits = c(0, 1.00001), 
                            # expand = c(0, 0), 
                            breaks = seq(0, 0.9, by = 0.1))
p <- p + theme(text = element_text(size=9),
               strip.text = element_text(size = 8),
               axis.text.x = element_text(size=7),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank(),
               legend.text=element_text(size=6),
               legend.title=element_text(size=7),
               legend.key.size = unit(0.35, "cm"),
               plot.margin=unit(c(0,0.1,0,0),"mm"))
p <- p + guides(fill=guide_legend(title="Genus"))
p

ggsave(filename = "figures/core_microbiome.png", plot = p, device = "png", width = 6, height = 8)
ggsave(filename = "figures/Figure4.png", plot = p, device = "png", width = 6, height = 8)

##############################
## FIGURE - alpha diversity
##############################

alpha <- read.table("results/alpha_diversity/alpha.txt", header = TRUE)
alpha$sample <- row.names(alpha)
alpha$observed_species <- NULL

mAlpha <- melt(alpha, id.vars = "sample", variable.name = "metric", value.name = "value")

mAlpha$group <- metadata$anatomic_portion[match(mAlpha$sample,metadata$`#SampleID`  )]
mAlpha$sex <- metadata$sex[match(mAlpha$sample,metadata$`#SampleID`)]
mAlpha <- na.omit(mAlpha)
mAlpha$group <- c("colon","caecum","stomach","ileum","duodenum","jejunum")[match(mAlpha$group,unique(mAlpha$group))]
mAlpha$group <- factor(mAlpha$group, levels = c("stomach","duodenum","jejunum","ileum","caecum","colon"))

p <- ggplot(mAlpha, aes(x=group,y=value)) 
p <- p + geom_boxplot(aes(fill=group))
p <- p + facet_wrap(~metric, scales = "free_y")
p <- p + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
# p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p <- p + xlab("anatomy")
print(p)  

ggsave(filename = "figures/alpha.png", plot = p, device = "png")
ggsave(filename = "figures/Figure5.png", plot = p, device = "png")

mAlpha$group <- factor(mAlpha$group, levels = c("caecum","colon","duodenum","jejunum","ileum","stomach"))
D <- mAlpha %>%
  group_by(metric) %>%
  do(tidy(lm(value ~ group, data = .))) %>%
  filter(term != "(Intercept)") %>%
  mutate(term = gsub("group","",term))



p <- ggplot(D, aes(x=term, y=p.value))
p <- p + geom_line(aes(group=metric, colour=metric))
p <- p + xlab("group")
p

ggsave(filename = "figures/alpha_significance.png", plot = p, device = "png")
ggsave(filename = "figures/Figure6.png", plot = p, device = "png")


##############################
## FIGURE - beta diversity
##############################
matrice= read.table("results/beta_diversity/bray_curtis_CSS_normalized_otu_table.txt", row.names=1, header=T)
names(matrice) <- gsub("X","",names(matrice))
# 
matrice$group <- as.character(metadata$anatomic_portion[match(row.names(matrice),metadata$`#SampleID`)])
matrice$group <- c("colon","caecum","stomach","ileum","duodenum","jejunum")[match(matrice$group,unique(matrice$group))]

matx= data.matrix(matrice[,-ncol(matrice)])

udder.mds= metaMDS(matx, k=3) #function metaMDS in Vegan

hull_f <- function(df) {
  
  temp <- data.frame(NULL)
  for (ll in unique(df$group)) {
    
    nn <- df[df$group == ll,][chull(df[df$group == ll, c("NMDS1","NMDS2","NMDS3")]),]
    temp <- rbind.data.frame(temp,nn)
  }
  return(temp)
}

udder.scores <- as.data.frame(scores(udder.mds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
udder.scores$group <- matrice$group #  add the grp variable created earlier

hull.data <- hull_f(udder.scores)

g <- ggplot(data=udder.scores, aes(x=NMDS1,y=NMDS2))
g <- g + coord_equal()
g = g + geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=group,group=group),alpha=0.30)
g = g + geom_point(data=udder.scores,aes(x=NMDS1,y=NMDS2,shape=group,colour=group),size=4)
g = g + theme_bw() + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
# g <- g + ggtitle("Anatomy")
print(g)

recoded_labs <- c("large_intestine", "large_intestine", "stomach", "ileum", "small_intestinse", "small_intestinse") 
matrice$group <- recoded_labs[match(matrice$group,unique(matrice$group))]

udder.scores$group <- matrice$group #  add the grp variable created earlier
hull.data <- hull_f(udder.scores)

g1 <- ggplot(data=udder.scores, aes(x=NMDS1,y=NMDS2))
g1 <- g1 + coord_equal()
g1 = g1 + geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=group,group=group),alpha=0.30)
g1 = g1 + geom_point(data=udder.scores,aes(x=NMDS1,y=NMDS2,shape=group,colour=group),size=4)
g1 = g1 + theme_bw() + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
# g1 <- g1 + ggtitle("Anatomy")
print(g1)

p <- ggarrange(g, g1, 
               # labels = c("A","G"),
               # heights = c(2,2),
               
               legend = "right",
               common.legend = FALSE,
               widths = c(10,11),
               ncol = 2, nrow = 1)
p

ggsave(filename = "figures/beta.png", plot = p, device = "png", width = 10, height = 6)
ggsave(filename = "figures/Figure7.png", plot = p, device = "png", width = 10, height = 6)

#########################################
## FIGURE XX: HEATMAP OF SIGNIFICANT OTUs
#########################################
load("taxonomy_timepoint.RData")
D <- to_save[[1]]
DX <- to_save[[2]]
D0 <- to_save[[3]]

D0 <- spread(D0, key = treatment, value = avg_counts)

D1 <- DX %>%
  inner_join(D0, by = c("level" = "level", "new_taxa" = "new_taxa", "timepoint" = "timepoint")) %>%
  mutate(p.value = -log10(p.value), Control = ifelse(log(Control) < 1,0,log10(Control)), Treated = ifelse(log(Treated) < 1, 0, log10(Treated))) %>%
  gather(key = "treatment", value = "counts", -c(level,new_taxa,timepoint,p.value))

D1$level <- factor(D1$level, levels = c("phylum","class","order","family","genus"))

p <- ggplot(D1, aes(x = treatment, y = new_taxa, group=level))
p <- p + geom_tile(aes(fill = counts), colour = "white")
p <- p + facet_grid(level~timepoint, space="free", scales = "free_y")
p <- p + scale_fill_gradient(low = "white", high = "red")
p <- p + theme(strip.text.y = element_text(size = 5), 
               strip.text.x = element_text(size = 6),
               # axis.text.y = element_text(size = 4),
               axis.text.x = element_text(size = 6),
               axis.title = element_text(size = 6))
p <- p + guides(fill=FALSE) + theme(axis.title.y = element_blank(),
                                    axis.text.y = element_blank(),
                                    axis.ticks.y = element_blank())
p
dd <- filter(D1, treatment == "Control") %>% mutate(variable = "p-value")

q <- ggplot(dd, aes(x = factor(1), y = new_taxa, group=level))
q <- q + geom_tile(aes(fill = p.value), colour = "white")
q <- q + facet_grid(level~variable, space="free", scales = "free_y")
q <- q + scale_fill_gradient(low = "yellow", high = "darkred")
q <- q + theme(strip.text = element_text(size = 4), 
               strip.text.x = element_text(size = 6),
               axis.text.y = element_text(size = 5),
               axis.title = element_text(size = 6))
q <- q + guides(fill=FALSE) + theme(
  # axis.title.x = element_blank(),
  # axis.text.x=element_blank(),
  # axis.ticks.x=element_blank(),
  strip.text.y = element_blank(),
  # axis.text.x = element_blank()
  axis.text.x = element_text(size = 6)
)
q <- q + xlab("") 
q

figure_final <- ggarrange(q, p, widths=c(0.25, 0.75), 
                          labels=c("A", "B"))

ggsave(filename = "figures/heatmap.png", plot = figure_final, device = "png", width = 8, height = 12)
ggsave(filename = "figures/heatmap.pdf", plot = figure_final, device = "pdf", width = 8, height = 12)


##########################
## FIGURE 1S (rarefaction)
#########################

library("ggsci")
library("scales")

#######################################################
## recursive function to do the same loop as above
getIncrementalOtus <- function(tab) {
  
  indx <- list("ind"=which.max(colSums(tab>0)),"max"=max(colSums(tab>0)))
  steps <- indx$max
  
  tab <- tab[tab[,indx$ind]==0,]
  
  if(nrow(tab)>1) {
    
    steps <- c(steps,getIncrementalOtus(tab))
  } else {
    
    res <- ifelse(sum(steps)<nrow(tab),nrow(tab)-sum(steps),0)
    if(res==0) res <- NULL
    steps <- c(steps,res)
  }
  return(steps)
}
#######################################################

## otu table from closed_otupicking (unfiltered)
otu <- fread("results/otu_table/otu_table.csv",header=TRUE)
otu <- otu[,-c(1,ncol(otu)),with=FALSE]

seqs <- getIncrementalOtus(as.data.frame(otu))
dd <- cbind.data.frame("sample"=seq(1,ncol(otu)),"otus"=cumsum(c(seqs,rep(0,ncol(otu)-length(seqs)))))

p <- ggplot(dd,aes(x=factor(sample),y=otus, group=1)) + geom_line(size=0.5)
p <- p + xlab("N. of samples") + ylab("")
p <- p + scale_x_discrete(breaks=seq(1, 96, 2))
p <- p + theme(text = element_text(size=5),
               axis.text.x = element_text(size=3))
p


obsOtu <- read.table("results/alpha_diversity_rarefaction/seqs_rarefaction.csv",sep=",", header = TRUE, na.strings = "n/a")
names(obsOtu) <- c("iteration","id","nOTUs")

obsOtu$anatomic_portion <- metadata$anatomic_portion[match(obsOtu$id,metadata$`#SampleID`)]
obsOtu$anatomic_portion <- c("stomach","duodenum","jejunum","ileum","colon","caecum")[match(obsOtu$anatomic_portion,unique(obsOtu$anatomic_portion))]
unique(obsOtu$anatomic_portion)


p1 <- ggplot(obsOtu,aes(x=iteration,y=nOTUs,group=id)) + geom_line(aes(colour=anatomic_portion),size=0.5)
# p1 <- p1 + theme(legend.position="none") 
p1 <- p1 +  scale_color_tron() + xlab("N. of sequences") + ylab("N. of OTUs")
# p1 <- p1 + scale_color_manual(values = c("darkslateblue","deepskyblue","goldenrod1","firebrick","darkseagreen","lightpink"))
p1 <- p1 + theme(text = element_text(size=5), legend.title = element_text(size=5), legend.text = element_text(size=4),
                 legend.key.width = unit(0.3,"cm"))
p1 <- p1 + scale_x_continuous(label=comma)
p1

pp <- ggarrange(p1, p, ncol = 2, widths = c(3,2))
ggsave(filename = "figures/Figure1S.pdf", 
       plot = pp, device = "pdf", width = 7, height = 3.5)
ggsave(filename = "figures/Figure1S.png", 
       plot = pp, device = "png", width = 7, height = 3.5)

obsOtu %>%
  group_by(anatomic_portion) %>%
  summarise(avg_otus = mean(nOTUs), std = sd(nOTUs))
