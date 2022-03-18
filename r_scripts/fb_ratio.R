
### script to calculate the F:B ratio from microbiota samples

library("broom")
library("vegan")
library("ggpubr")
library("ggplot2")
library("reshape2")
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

## OTU - phylum
otu <- fread(file.path(main_path, path_to_results, "taxa_summary_abs/CSS_normalized_otu_table_L2.txt"), header = TRUE, skip = 1)
otu$`#OTU ID` <- gsub("^.*;","",otu$"#OTU ID")
otu <- gather(otu, key = "sample", value ="counts", -`#OTU ID`) %>% spread(key = `#OTU ID`, value = counts)
otu$treatment = metadata$treatment[match(otu$sample,metadata$sample)]
otu <- filter(otu, sample %in% meta_subset$sample)

## F:B ratio
mO <- gather(otu, key = "phylum", value = "counts", -c(sample,treatment))
mO <- filter(mO, phylum %in% c("Bacteroidetes", "Firmicutes")) %>% spread(key = "phylum", value = "counts")

D <- mO %>%
  select(sample, treatment, Bacteroidetes, Firmicutes) %>%
  group_by(treatment) %>%
  mutate(ratio=Firmicutes/Bacteroidetes) %>%
  summarize("F/B_avg"=mean(ratio),"B_avg"=mean(Bacteroidetes),"F_avg"=mean(Firmicutes),
            "F/B_med"=median(ratio),"B_med"=median(Bacteroidetes),"F_med"=median(Firmicutes))

temp <- mO %>%
  mutate(ratio=Firmicutes/Bacteroidetes) 

p <- ggplot(temp, aes(x=treatment,y=ratio)) + geom_boxplot(aes(fill=treatment))
p

## significance of differences
mO %>%
  select(sample, treatment, Bacteroidetes, Firmicutes) %>%
  mutate(ratio=Firmicutes/Bacteroidetes) %>%
  do(glance(lm(ratio ~ treatment, .)))

mO %>%
  select(sample, treatment, Bacteroidetes, Firmicutes) %>%
  mutate(ratio=Firmicutes/Bacteroidetes) %>%
  do(tidy(lm(ratio ~ treatment, .)))

### BOOTSTRAPPING

medBoot <- function(x) {
  
  names(x) <- c("sample","group","Bacteroidetes","Firmicutes")
  labels <- unique(x$group)
  
  ind <- sample(nrow(x),nrow(x),replace = TRUE)
  x <- x[ind,]
  x$ratio <- x$Firmicutes/x$Bacteroidetes
  meds <- tapply(x$ratio, x$group, median)
  naam <- labels[!(labels %in% names(meds))]
  
  ## fill in voids (zero-sized groups due to bootstrapping)
  if(length(naam)>0) {
    z <- rep(NA, length(naam))
    names(z) <- naam
    meds <- c(meds,z)
  }
  
  meds <- meds[order(names(meds))]
  return(meds)
}

group_names <- names(medBoot(mO))
# res <- replicate(10,medBoot(dx1),simplify = TRUE)

n <- 10000
res <- replicate(n,medBoot(mO),simplify = FALSE)
res <- matrix(unlist(res), ncol = n, byrow = FALSE)
res <- t(res)
res <- as.data.frame(res)
names(res) <- group_names
save(res,file = "boot_res.RData")

mR <- reshape2::melt(res, variable.name = "treatment", value.name = "ratio")

D <- mR %>%
  group_by(treatment) %>%
  summarise(med=median(ratio,na.rm=TRUE))

fwrite(D, file = "fb_ratio.csv")


p <- ggplot(mR, aes(x=ratio, fill=treatment))
p <- p + geom_density(alpha=0.25)
p

to_save = list(temp,mR)
save(to_save, file = "fb.RData")


