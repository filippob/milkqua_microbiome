otu <- fread("/home/mycelium/Results/PIPELINE_COMPARISON/temp-mockprodandataset/Micca/otutable-Micca_x_16SITGDB.97-mockprodandataset.txt", header = TRUE, sep = ",")
otu <- filter(otu, Family !="Mitochondria")
otu <- filter(otu, Class !="Chloroplast")
otu <- filter(otu, Order !="Chloroplast")
otu <- otu %>% select(-c(`sample-2`, `sample-5`, `sample-8`, `sample-12`, `sample-15`, `sample-17`, `sample-19`, `sample-21`, `sample-24`, `sample-25`, `sample-26`, `sample-28`, `sample-31`, `sample-30`, `sample-32`))
otu <- select(otu, 2:33)
otu$sum <- rowSums(otu)
0.0 %in% otu$sum

library(stringr)

otu <- read.csv("~/Results/PIPELINE_COMPARISON/temp-mockprodandataset/Micca/otutable-Micca_x_16SITGDB.97-mockprodandataset.txt", sep = "\t", header = T)
taxa <- read.csv("~/Results/PIPELINE_COMPARISON/temp-mockprodandataset/Micca/taxa_Micca_x_16SITGDB.97-mockprodandataset.txt", sep = "\t", header = F)
colnames(taxa)[1] ="OTU"
colnames(taxa)[2] = "taxonomy"
taxonomy <- merge(otu, taxa, by=c("OTU"),all.x=TRUE)
taxonomy[c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')] <- str_split_fixed(taxonomy$taxonomy, ';', 7)
taxonomy$taxonomy <- NULL
