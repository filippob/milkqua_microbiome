library("tidyverse")
library("data.table")

## PARAMETERS
HOME <- Sys.getenv("HOME")
prj_folder = file.path(HOME, "Results")
analysis_folder = "SKINSWABS"
fname = "results/alpha.csv"
conf_file = "mapping_milkqua_skinswabs.csv"
outdir = file.path(analysis_folder)

metadata = fread (file.path(prj_folder, analysis_folder, conf_file))
metadata$`sample-id` <- gsub('-', '.', metadata$`sample-id`)
alpha = fread(file.path(prj_folder, analysis_folder, fname))

## Correction for baseline

mAlpha <- reshape2::melt(alpha, id.vars = "sample-id", variable.name = "metric", value.name = "value")
mAlpha$timepoint <- metadata$timepoint[match(mAlpha$`sample-id`,metadata$`sample-id`)]
mAlpha$treatment <- metadata$treatment[match(mAlpha$`sample-id`,metadata$`sample-id`)]
mAlpha$timepoint[mAlpha$timepoint == "before_oil" ] <- "T0" 
mAlpha$timepoint[mAlpha$timepoint == "after_oil" ] <- "T1" 
mAlpha$timepoint[mAlpha$timepoint == "8" ] <- "T2"

bl_medie <- mAlpha %>%
  filter(timepoint=="T0") %>%
  group_by(metric,treatment) %>%
  summarize(media_bl=mean(value))

bl_counts <- mAlpha %>%
  group_by(metric,treatment) %>%
  filter(timepoint=="T0") %>%
  arrange(metric) %>%
  rename(value.bl = value)

M <- merge(mAlpha,bl_counts[,c(1,2,3)],by=c("metric","sample-id"),all.x = TRUE)
M1 <- M%>%
  filter(!is.na(M$value.bl))

M2 <- M %>%
  filter(is.na(value.bl)) %>%
  mutate(value.bl=replace(value.bl,is.na(value.bl),right_join(bl_medie, ., by =c("metric","treatment"))$media_bl))

M <- rbind.data.frame(M1,M2)  

M <- M %>%
  mutate(corrected_counts=value-value.bl) %>%
  filter(value!=0) %>%
  arrange(`sample-id`)

M <- M %>%
  group_by(metric) %>%
  mutate(scaled_counts=scales::rescale(corrected_counts,c(0,100)))

# Balpha <- M %>% filter(timepoint == "T1") # options are T1 and T2
# Balpha <- Balpha %>% filter (metric == "Shannon") # options are "Observed", "Chao1", "se.chao1","ACE","se.ACE","Shannon","Simpson","InvSimpson","Fisher"
# Balpha <- Balpha %>% spread(metric, corrected_counts)
M <- M[-c(3,6,8)]

alpe <- M %>% spread(metric, corrected_counts)
alpe <- alpe[alpe$timepoint != "T0",]

Indexes <- colnames(alpe[4:12])
Timepoints <- c("T1", "T2")


## bootstrapping

boot_sample = function(data,index) {
  
  n = nrow(data)
  vec = sample(1:n, n, replace = TRUE)
  temp = data[vec,]
  vec = c("sample-id", index)
  temp = dplyr::select(temp, all_of(vec))
  temp$treatment <- mAlpha$treatment[match(temp$`sample-id`,mAlpha$`sample-id`)]
  temp$timepoint <- mAlpha$timepoint[match(temp$`sample-id`,mAlpha$`sample-id`)]
  return(temp)
  
}

## make comparison

make_comparison = function(btstr_data,index) {
  
  ## we rely on the assumption that the 2nd column is the alpha index 
  ## and the 3rd column is the experimental group
  names(btstr_data)[2] = "index"
  names(btstr_data)[3] = "treatment"
  
  fit = lm(index ~ treatment, data = btstr_data)
  
  ## retrieve results
  coefficient = as.numeric(coef(fit)[2])
  stat = anova(fit)[,"F value"][1]
  pval = anova(fit)[,"Pr(>F)"][1]
  
  temp = data.frame("index"=index, "stat"=stat, "pval"=pval, "coef"=coefficient)
  return(temp)
}

res = data.frame("index"=NULL, "stat"=NULL, "pval"=NULL, "coef"=NULL, "timepoint"=NULL)

for (k in Indexes) {

for (i in 1:1000) {
  
  print(paste("bootstrap replicate n.", i))
  temp = boot_sample(alpe, k)
  
  temp1 = temp[temp$timepoint == "T1", ] 
  temp1 = make_comparison(temp1, k)
  temp1['Timepoint'] = 'T1'
  res = bind_rows(res, temp1)
  
  temp2 = temp[temp$timepoint == "T2", ] 
  temp2 = make_comparison(temp2, k)
  temp2['Timepoint'] = 'T2'
  res = bind_rows(res, temp2)

}

write.csv(res, "~/Results/SKINSWABS/results/Bootstrap_res.csv")

}

## create graph and observations

# res = res[order(res$pval, decreasing = TRUE),]
# print(paste("p-value median is", median(res$pval)))
# print(paste("p-value statistic is", median(res$stat)))

# png("barplot.png")
# barplot(height = res$pval)
# dev.off()

