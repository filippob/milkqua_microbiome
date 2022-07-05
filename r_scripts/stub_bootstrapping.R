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
#alpha = rename(alpha, id = `sample-id`)
metadata$timepoint[metadata$timepoint == "before_oil" ] <- "T0" 
metadata$timepoint[metadata$timepoint == "after_oil" ] <- "T1" 
metadata$timepoint[metadata$timepoint == "8" ] <- "T2"

boot_sample = function(data,index) {
  
  n = nrow(data)
  vec = sample(1:n, n, replace = TRUE)
  temp = data[vec,]
  vec = c("sample-id", index)
  temp = dplyr::select(temp, all_of(vec))
  
  return(temp)
}

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

res = data.frame("index"=NULL, "stat"=NULL, "pval"=NULL, "coef"=NULL)

for (i in 1:100) {
  
  print(paste("bootstrap replicate n.", i))
  temp = boot_sample(alpha, "ACE")
  temp$treatment <- metadata$treatment[match(temp$`sample-id`,metadata$`sample-id`)]
  temp$timepoint <- metadata$timepoint[match(temp$`sample-id`,metadata$`sample-id`)]
  # temp = baseline_correction(temp, "ACE")
  
  temp=as.data.frame(temp)
  names(temp)[names(temp) == "ACE" ] <- "value" 
  temp <- temp %>%
    add_column(metric = "ACE")
  
  bl_medie <- temp %>%
    filter(timepoint=="T0") %>%
    group_by(metric, treatment) %>%
    summarize(media_bl=mean(value))

  bl_counts <- temp %>%
    group_by(metric, treatment) %>%
    filter(timepoint=="T0") %>%
    arrange(metric) %>%
    rename(value.bl = value)

  # M <- merge(temp,bl_counts, by= "sample-id")
  # names(M)[names(M) == "treatment.x"] <- "treatment" 
  # names(M)[names(M) == "timepoint.x"] <- "timepoint"  
  # M$timepoint.y <- NULL
  # M$treatment.y <- NULL
  # M <- merge(M, bl_medie, by="treatment") 
  
  M <- merge(temp,bl_counts[,c(1,2)],by="sample-id", all.x = T)
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

  
  temp = make_comparison(temp, "ACE")
  res = bind_rows(res, temp)
}

barplot(height = res$pval)

