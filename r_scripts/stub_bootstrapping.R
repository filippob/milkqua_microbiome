library("tidyverse")
library("data.table")

alpha = fread("results/alpha_diversity/alpha.txt")
alpha = rename(alpha, id = V1)

boot_sample = function(data,index) {
  
  n = nrow(data)
  vec = sample(1:n, n, replace = TRUE)
  temp = data[vec,]
  vec = c("id", index)
  temp = dplyr::select(temp, all_of(vec))
  
  return(temp)
}

make_comparison = function(btstr_data,index) {
  
  ## we rely on the assumption that the 2nd column is the alpha index 
  ## and the 3rd column is the experimental group
  names(btstr_data)[2] = "index"
  names(btstr_data)[3] = "group"
  
  fit = lm(index ~ group, data = btstr_data)
  
  ## retrieve results
  coefficient = as.numeric(coef(fit)[2])
  stat = anova(fit)[,"F value"][1]
  pval = anova(fit)[,"Pr(>F)"][1]
  
  temp = data.frame("index"=index, "stat"=stat, "pval"=pval, "coef"=coefficient)
  return(temp)
}


res = data.frame("index"=NULL, "stat"=NULL, "pval"=NULL, "coef"=NULL)

for (i in 1:10) {
  
  print(paste("bootstrap replicate n.", i))
  temp = boot_sample(alpha, "ace")
  temp$group = sample(c("high","low"), nrow(temp), replace = TRUE)
  temp = make_comparison(temp, "ace")
  res = bind_rows(res, temp)
}

