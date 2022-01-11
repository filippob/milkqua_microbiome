library("tidyverse")

mapping <- readxl::read_excel("milkqua_rumen_samples.xlsx")

## renaming
names(mapping) <- c("nid","exp","vial","treatment","dosis","inoculation","replicate")

## create duplicated id and reorder
mapping$`#sampleID` <- paste(mapping$nid,mapping$nid,sep="")
mapping <- relocate(mapping, `#sampleID`, .before = nid)

## write out mapping file
fwrite(x = mapping, file = "mapping_file.csv", sep = "\t", row.names = FALSE, col.names = TRUE)
