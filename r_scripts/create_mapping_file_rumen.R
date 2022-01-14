library("tidyverse")
library("data.table")

mapping <- readxl::read_excel("mapping_file_rumen.xlsx")

## renaming
mapping = rename(mapping, `#SampleID` = `#sampleID`)

## create duplicated id and reorder
# mapping$`#sampleID` <- paste(mapping$nid,mapping$nid,sep="")
# mapping <- relocate(mapping, `#sampleID`, .before = nid)

## write out mapping file
fwrite(x = mapping, file = "mapping_file.csv", sep = "\t", row.names = FALSE, col.names = TRUE)
