library("tidyverse")
library("data.table")

#########################
## PARAMETERS   #########
#########################
mapping_file = "~/Documents/cremonesi/milkqua/mapping_file.csv.xlsx"
type = "milk"
abs_path = "/home/users/filippo.biscarini.est/MILKQUA/data/210902_M04028_0139_000000000-JRGYP_milk"
n_samples = 10 ## controls the number of samples randomly selected
outdir = "Config/"

## reading the data
writeLines(" - reading the mapping file")
mapping <- readxl::read_xlsx(mapping_file)
mapping$`sample-id` = paste("sample",mapping$Sample_ID, sep="-")
mapping = rename(mapping, sample_2n = `#SampleID`) %>% relocate(`sample-id`, .before = sample_2n)

writeLines(" - filtering samples")
print(paste("selecting", type, "samples"))
mapping <- filter(mapping, sample_type == type)

## create manifest file
writeLines(" - making the manifest file")
ext="fastq.gz"
int=paste("S",mapping$Sample_ID, sep="")

# filenames = c(
#   paste(mapping$Sample_ID,int,"L001","R1","001",sep="_"),
#   paste(mapping$Sample_ID,int,"L001","R2","001",sep="_")
# )

forward_filenames = c(
  paste(mapping$Sample_ID,int,"L001","R1","001",sep="_")
)
reverse_filenames = c(
  paste(mapping$Sample_ID,int,"L001","R2","001",sep="_")
)

manifest <- data.frame(
  "sample_id" = paste("sample",gsub("_.*$","", forward_filenames), sep="-"),
  "forward_absolute_filepath"=paste(forward_filenames,ext,sep="."),
  "reverse_absolute_filepath"=paste(reverse_filenames,ext,sep="."))


manifest <- manifest %>%
  mutate(# direction = ifelse(grepl("R1",manifest$absolute_filepath),"forward","reverse"),
    forward_absolute_filepath = paste(abs_path,forward_absolute_filepath,sep="/"),
    reverse_absolute_filepath = paste(abs_path,reverse_absolute_filepath,sep="/")) %>%
  arrange(sample_id) %>%
  slice_sample(n = n_samples) %>% ## !! PAY ATTENTION: we are taking a sample of the data !!
  rename(`sample-id` = sample_id,
         `forward-absolute-filepath` = forward_absolute_filepath,
         `reverse-absolute-filepath` = reverse_absolute_filepath)

writeLines(" - writing out the manifest file")
fname = paste(outdir, "manifest_",type,".csv", sep="")
print(paste("writing to file", fname))
fwrite(manifest, file = fname, sep="\t", col.names = TRUE)

print("DONE!!")
