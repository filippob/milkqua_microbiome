library("dplyr")
library("data.table")

#########################
## PARAMETERS   #########
#########################
homedir = "/home/users/chiara.gini/MILKQUA" ## home to where the original mapping file is located (from gdrive()
# mapping_file = "~/Documents/cremonesi/milkqua/mapping_file.csv.xlsx"
mapping_file = "Config/milkqua_stools_swabs.csv"
sample_sheet = "SampleSheet.csv"
nskip = 16 ## n of rows to skip in the sample sheet (header)
outdir = "Config"
type = "faeces"
orig_data_folder = "220225_M04028_0144_000000000-K6CMG"
abs_path = "/home/users/chiara.gini/MILKQUA/data"
n_samples = 64 ## controls the number of samples randomly selected

## reading the data
writeLines(" - reading the mapping file")
# mapping <- readxl::read_xlsx(mapping_file)
fname = file.path(homedir, mapping_file)
mapping <- fread(fname)
fname = file.path(homedir,"data",orig_data_folder,sample_sheet)
samplesheet = fread(fname, skip=nskip)

writeLines(" - data preprocessing")
## renaming columns
mapping$`sample-id` = paste("sample",mapping$`#SampleID`, sep="-")
mapping = rename(mapping, sample_n = `#SampleID`) %>% relocate(`sample-id`, .before = sample_n)

## sample sheet
samplesheet = mutate(samplesheet, `sample-id` = paste("sample", Sample_ID, sep="-")) %>%
	select(`sample-id`, index, index2)

## extracting file names
mapping$sample_name_R1 = gsub("^/.*/","",mapping$absolute_path_forward)
mapping$sample_name_R2 = gsub("^/.*/","",mapping$absolute_path_reverse)

## join mapping file and sample sheet
mapping <- inner_join(mapping, samplesheet, by="sample-id")

writeLines(" - filtering samples")
print(paste("selecting", type, "samples"))
mapping <- filter(mapping, sample_type == type)

## create manifest file
writeLines(" - making the manifest file")
# ext="fastq.gz"
# int=paste("S",mapping$Sample_ID, sep="")

# forward_filenames = c(
#   file.path(abs_path, orig_data_folder, mapping$sample_name_R1)
# )
# reverse_filenames = c(
#   file.path(abs_path, orig_data_folder, mapping$sample_name_R2)
# )

manifest <- mutate(mapping, 
       forward_absolute_filepath = file.path(abs_path, orig_data_folder, mapping$sample_name_R1), 
       reverse_absolute_filepath = file.path(abs_path, orig_data_folder, mapping$sample_name_R2)
       ) %>% 
  select(`sample-id`, forward_absolute_filepath, reverse_absolute_filepath, index, index2, project, sample_type, samplename, treatment, timepoint)


manifest <- manifest %>%
  arrange(`sample-id`) %>%
#  slice_sample(n = n_samples) %>% ## !! PAY ATTENTION: we are taking a sample of the data !!
  rename(`forward-absolute-filepath` = forward_absolute_filepath,
         `reverse-absolute-filepath` = reverse_absolute_filepath)

writeLines(" - writing out the manifest file")
fname = file.path(outdir, paste("manifest_",type,"2.csv", sep=""))
print(paste("writing to file", fname))
fwrite(manifest, file = fname, sep="\t", col.names = TRUE)

print("R finished")

