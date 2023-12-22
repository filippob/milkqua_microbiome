homedir="~/Results/PIPELINE_COMPARISON/"
outputfile="~/Results/PIPELINE_COMPARISON/UofGuelph/Mockrobiota/otutab-Vsearch_x_Silva_138-mock12.txt"
taxfile="~/Downloads/taxa_itgdb_taxa_mothur.txt"
output_dir=homedir

library("dplyr")
library("data.table")

allowed_parameters = c(
  'homedir',
  'outputfile',
  'taxfile',
  'output_dir'
)

args <- commandArgs(trailingOnly = TRUE)

print(args)
for (p in args){
  pieces = strsplit(p, '=')[[1]]
  #sanity check for something=somethingElse
  if (length(pieces) != 2){
    stop(paste('badly formatted parameter:', p))
  }
  if (pieces[1] %in% allowed_parameters)  {
    assign(pieces[1], pieces[2])
    next
  }
  
  #if we get here, is an unknown parameter
  stop(paste('bad parameter:', pieces[1]))
}

print(paste("home directory:",homedir))
print(paste("Otu table w/o taxonomy:",outputfile))
print(paste("taxonomy:",taxfile))
print(paste("Otu table with taxonomy is going to be saved  in ",output_dir))

## reading the data
writeLines(" - reading the files")
fname1 = file.path(outputfile)
df1  <- fread(fname1)
fname2 = file.path(taxfile)
df2 <- fread(fname2, header = F, sep = ";", fill = T) #header False for s138, for others T. Silva_138 requires "\t" as sep, the others ";". s138 requires fill = T
# df2$Confidence <- NULL #only for Silva_138

writeLines(" - merging the tables")
colnames(df1)[1] <- "#OTU_ID"
colnames(df2)[1] <- "#OTU_ID"
df2$`#OTU_ID` <- as.character(df2$`#OTU_ID`)
df1$`#OTU_ID` <- as.character(df1$`#OTU_ID`)
data <- left_join(df1,df2, by="#OTU_ID")
write.table(data, file = "~/Results/PIPELINE_COMPARISON/UofGuelph/Mockrobiota/otu_table-Vsearch_x_Silva_138-mock12.csv")
print("R finished")

