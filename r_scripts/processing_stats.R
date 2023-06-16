## script to obtain statistics on the reads processed through the bioinformatics pipeline
## this is currently designed on the output from Qiime 1.9

## 1. SET UP
library("tidyverse")
library("data.table")

## 2. PARAMETERS
project_folder = "/home/mycelium/Results/Milkqua_skinswab_bootstrapping"
path_to_file = "results"
fname1 = "readsPerSample.tsv" ## n. of reads after joining
fname2 = "reads_after_filter.tsv" ## n. of reads after filtering
out_fname = "processing_stats.tsv"
sep_line = list("\n--------------------------------------\n")
header = list("sample \t input_seqs \t output_seqs \t loss\n")

## 3. READS AFTER JOINING
writeLines(" - 1) reads after joining")
reads = fread(file.path(project_folder,path_to_file,fname1))
names(reads) <- c("sample","n_reads")
reads$sample <- gsub("/","",reads$sample)

D <- reads %>%
  dplyr::summarize(
    n=n(),
    avgReads=mean(n_reads),
    sdReads = sd(n_reads),
    maxReads=round(max(n_reads),3),
    minReads=round(min(n_reads),3)
  )


fwrite(x = D, file = file.path(project_folder,path_to_file,out_fname), sep = "\t")
fwrite(x = sep_line, file = file.path(project_folder,path_to_file,out_fname), sep = "\t", append = TRUE)
fwrite(x = header, file = file.path(project_folder,path_to_file,out_fname), sep = "\t", append = TRUE)

p <- ggplot(reads,aes(x=sample,y=n_reads)) + geom_bar(aes(fill=sample), stat = "identity")
p <- p + guides(fill='none')
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=3))
ggsave(filename = file.path(project_folder,path_to_file, "figures", "1.reads_per_sample.png"), plot = p, device = "png")

## 3. after filtering
writeLines(" - 2) reads after filtering")
after_reads = fread(file.path(project_folder,path_to_file,fname2), fill = TRUE, sep = ":", header = FALSE)

prefix = "/gpfs/home/projects/MILKQUA/Analysis/milkqua_skinswab/qiime1.9/join_paired_ends/"
suffix = "reads/seqprep_assembled.fastq.gz \\(md5" ## escaping barckets
samples = filter(after_reads, grepl("Sequence read filepath", V1)) %>% select(V2) %>% rename(sample_id = V2)
samples$sample_id <- gsub(prefix,"",samples$sample_id)
samples$sample_id <- gsub(suffix,"",samples$sample_id)

input_seqs = filter(after_reads, grepl("Total number of input sequences", V1)) %>% select(V2) %>% rename(inp_seqs = V2)

output_seqs = filter(after_reads, grepl("Total number seqs written", V1)) %>% select(V1)
output_seqs <- output_seqs %>% separate(col = V1, into = c("V1","out_seqs"), sep = "\t") %>% select(out_seqs)

after_reads = cbind.data.frame(samples, input_seqs, output_seqs)
after_reads = mutate(after_reads, 
       inp_seqs = as.numeric(inp_seqs),
       out_seqs = as.numeric(out_seqs)
       )

writeLines(" - 2) loss due to filtering")
seqs <- after_reads %>%
  mutate(loss=(inp_seqs-out_seqs)/inp_seqs) %>%
  arrange(loss)

fwrite(x = seqs, file = file.path(project_folder,path_to_file,out_fname), sep = "\t", append = TRUE)
fwrite(x = sep_line, file = file.path(project_folder,path_to_file,out_fname), sep = "\t", append = TRUE)

D <- seqs %>%
  # .(group),
  summarize(
    "inputSeq"=sum(inp_seqs),
    "outputSeqs"=sum(out_seqs),
    "avgInp"=mean(inp_seqs),
    "avgOutput"=mean(out_seqs),
    "stdInp"=sd(inp_seqs),
    "stdOutput"=sd(out_seqs,na.rm=TRUE)
  )

header = list("tot input_seqs \t tot output_seqs \t avg input seqs \t avg output_seqs \t std inp_seqs \t std out_seqs\n")
fwrite(x = header, file = file.path(project_folder,path_to_file,out_fname), sep = "\t", append = TRUE)
fwrite(x = D, file = file.path(project_folder,path_to_file,out_fname), sep = "\t", append = TRUE)
fwrite(x = sep_line, file = file.path(project_folder,path_to_file,out_fname), sep = "\t", append = TRUE)


d1 <- seqs %>%
  mutate(retained=out_seqs/inp_seqs) %>%
  summarize(maxRetained=max(retained),minRetained=min(retained),avgRetained=mean(retained))

header = list("max retained \t min retained \t avg retained\n")
fwrite(x = header, file = file.path(project_folder,path_to_file,out_fname), sep = "\t", append = TRUE)
fwrite(x = d1, file = file.path(project_folder,path_to_file,out_fname), sep = "\t", append = TRUE)
fwrite(x = sep_line, file = file.path(project_folder,path_to_file,out_fname), sep = "\t", append = TRUE)


mS <- reshape2::melt(seqs,id.vars = c("sample_id","loss"), variable.name = "seq", value.name = "reads")
mS <- mS %>%
  filter(!grepl("Undetermined",sample_id)) %>%
  arrange(reads)

mS$sample_id <- factor(mS$sample_id, levels = seqs$sample_id)

p <- ggplot(mS, aes(x=sample_id,y=reads,group=seq)) + geom_line(aes(colour=seq))
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=4))

ggsave(filename = file.path(project_folder,path_to_file,, "figures", 2.reads_after_filtering.png"), plot = p, device = "png")

writeLines("DONE!!")


