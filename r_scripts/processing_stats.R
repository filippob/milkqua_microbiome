library("tidyverse")
library("data.table")

project_folder = "/home/filippo/Documents/MILKQUA"
path_to_file = "Analysis/milkqua_skinswab/qiime1.9"
fname1 = "2.join_reads/readsPerSample.tsv" ## n. of reads after joining
fname2 = "3.quality_filtering/reads_after_filter.tsv" ## n. of reads after filtering

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

D


p <- ggplot(reads,aes(x=sample,y=n_reads)) + geom_bar(aes(fill=sample), stat = "identity")
p <- p + guides(fill='none')
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=3))
print(p)

## after filtering

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

seqs <- after_reads %>%
  mutate(loss=(inp_seqs-out_seqs)/inp_seqs) %>%
  arrange(loss)

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
D

d1 <- seqs %>%
  mutate(retained=out_seqs/inp_seqs) %>%
  summarize(maxRetained=max(retained),minRetained=min(retained),avgRetained=mean(retained))

d1

mS <- reshape2::melt(seqs,id.vars = c("sample_id","loss"), variable.name = "seq", value.name = "reads")
mS <- mS %>%
  filter(!grepl("Undetermined",sample_id)) %>%
  arrange(reads)

mS$sample_id <- factor(mS$sample_id, levels = seqs$sample_id)

p <- ggplot(mS, aes(x=sample_id,y=reads,group=seq)) + geom_line(aes(colour=seq))
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=4))
p
