#!/usr/bin/R
# Downsample to a max number of seqs per V germline
# Taken from the VDJbase pipeline (Ayelet Peres) - https://bitbucket.org/yaarilab/processpipeline/src/master/pipelines/vdjbase-pipeline.R

library(tigger)
library(alakazam)


# testing - source_with_args('../../../python/downsample.R', c('../tigger/P22_I10_S1_changeo_genotyped.tab', 'P22_I10_S1_changeo_genotyped_downsampled.tab', '../ref/igblast/v_gapped.fasta','IGH', '500'))

args = commandArgs(trailingOnly=TRUE)
#args = c('../tigger/P22_I10_S1_changeo_genotyped.tab', 'P22_I10_S1_changeo_genotyped_downsampled.tab', '../ref/igblast/v_gapped.fasta','IGH', '500')

if (length(args)!=5) {
  stop("Usage: downsample infile outfile v_reference chain sample_size", call.=FALSE)
}

infile = args[1]
outfile = args[2]
v_ref_file = args[3]
chain = args[4]
sample_size = as.integer(args[5])

VGERM = readIgFasta(v_ref_file)
DATA = read.csv(infile, sep='\t')

### make sure that all sequences have V, D, and J call
calls <-
  unlist(strsplit(
    ifelse(
      chain %in% c("IGL", "IGK"),
      "v_call,j_call",
      "v_call,j_call,d_call"
    ),
    ","
  ))

### read collapsed file and filter
for (call in calls) {
  DATA <- DATA[grepl(toupper(substr(call, 1, 1)), DATA[[call]]),]
}

### filter productive sequences
DATA <- DATA[DATA$productive %in% c(TRUE, "T", "TRUE"),]

### sample X sequences per V gene.
DATA$v_gene <-
  getGene(
    DATA$v_call,
    first = F,
    collapse = T,
    strip_d = T
  )

genes <- unique(getGene(names(VGERM), strip_d = T))

ids <- unlist(sapply(genes, function(g) {
  id <- grep(g, DATA$v_gene, fixed = T)

    id_l <-
    if (sample_size <= 1)
      length(id) * sample_size
  else if (sample_size > length(id))
    length(id)
  else
    sample_size
  sample(id, id_l)
}))

DATA <- DATA[ids, ]

write.table(DATA, outfile, row.names = F, quote = F, sep = '\t')

