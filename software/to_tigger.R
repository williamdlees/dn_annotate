library(tigger)
library(dplyr)
library(stringr)


args = commandArgs(trailingOnly=TRUE)

if (length(args)!=6) {
  stop("Usage: to_tigger infile infer_alleles v_reference d_reference j_reference outfile_prefix", call.=FALSE)
}

infile = args[1]
find_unmutated = (args[2] == "true")
v_reference_file = args[3]
d_reference_file = args[4]
j_reference_file = args[5]
outfile = args[6]

v_ref = readIgFasta(v_reference_file)
d_ref = readIgFasta(d_reference_file)
j_ref = readIgFasta(j_reference_file)
seqs = read.csv(infile, sep='\t')

# get the V genotype
geno = inferGenotypeBayesian(seqs, germline_db=v_ref, find_unmutated=find_unmutated, v_call='v_call')
all_geno = geno

# get the D genotype
geno = inferGenotypeBayesian(seqs, germline_db=d_ref, find_unmutated=F, v_call='d_call')
all_geno = rbind(all_geno, geno)

# get the D genotype
geno = inferGenotypeBayesian(seqs, germline_db=j_ref, find_unmutated=F, v_call='j_call')
all_geno = rbind(all_geno, geno)

GENOTYPED_ALLELES <- function(y) {
  m <- which.max(as.numeric(y[2:5]))

  paste0(unlist(strsplit((y[1]), ','))[1:m], collapse = ",")
}

all_geno$GENOTYPED_ALLELES <-
  apply(all_geno[, c(2, 6:9)], 1, function(y) {
    GENOTYPED_ALLELES(y)
  })

write.table(all_geno, paste(outfile, 'genotype.tsv', sep=''), sep='\t', row.names=F)

# fudge the genotyped alleles into column 2, where the remaining functions expect them

all_geno$alleles = all_geno$GENOTYPED_ALLELES

genotype_db <- genotypeFasta(all_geno[grep('IGHV', substr(all_geno$gene, 1, 4)), ], v_ref, c())
writeFasta(genotype_db, paste(outfile, 'v_germline_gapped.fasta', sep=''))
genotype_db = sapply(genotype_db, function(x) {str_replace_all(x, '\\.', '')})
writeFasta(genotype_db, paste(outfile, 'v_germline.fasta', sep=''))

genotype_db <- genotypeFasta(all_geno[grep('IGHD', substr(all_geno$gene, 1, 4)), ], d_ref, c())
writeFasta(genotype_db, paste(outfile, 'd_germline.fasta', sep=''))

genotype_db <- genotypeFasta(all_geno[grep('IGHJ', substr(all_geno$gene, 1, 4)), ], j_ref, c())
writeFasta(genotype_db, paste(outfile, 'j_germline.fasta', sep=''))

pdf(file=paste(outfile, 'genotype.pdf', sep=''))
plotGenotype(all_geno, text_size=10)
dev.off()


