library(rabhit)
library(tigger)

# testing - source_with_args('../../../python/to_rabhit.R', c('P22_I10_S1_changeo_genotyped_downsampled_clone-pass_singlerep.tsv', '../ref/igblast/V_gapped.fasta', '../ref/igblast/d.fasta','../ref/igblast/j.fasta', 'P22_I10_S1', 'IGHJF2-G2'))

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 8) {
  stop("Usage: to_rabhit infile outfile_prefix chain haplo_gene a1 a2 v_reference [d_reference]", call.=FALSE)
}

infile = args[1]
outfile = args[2]
chain = args[3]
haplo_gene = args[4]
a1 = args[5]
a2 = args[6]
v_reference_file = args[7]
d_reference_file = args[8]


v_ref = readIgFasta(v_reference_file)
seqs = read.csv(infile, sep='\t')

if (chain == 'IGH' || chain == 'TRB' || chain == 'TRD') {
    d_ref = readIgFasta(d_reference_file)
    haplo_db = createFullHaplotype(seqs, toHap_col=c("v_call","d_call"), hapBy_col="j_call", hapBy=haplo_gene, toHap_GERM=c(v_ref, d_ref))
} else {
    haplo_db = createFullHaplotype(seqs, toHap_col=c("v_call"), hapBy_col="j_call", hapBy=haplo_gene, toHap_GERM=c(v_ref))
}
write.table(haplo_db, paste(outfile, '_', haplo_gene, '-', a1, '_', a2, '_haplotype.tab', sep=''), sep='\t', row.names=F)

#pdf(file=paste(outfile, '_', haplo_gene, '-', a1, '_', a2, '_haplotype.pdf', sep=''))
#plotHaplotype(haplo_db)
#dev.off()
