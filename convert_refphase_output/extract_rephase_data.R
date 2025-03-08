library(argparse)
parser = ArgumentParser(description = "Extract tsv files from refphase output")
parser$add_argument("--refphase_rData",type = "character",help = "path to refphase .RData file")
parser$add_argument("--output_dir",type = "character",help = "output directory")
args = parser$parse_args()

load(args$refphase_rData)
segments = as.data.frame(refphase_results$phased_segs)
snps = as.data.frame(refphase_results$phased_snps)
purity_ploidy = refphase_results$sample_data
write.table(segments, file = paste0(args$output_dir, "/phased_segs.tsv"), sep = "\t", row.names = FALSE)
write.table(snps, file = paste0(args$output_dir, "/phased_snps.tsv"), sep = "\t", row.names = FALSE)
write.table(purity_ploidy, file = paste0(args$output_dir, "/purity_ploidy.tsv"), sep = "\t", row.names = FALSE)