library(argparse)
library(fst)
parser = ArgumentParser(description = "Extract tsv files from refphase output")
parser$add_argument("--refphase_rData",type = "character",help = "path to refphase .RData file")
parser$add_argument("--ascat_rds",type = "character",help = "path to ASCAT output")
parser$add_argument("--output_dir",type = "character",help = "output directory")
args = parser$parse_args()

load(args$refphase_rData)
phased_segments = as.data.frame(refphase_results$phased_segs)
snps = as.data.frame(refphase_results$phased_snps)
purity_ploidy = refphase_results$sample_data
write.table(phased_segments, file = paste0(args$output_dir, "/phased_segs.tsv"), sep = "\t", row.names = FALSE)
write.table(snps, file = paste0(args$output_dir, "/phased_snps.tsv"), sep = "\t", row.names = FALSE)
write.table(purity_ploidy, file = paste0(args$output_dir, "/purity_ploidy.tsv"), sep = "\t", row.names = FALSE)

ascat_output = readRDS(args$ascat_rds)
df_list <- list()
for (sample_name in names(ascat_output)) {
    df <- ascat_output[[sample_name]][['segments_raw']]
    df_list[[sample_name]] <- df
}
ascat_segments <- do.call(rbind, df_list)
write.table(ascat_segments, file = paste0(args$output_dir, "/ascat_segments.tsv"), sep = "\t", row.names = FALSE)