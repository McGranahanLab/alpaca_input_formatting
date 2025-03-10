import pandas as pd
import argparse
import os
from functions import calculate_confidence_intervals_logr
# arguments
parser = argparse.ArgumentParser(
    description="Calculate confidence intervals from refphase output"
)
parser.add_argument(
    "--tumour_id", type=str, help="Unique identifier for the tumour", required=True
)
parser.add_argument("--output_dir", type=str, help="Output directory", required=True)
parser.add_argument(
    "--refphase_segments",
    type=str,
    help="Location of refphase segments file",
    required=True,
)
parser.add_argument(
    "--refphase_snps", type=str, help="Location of refphase snps file", required=True
)
parser.add_argument(
    "--refphase_purity_ploidy",
    type=str,
    help="Location of refphase purity ploidy file",
    required=True,
)

# options
parser.add_argument(
    "--heterozygous_SNPs_threshold",
    type=int,
    default=5,
    help="Minimum number of heterozygous SNPs to consider a segment. Segments with fewer heterozygous SNPs will be discarded.",
)
parser.add_argument(
    "--ci_value", type=float, default=0.5, help="Confidence interval value"
)
parser.add_argument(
    "--n_bootstrap", type=int, default=100, help="Number of bootstrap samples"
)

args = parser.parse_args()
tumour_id = args.tumour_id
output_dir = args.output_dir
ci_value = args.ci_value
n_bootstrap = args.n_bootstrap
# create output directory:
output_dir_segments = f"{output_dir}/segments"
os.makedirs(output_dir_segments, exist_ok=True)

# read data
refphase_segments = pd.read_csv(args.refphase_segments, sep="\t")
refphase_snps = pd.read_csv(args.refphase_snps, sep="\t")
refphase_purity_ploidy = pd.read_csv(args.refphase_purity_ploidy, sep="\t")

# load ASCAT output for segments which were not updated by refPhase:


# remove segments with fewer than 5 heterozygous SNPs:
refphase_segments = refphase_segments[
    refphase_segments["heterozygous_SNP_number"] >= args.heterozygous_SNPs_threshold
]

# rename columns:

refphase_segments = refphase_segments.rename(
    columns={
        "group_name": "sample",
        "seqnames": "chr",
        "patient_tumour": "tumour_id",
    }
)

refphase_snps = refphase_snps.rename(
    columns={
        "group_name": "sample",
        "seqnames": "chr",
        "patient_tumour": "tumour_id",
    }
)

# create segment column by combining chromosome, start and end:
refphase_segments["segment"] = (
    refphase_segments["chr"].astype(str)
    + "_"
    + refphase_segments["start"].astype(str)
    + "_"
    + refphase_segments["end"].astype(str)
)


## calculate confidence intervals:
print(f"Calculating confidence intervals for {tumour_id}")
# assign SNPS to segments:
snps_with_segments = refphase_snps.merge(
    refphase_segments[["chr", "start", "end", "segment", "sample", "cn_a", "cn_b"]],
    left_on=["sample", "chr"],
    right_on=["sample", "chr"],
    how="inner",
)
snps_with_segments = snps_with_segments[
    (snps_with_segments["pos"] >= snps_with_segments["start"])
    & (snps_with_segments["pos"] <= snps_with_segments["end"])
]
# add purity and ploidy information
snps_with_segments_purity_ploidy = snps_with_segments.merge(
    refphase_purity_ploidy, left_on="sample",right_on='sample_id', how="inner"
)

# estimate the confidence intervals:
confidence_intervals = (
    snps_with_segments_purity_ploidy.groupby(["segment", "sample"])
    .apply(
        calculate_confidence_intervals_logr, ci_value=ci_value, n_bootstrap=n_bootstrap
    )
    .reset_index()
)
ci_table = confidence_intervals.merge(refphase_segments)[
    [
        "segment",
        "sample",
        "cn_a",
        "cn_b",
        "lower_CI_A",
        "upper_CI_A",
        "lower_CI_B",
        "upper_CI_B",
        "was_cn_updated",
    ]
]
ci_table["tumour_id"] = tumour_id
ci_table["ci_value"] = ci_value

for allele in ["a", "b"]:
    assert all(
        ci_table[f"cn_{allele}"] >= ci_table[f"lower_CI_{allele}"]
    ), f"cn_{allele} >= lower_CI_{allele}"
    assert all(
        ci_table[f"cn_{allele}"] < ci_table[f"upper_CI_{allele}"]
    ), f"cn_{allele} < upper_CI_{allele}"

ci_table.drop(columns=["cn_a", "cn_b", "was_cn_updated"], inplace=True)
ci_table.to_csv(f"{output_dir}/ci_table.csv", index=False)
print(f'{tumour_id} done')

# TODO create input based on CI table instead
# keep only relevant columns:
print(f"Creating ALPACA input table for {tumour_id}")
alpaca_input = refphase_segments[["sample", "chr", "segment", "cn_a", "cn_b"]].copy()
# rename columns:
alpaca_input = alpaca_input.rename(
    columns={
        "cn_a": "cpnA",
        "cn_b": "cpnB",
    }
)

alpaca_input["tumour_id"] = tumour_id
# write to file:

alpaca_input.to_csv(f"{output_dir}/ALPACA_input_table.csv", index=False)

# split input into separate files for each segment to faciliate parallel processing:
for segment in alpaca_input["segment"].unique():
    alpaca_input[alpaca_input["segment"] == segment].to_csv(
        f"{output_dir_segments}/ALPACA_input_table_{tumour_id}_{segment}.csv",
        index=False,
    )
