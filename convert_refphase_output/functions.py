import pandas as pd
import numpy as np


def estimate_cn_ascat(baf, logr, purity, ploidy, logr_compaction=1.0):
    cn = (
        purity
        - 1
        + baf * 2 ** (logr / logr_compaction) * ((1 - purity) * 2 + purity * ploidy)
    ) / purity
    return cn


def bootstrap_sample(data, i):
    return data.sample(frac=1, replace=True, random_state=i)


def calculate_final_value_cn_tot(seg_sample_df, logr_shift=0, logr_scale=1):
    baf = 1
    mean_logr = seg_sample_df["logr"].mean()
    purity = seg_sample_df["purity"].unique()[0]
    ploidy = seg_sample_df["ploidy"].unique()[0]
    final_value = estimate_cn_ascat(
        baf, logr_shift + (mean_logr * logr_scale), purity, ploidy
    )
    return final_value


def calculate_confidence_intervals_logr(seg_sample_df, ci_value=0.95, n_bootstrap=1000):
    cn_tot = seg_sample_df["cn_a"] + seg_sample_df["cn_b"]
    assert len(cn_tot.unique()) == 1
    cn_tot = cn_tot.unique()[0]
    bootstrap_values = []
    for i in range(n_bootstrap):
        bootstrap_sample_df = bootstrap_sample(seg_sample_df, i)
        bootstrap_value = calculate_final_value_cn_tot(bootstrap_sample_df)
        bootstrap_values.append(bootstrap_value)
    # remove nans from the bootstrap values:
    bootstrap_values = [x for x in bootstrap_values if not np.isnan(x)]
    lower_bound = np.percentile(bootstrap_values, (1 - ci_value) / 2 * 100)
    upper_bound = np.percentile(bootstrap_values, (1 + ci_value) / 2 * 100)
    ci_span = upper_bound - lower_bound
    # center around the original cn_tot value:
    lower_bound = cn_tot - ci_span / 2
    upper_bound = cn_tot + ci_span / 2
    if cn_tot > 0:
        a_frac = seg_sample_df["cn_a"].values[0] / cn_tot
        b_frac = seg_sample_df["cn_b"].values[0] / cn_tot
    else:
        a_frac, b_frac = 0.5, 0.5
    # apply the ratio to bounds:
    ## prevent negative values
    lower_CI_A = max(lower_bound * a_frac, 0)
    lower_CI_B = max(lower_bound * b_frac, 0)
    ## ci range cannot be 0, this ensures that minimum span will be 0.001
    upper_CI_A = max(upper_bound * a_frac, 0.001)
    upper_CI_B = max(upper_bound * b_frac, 0.001)

    return pd.DataFrame(
        {
            f"lower_CI_A": lower_CI_A,
            f"upper_CI_A": upper_CI_A,
            f"lower_CI_B": lower_CI_B,
            f"upper_CI_B": upper_CI_B,
        },
        index=[0],
    )
