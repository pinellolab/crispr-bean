import pandas as pd
from perturb_tools._qc.qc import get_outlier_guides


def get_outlier_guides_and_mask(
    bdata,
    condit_col: str,
    replicate_col: str = "replicate",
    mad_z_thres: float = 5,
    abs_RPM_thres: float = 10000,
):
    """Obtains outlier guides and per-(replicate, guide) mask to mask the outlier guides.
    For each experimental condition in `samples`, find outlier guides that shows extreme counts compared to guides.
    For each condition in `condit_col`, guides both with median absolute deviation z > `mad_z_thres` AND the RPM larger than `abs_RPM_thres` will be defined as outlier.

    Args
    --
    bdata: ReporterScreen object
    condit_col: column in bdata.samples with experimental samples
    replicate_col: column in bdata.samples with experimental replicate
    mad_z_thres: median absolute deviation threshold used to define outlier guides.
    abs_RPM_thres: RPM threshold value that will be used to define outlier guides.
    """
    outlier_guides = get_outlier_guides(bdata, condit_col, mad_z_thres, abs_RPM_thres)
    outlier_guides[replicate_col] = bdata.samples.loc[
        outlier_guides["sample"], replicate_col
    ].values
    mask = pd.DataFrame(
        index=bdata.guides.index, columns=bdata.samples[replicate_col].unique()
    ).fillna(1)
    print(outlier_guides)
    for _, row in outlier_guides.iterrows():
        mask.loc[row["name"], row[replicate_col]] = 0
    return outlier_guides, mask


def filter_no_info_target(
    bdata,
    condit_col: str,
    control_condition: str,
    target_col: str = "target",
    write_no_support_targets: bool = False,
    no_support_target_write_path: str = None,
):
    """Remove gRNAs where its target has 0 gRNA counts over all samples and targeting gRNAs"""
    gRNA_counts_df = pd.concat(
        [
            bdata.guides[[target_col]],
            pd.DataFrame(
                bdata.X, index=bdata.guides.index, columns=bdata.samples.index
            ),
        ],
        axis=1,
    )
    zero_support_targets = gRNA_counts_df.groupby(target_col).sum().sum(axis=1) == 0
    zero_support_targets = zero_support_targets[zero_support_targets].index
    if write_no_support_targets:
        zero_support_targets.to_series().to_csv(
            no_support_target_write_path, index=False
        )
    bdata = bdata[~bdata.guides[target_col].isin(zero_support_targets), :].copy()
    return len(zero_support_targets), bdata
