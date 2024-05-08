import sys
import abc
import logging
from dataclasses import dataclass
from typing import Optional, Dict, Tuple, List
from xmlrpc.client import Boolean
from copy import deepcopy
import torch
import numpy as np
import pandas as pd
import bean as be
from .get_alpha0 import get_fitted_alpha0, get_pred_alpha0
from .get_pi_alpha0 import get_fitted_alpha0 as get_fitted_pi_alpha0
from .get_pi_alpha0 import get_pred_alpha0 as get_pred_pi_alpha0
from .utils import (
    get_accessibility_guides,
    get_edit_to_index_dict,
    _assign_rep_ids_and_sort,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)-5s @ %(asctime)s:\n\t %(message)s \n",
    datefmt="%a, %d %b %Y %H:%M:%S",
    stream=sys.stderr,
    filemode="w",
)
error = logging.critical
warn = logging.warning
debug = logging.debug
info = logging.info


@dataclass
class ScreenData(abc.ABC):
    def __init__(
        self,
        screen: be.ReporterScreen,
        repguide_mask: str = None,
        sample_mask_column: str = None,
        shrink_alpha: bool = False,
        condition_column: str = "sort",
        sample_covariate_column: List[str] = [],
        control_condition: str = "bulk",
        accessibility_col: str = None,
        accessibility_bw_path: str = None,
        device: str = None,
        replicate_column: str = "replicate",
        popt: Optional[Tuple[float]] = None,
        pi_popt: Optional[Tuple[float]] = None,
        control_can_be_selected: bool = False,
        **kwargs,
    ):
        """
        Args
        condition_column: By default, a single condition column, but you can optionally inlcude sample covariate column
        control_can_be_selected: If True, screen.samples[condition_column] == control_condition can also be included in effect size inference if its condition column is not NA (Currently only suppoted for prolifertion screens).
        """
        # TODO: remove replicate with too small number of (ex. only 1) sorting bin
        self.condition_column = condition_column
        self.device = device
        screen.samples["size_factor"] = self.get_size_factor(screen.X)
        if not (
            replicate_column in screen.samples.columns
            and condition_column in screen.samples.columns
        ):
            screen.samples[replicate_column], screen.samples[condition_column] = zip(
                *screen.samples.index.map(lambda s: s.rsplit("_", 1))
            )
        if condition_column not in screen.samples.columns:
            screen.samples[condition_column] = screen.samples["index"].map(
                lambda s: s.split("_")[-1]
            )
        if "sample_covariates" in screen.uns:
            self.sample_covariates = screen.uns["sample_covariates"]
            self.n_sample_covariates = len(self.sample_covariates)
            screen.samples["_rc"] = screen.samples[
                [replicate_column] + self.sample_covariates
            ].values.tolist()
            screen.samples["_rc"] = screen.samples["_rc"].map(
                lambda slist: ".".join(slist)
            )
            self.rep_by_cov = torch.as_tensor(
                (
                    screen.samples[["_rc"] + self.sample_covariates]
                    .drop_duplicates()
                    .set_index("_rc")
                    .values.astype(int)
                )
            )
            replicate_column = "_rc"
        self.screen = screen
        if not control_can_be_selected:
            self.screen_selected = screen[
                :, screen.samples[condition_column] != control_condition
            ]
        else:
            self.screen_selected = screen[:, ~screen.samples[condition_column].isnull()]

        self.n_condits = len(
            self.screen_selected.var[condition_column].unique()
        )  # excluding bulk
        self.screen_control = screen[
            :, screen.samples[condition_column] == control_condition
        ]
        self.n_samples = len(screen.samples)  # 8
        self.n_guides = len(screen.guides)
        self.n_reps = len(screen.samples[replicate_column].unique())
        self.accessibility_col = accessibility_col
        self.accessibility_bw_path = accessibility_bw_path
        self.replicate_column = replicate_column
        self.condition_column = condition_column
        self.control_condition = control_condition
        self.sample_mask_column = sample_mask_column
        self.repguide_mask = repguide_mask
        self.shrink_alpha = shrink_alpha
        self.popt = popt

    def _post_init(self):
        # Assign accessibility info
        if self.accessibility_col is not None:
            self.guide_accessibility = torch.as_tensor(
                self.screen.guides[self.accessibility_col].values
            )
        elif self.accessibility_bw_path is not None:
            self.guide_accessibility = get_accessibility_guides(
                self.accessibility_bw_path, self.screen.guides
            )
        else:
            self.guide_accessibility = None

        # Mask counts
        if self.sample_mask_column is not None:
            self.sample_mask = torch.as_tensor(
                self.screen_selected.samples[self.sample_mask_column].to_numpy()
            ).reshape(self.n_reps, self.n_condits)
            self.bulk_sample_mask = torch.as_tensor(
                self.screen_control.samples[self.sample_mask_column].to_numpy()
            )
        else:
            self.sample_mask = torch.ones((self.n_reps, self.n_condits), dtype=Boolean)
            self.bulk_sample_mask = torch.ones(self.n_reps, dtype=Boolean)
        self.X = self.transform_data(
            self.screen_selected.X
        )  # (n_reps, n_bins, n_guides)
        self.X_masked = self.X * self.sample_mask[:, :, None]
        self.X_control = self.transform_data(self.screen_control.X, 1)
        self.X_control_masked = self.X_control * self.bulk_sample_mask[:, None, None]
        if self.repguide_mask is None:
            self.repguide_mask = ~(self.X == 0).any(axis=1)
        else:
            info(
                f"Using replicate x guide mask in ReporterScreen.uns['{self.repguide_mask}'] to filter out outlier guides."
            )
            assert (
                self.repguide_mask in self.screen.uns.keys()
            ), f"{self.repguide_mask} not in screen.uns"
            assert self.screen_selected.uns[self.repguide_mask].shape == (
                self.n_guides,
                self.n_reps,
            )
            assert (
                self.screen_selected.uns[self.repguide_mask].index
                == self.screen_selected.guides.index
            ).all()
            assert (
                self.screen_selected.uns[self.repguide_mask].columns
                == self.screen_selected.samples[self.replicate_column].unique()
            ).all()
            self.repguide_mask = (
                torch.as_tensor(self.screen_selected.uns[self.repguide_mask].values.T)
                > 0
            )
            self.repguide_mask = torch.logical_and(
                self.repguide_mask, ~(self.X == 0).any(axis=1)
            )

        self.size_factor = torch.as_tensor(
            self.screen_selected.samples["size_factor"].to_numpy()
        ).reshape(self.n_reps, self.n_condits)
        self.size_factor_control = torch.as_tensor(
            self.screen_control.samples["size_factor"].to_numpy()
        ).reshape(self.n_reps, 1)
        # Get a0
        fitted_a0, self.popt = get_fitted_alpha0(
            self.X.clone().cpu(),
            self.size_factor.clone().cpu(),
            self.sample_mask.cpu(),
            shrink=self.shrink_alpha,
            popt=self.popt,
        )
        fitted_a0 = torch.as_tensor(fitted_a0)
        a0 = fitted_a0
        self.a0 = torch.as_tensor(a0)

    def __getitem__(self, guide_idx):
        ndata = deepcopy(self)
        ndata.screen_selected = ndata.screen_selected[guide_idx, :]
        ndata.screen_control = ndata.screen_control[guide_idx, :]
        ndata.n_guides = len(guide_idx)
        ndata.X = ndata.X[:, :, guide_idx]  # (n_reps, n_bins, n_guides)
        ndata.X_masked = ndata.X_masked[:, :, guide_idx]
        ndata.X_control = ndata.X_control[:, :, guide_idx]
        ndata.repguide_mask = ndata.repguide_mask[:, guide_idx]
        if hasattr(ndata, "a0") and self.a0 is not None:
            ndata.a0 = ndata.a0[guide_idx]
        return ndata

    def transform_data(self, X, n_bins=None):
        if n_bins is None:
            n_bins = self.n_condits
        try:
            x = (
                torch.as_tensor(X)
                .T.reshape((self.n_reps, n_bins, self.n_guides))
                .float()
            )
        except RuntimeError:
            print((self.n_reps, n_bins, self.n_guides))
            print(X.shape)
            exit(1)
        if self.device is not None:
            x = x.cuda()
        return x

    def get_size_factor(self, X: np.array):
        """
        Get size factor for samples.
        """
        n_guides, n_samples = X.shape
        geom_mean_x = np.exp((1 / n_samples) * np.log(X + 0.5).sum(axis=1))
        assert geom_mean_x.shape == (n_guides,)
        norm_count = X / geom_mean_x[:, None]
        size_factor = np.median(norm_count, axis=0)
        if any(size_factor == 0):
            size_factor = np.mean(norm_count, axis=0)
        assert size_factor.shape == (n_samples,)
        return size_factor

    def get_sample_ids_and_sort(self, screen: be.ReporterScreen, condit_id_col: str):
        screen.samples["replicate_id"] = -1
        for i, rep in enumerate(sorted(screen.samples.replicate.unique())):
            screen.var.loc[screen.samples.replicate == rep, "replicate_id"] = i
        if condit_id_col in screen.samples.columns:
            screen = screen[
                :, screen.samples.sort_values(["replicate_id", condit_id_col]).index
            ]
        else:
            raise ValueError("AnnData doesn't have condition id provided.")
        return screen


@dataclass
class ReporterScreenData(ScreenData):
    X_bcmatch: torch.Tensor
    size_factor_bcmatch: torch.Tensor
    X_bcmatch_control: torch.Tensor
    size_factor_bcmatch_control: torch.Tensor
    allele_counts: torch.Tensor
    allele_counts_control: torch.Tensor
    a0_bcmatch: torch.Tensor
    a0_allele: torch.Tensor
    pi_a0: torch.Tensor

    def __init__(
        self,
        screen: be.ReporterScreen,
        repguide_mask: str = None,
        sample_mask_column: str = None,
        pi_popt: Tuple[float] = None,
        impute_pi_popt: bool = False,
        shrink_alpha: bool = False,
        condition_column: str = "sort",
        control_condition: str = "bulk",
        accessibility_bw_path: str = None,
        accessibility_col: str = None,
        use_const_pi: bool = False,
        pi_prior_count: int = 10,
        **kwargs,
    ):
        ScreenData.__init__(
            self,
            screen=screen,
            repguide_mask=repguide_mask,
            sample_mask_column=sample_mask_column,
            shrink_alpha=shrink_alpha,
            condition_column=condition_column,
            control_condition=control_condition,
            pi_popt=pi_popt,
            **kwargs,
        )
        self._post_init(
            screen,
            use_const_pi,
            impute_pi_popt,
            pi_prior_count,
            shrink_alpha,
            pi_popt,
        )

    def _post_init(
        self,
        screen,
        use_const_pi,
        impute_pi_popt,
        pi_prior_count,
        shrink_alpha,
        pi_popt,
    ):
        screen.samples["size_factor_bcmatch"] = self.get_size_factor(
            screen.layers["X_bcmatch"]
        )
        self.screen_selected.samples["size_factor_bcmatch"] = screen.samples.loc[
            self.screen_selected.samples.index, "size_factor_bcmatch"
        ]
        self.screen_control.samples["size_factor_bcmatch"] = screen.samples.loc[
            self.screen_control.samples.index, "size_factor_bcmatch"
        ]
        self.X_bcmatch = self.transform_data(self.screen_selected.layers["X_bcmatch"])
        self.X_bcmatch_masked = self.X_bcmatch * self.sample_mask[:, :, None]
        self.X_bcmatch_control = self.transform_data(
            self.screen_control.layers["X_bcmatch"], 1
        )
        self.X_bcmatch_control_masked = (
            self.X_bcmatch_control * self.bulk_sample_mask[:, None, None]
        )
        self.size_factor_bcmatch = torch.as_tensor(
            self.screen_selected.samples["size_factor_bcmatch"].to_numpy()
        ).reshape(self.n_reps, self.n_condits)
        self.size_factor_bcmatch_control = torch.as_tensor(
            self.screen_control.samples["size_factor_bcmatch"].to_numpy()
        ).reshape(self.n_reps, 1)

        edited_control = self.transform_data(
            self.screen_control.layers["edits"], n_bins=1
        )
        nonedited_control = self.X_bcmatch_control - edited_control
        nonedited_control[nonedited_control < 0] = 0
        if (
            not hasattr(self, "allele_counts_control")
            or self.allele_counts_control is None
        ):
            self.allele_counts_control = torch.stack(
                [nonedited_control, edited_control], axis=-1
            )  # (n_reps, n_bins, n_guides, n_alleles)
        # assert (self.allele_counts_control.sum(axis=-1) == self.X_bcmatch_control).all()
        a0_bcmatch = get_pred_alpha0(
            self.X_bcmatch.clone().cpu(),
            self.size_factor_bcmatch.clone().cpu(),
            self.popt,
            self.sample_mask.cpu(),
        )
        self.a0_bcmatch = torch.as_tensor(a0_bcmatch)

        if use_const_pi:
            self.pi = (self.allele_counts.sum(axis=(0, 1)) + pi_prior_count * 0.5) / (
                self.allele_counts.sum(axis=(0, 1)).sum(axis=-1)[..., None]
                + pi_prior_count
            )  # (n_guides, n_alleles)
            self.allele_counts_control = (
                self.X_bcmatch_control[:, :, :, None] * self.pi[None, None, :, :]
            )
            impute_pi_popt = True
        if impute_pi_popt:
            pi_popt = self.popt
        self.set_fitted_pi_a0(pi_popt, shrink_alpha, use_fitted_value=True)

    def __getitem__(self, guide_idx):
        ndata = super().__getitem__(guide_idx)
        ndata.X_bcmatch = ndata.X_bcmatch[:, :, guide_idx]
        ndata.X_bcmatch_masked = ndata.X_bcmatch_masked[:, :, guide_idx]
        ndata.X_bcmatch_control = ndata.X_bcmatch_control[:, :, guide_idx]
        if hasattr(ndata, "allele_counts"):
            ndata.allele_counts = ndata.allele_counts[:, :, guide_idx, :]
        if hasattr(ndata, "a0_allele"):
            ndata.a0_allele = ndata.a0_allele[guide_idx, :]
        if hasattr(ndata, "pi_a0") and self.pi_a0 is not None:
            ndata.pi_a0 = ndata.pi_a0[guide_idx]
        if hasattr(ndata, "a0_bcmatch") and self.a0_bcmatch is not None:
            ndata.a0_bcmatch = ndata.a0_bcmatch[guide_idx]
        if hasattr(ndata, "pi_a0_bcmatch") and self.pi_a0_bcmatch is not None:
            ndata.pi_a0_bcmatch = ndata.pi_a0_bcmatch[guide_idx]
        return ndata

    def set_fitted_pi_a0(
        self, pi_popt: Tuple, shrink_alpha: bool = False, use_fitted_value: bool = True
    ):
        """Set pi_a0 values

        Args
        pi_popt: pre-fitted values
        shrink_alpha: Return shrinked alpha based on fitted trend, instead of fitted value itself.
        use_fitted_value: use fitted value. If False, use raw and nan values are substituted to fitted values.
        """
        if pi_popt is not None:
            pi_a0 = get_pred_pi_alpha0(
                self.allele_counts_control.clone().cpu(),
                self.size_factor_control.clone().cpu(),
                pi_popt,
            )
        else:
            pi_a0, _pi_popt = get_fitted_pi_alpha0(
                self.allele_counts_control.clone().cpu(),
                self.size_factor_control.clone().cpu(),
                fit_quantile=None,
                shrink=shrink_alpha,
            )
            self._pi_popt = _pi_popt

        if not use_fitted_value:
            pi_a0 = torch.as_tensor(
                get_fitted_pi_alpha0(
                    self.allele_counts_control.clone().cpu(),
                    self.size_factor_control.clone().cpu(),
                    fit=False,
                    fit_quantile=None,
                    shrink=shrink_alpha,
                )[0]
            )
            pi_a0[torch.isnan(pi_a0)] = self.pi_a0[torch.isnan(pi_a0)]
            self.pi_a0 = pi_a0
        self.pi_a0 = torch.as_tensor(pi_a0)


@dataclass
class VariantScreenData(ScreenData):
    n_targets: int
    n_sgRNAs_per_target: int
    target_lengths: torch.Tensor

    def __init__(
        self,
        screen: be.ReporterScreen,
        repguide_mask: str = None,
        sample_mask_column: str = None,
        shrink_alpha: bool = False,
        condition_column: str = "sort",
        control_condition: str = "bulk",
        target_col: str = "target",
        **kwargs,
    ):
        max_n_guides_per_target = screen.guides.groupby(target_col).size().max()
        info(
            f"Data has max {max_n_guides_per_target} guides per target element in '{target_col}' column."
        )

        ScreenData.__init__(
            self,
            screen=screen,
            repguide_mask=repguide_mask,
            sample_mask_column=sample_mask_column,
            shrink_alpha=shrink_alpha,
            condition_column=condition_column,
            control_condition=control_condition,
            **kwargs,
        )
        self.n_targets = len(screen.guides[target_col].unique())
        self.n_sgRNAs_per_target = max_n_guides_per_target
        self.target_lengths = self.get_target_lengths(screen, target_col)
        self.target_col = target_col

    def _post_init(self, target_col):
        self.target_col = target_col
        max_n_guides_per_target = self.screen.guides.groupby(target_col).size().max()
        info(
            f"Data has max {max_n_guides_per_target} guides per target element in '{target_col}' column."
        )
        self.n_targets = len(self.screen.guides[target_col].unique())
        self.n_sgRNAs_per_target = max_n_guides_per_target
        self.target_lengths = self.get_target_lengths(self.screen, self.target_col)

    def __getitem__(self, guide_idx):
        ndata = super().__getitem__(guide_idx)
        ndata.n_targets = len(ndata.screen_selected.guides[self.target_col].unique())
        ndata.target_lengths = ndata.get_target_lengths(
            ndata.screen_selected, self.target_col
        )
        return ndata

    def get_target_lengths(self, screen, target_col="target"):
        target_len_list = []
        screen.guides[target_col] = screen.guides[target_col].astype("category")
        cur_item = screen.guides[target_col].cat.codes.iloc[0]
        cur_len = 0
        n_targets = 0
        for i in screen.guides[target_col].cat.codes:
            if cur_item == i:
                cur_len += 1
            else:
                n_targets += 1
                target_len_list.append(cur_len)
                cur_len = 1
                cur_item = i
        n_targets += 1
        if n_targets != len(screen.guides[target_col].unique()):
            raise ValueError(
                f"Input Screen object not sorted for target identity. Sort the screen object so that guides targeting the same object would occur as consecutive block by screen[screen.guides[{target_col}].argsort(),:]"
            )  # TODO: sort within the script
        target_len_list.append(cur_len)
        target_lengths = torch.tensor(target_len_list)
        return target_lengths


@dataclass
class TilingReporterScreenData(ReporterScreenData):
    def __init__(
        self,
        screen: be.ReporterScreen,
        repguide_mask: str = None,
        sample_mask_column: str = None,
        shrink_alpha: bool = False,
        condition_column: str = "sort",
        control_condition: str = "bulk",
        allele_df_key="allele_counts_filtered",
        edit_index=None,
        alleles_to_select: pd.MultiIndex = None,
        allele_col: str = None,
        **kwargs,
    ):
        """
        Args
            screen: An object of type be.ReporterScreen, representing the screen data.
            repguide_mask: A string specifying the mask for selecting guides x reps DataFrame in screen.uns[repguide_mask].
            sample_mask_column: A string specifying the mask column in screen.samples[sample_mask_column].
            shrink_alpha: A boolean flag that shrinks alpha towards fitted line if True (dev).
            condition_column: Column name in screen.samples that specifies condition that is selected on.
            control_condition: A string specifying the control condition in screen.samples[condition_column].
            allele_df_key: A string representing the key for accessing the allele count dataframe in screen.uns.
            edit_index: A Pandas DataFrame that specifies the edit sequences to index.
            alleles_to_select: A Pandas multi-index representing the alleles that need to be subsetted.
            allele_col: A string representing the name of the allele column in the allele count dataframe (screen.uns[allele_df_key]).
        """
        super().__init__(
            adata=screen,
            repguide_mask=repguide_mask,
            sample_mask_column=sample_mask_column,
            shrink_alpha=shrink_alpha,
            condition_column=condition_column,
            control_condition=control_condition,
            **kwargs,
        )

    def _post_init(
        self,
        allele_df_key: str = None,
        allele_col: str = None,
        alleles_to_select: pd.MultiIndex = None,
        edit_index: pd.DataFrame = None,
        control_guide_tag: str = None,
    ):
        """Called in initialization of child classes to add tiling-screen specific attributes.

        Args
            allele_df_key: A string representing the allele count dataframe key in screen.uns.
            allele_col: A string representing the name of the allele column in the allele count dataframe.
            alleles_to_select: A Pandas multi-index representing the alleles that need to be subsetted.
            edit_index: A Pandas DataFrame that specifies the edit sequences to index.
        """
        if allele_col is None:
            if "allele" in self.screen.uns[allele_df_key].reset_index():
                allele_col = "allele"
            if "aa_allele" in self.screen.uns[allele_df_key].reset_index():
                allele_col = "aa_allele"
        if alleles_to_select is not None:
            allele_counts_selected = (
                self.screen.uns[allele_df_key]
                .set_index(["guide", allele_col])
                .loc[alleles_to_select, :]
                .reset_index()
            )
        else:
            allele_counts_selected = self.screen.uns[allele_df_key]
        if control_guide_tag is not None:

            def _set_uid_to_row(row):
                if control_guide_tag in row.guide:
                    row[allele_col].set_uid(row.guide)

            allele_counts_selected.apply(lambda row: _set_uid_to_row(row), axis=1)
            assert (
                allele_counts_selected[allele_col].map(lambda a: "!" in str(a)).any()
            ), "uid not assinged."
        guide_to_allele, reindexed_df = self.reindex_allele_df(
            allele_counts_selected, allele_col
        )  # TODO: Fix this for combining different sorting scheme
        self.n_max_alleles = (
            reindexed_df.index.get_level_values("allele_id_for_guide").max() + 1
        )  # include no edit allele

        if edit_index is None:
            self.edit_index = get_edit_to_index_dict(guide_to_allele[allele_col])
        else:
            self.edit_index = edit_index
        self.n_edits = len(self.edit_index.keys())

        self.allele_to_edit = self.get_allele_to_edit_tensor(
            self.screen_control, self.edit_index, guide_to_allele, allele_col
        )
        assert self.allele_to_edit.shape == (
            self.n_guides,
            self.n_max_alleles - 1,
            self.n_edits,
        )

        # TODO: fix this for sample-separated data
        self.allele_counts = self.transform_allele(self.screen_selected, reindexed_df)
        assert self.allele_counts.shape == (
            self.n_reps,
            self.n_condits,
            self.n_guides,
            self.n_max_alleles,
        )
        self.allele_counts_control = self.transform_allele_control(
            self.screen_control, reindexed_df
        )
        self.allele_mask = self.get_allele_mask(self.screen_control, guide_to_allele)
        assert self.allele_mask.shape == (self.n_guides, self.n_max_alleles)
        assert self.allele_counts_control.shape == (
            self.n_reps,
            1,
            self.n_guides,
            self.n_max_alleles,
        )

    def get_allele_to_edit_tensor(
        self,
        screen,
        edits_to_index: Dict[str, int],
        guide_allele_id_to_allele_df: pd.DataFrame,
        allele_col="aa_allele",
    ) -> torch.Tensor:
        """
        Convert (guide, allele_id_for_guide) -> allele DataFrame into the tensor with shape (n_guides, n_max_alleles_per_guide, n_edits) tensor.

        Args
            edits_to_index: Dictionary from edit (str) to unique index (int)
            guide_allele_id_to_allele_df: pd.DataFrame of (guide(str), allele_id_for_guide(int)) -> CodingNoncodingAllele

        Returns
            allele_edit_assignment: Binary tensor of shape (n_guides, n_max_alleles_per_guide, n_edits. allele_edit_assignment(i, j, k) is 1 if jth allele of ith guide has kth edit.
        """
        if allele_col == "aa_allele":
            guide_allele_id_to_allele_df["edits"] = guide_allele_id_to_allele_df[
                allele_col
            ].map(lambda a: list(a.aa_allele.edits) + list(a.nt_allele.edits))
        else:
            guide_allele_id_to_allele_df["edits"] = guide_allele_id_to_allele_df[
                allele_col
            ].map(lambda a: list(a.edits))
        guide_allele_id_to_allele_df = guide_allele_id_to_allele_df.reset_index()
        guide_allele_id_to_allele_df["edit_idx"] = (
            guide_allele_id_to_allele_df.edits.map(
                lambda es: [edits_to_index[e.get_abs_edit()] for e in es]
            )
        )
        guide_allele_id_to_edit_df = guide_allele_id_to_allele_df[
            ["guide", "allele_id_for_guide", "edit_idx"]
        ].set_index(["guide", "allele_id_for_guide"])
        guide_allele_id_to_edit_df = guide_allele_id_to_edit_df.unstack(
            level=1, fill_value=[]
        ).reindex(screen.guides.index, fill_value=[])
        allele_edit_assignment = torch.zeros(
            (len(screen.guides), self.n_max_alleles - 1, len(edits_to_index.keys()))
        )
        for i in range(len(guide_allele_id_to_edit_df)):
            for j in range(len(guide_allele_id_to_edit_df.columns)):
                allele_edit_assignment[i, j, guide_allele_id_to_edit_df.iloc[i, j]] = 1
        return allele_edit_assignment

    def reindex_allele_df(self, alleles_df, allele_col):
        """
        Input: Dataframe of (guide, allele) -> (per sample count)
        Output:
            - guide_allele_id_to_allele
                DataFrame of (guide, allele_id_for_guide) -> (global_allele_id, aa_allele(str))
            - reindexed_allele_df
                Dataframe of (guide, allele_id_for_guide) -> (per sample count)

            allele_id_for_guide: order of allele within a guide.
            global_allele_id: global unique id for each (guide, allele) pair.
        """
        guide_to_allele = dict(
            list(alleles_df[["guide", allele_col]].groupby("guide")[allele_col])
        )
        dfs = []
        for k, s in guide_to_allele.items():
            df = s.reset_index()
            df = df.rename(columns={"index": "allele_id"})
            df.index = df.index + 1
            df = df.reset_index()
            df = df.rename(columns={"index": "allele_id_for_guide"})
            df["guide"] = k
            dfs.append(df)

        guide_to_allele_tbl = pd.concat(dfs)

        alleles_df = pd.merge(alleles_df, guide_to_allele_tbl, on=[allele_col, "guide"])
        reindexed_df = alleles_df.reset_index().set_index(
            ["guide", "allele_id_for_guide"]
        )
        guide_allele_id_to_allele = reindexed_df[["index", allele_col]]
        reindexed_allele_df = reindexed_df.drop([allele_col, "index"], axis=1)
        return (guide_allele_id_to_allele, reindexed_allele_df)

    def transform_allele(self, adata, reindexed_df):
        """
        Transform reindexed allele dataframe reindexed_df of (guide, allele_id_for_guide) -> (per sample count)
        to (n_reps, n_bins, n_guides, n_alleles) tensor.
        """
        allele_tensor = torch.empty(
            (self.n_reps, self.n_condits, self.n_guides, self.n_max_alleles),
        )
        if self.device is not None:
            allele_tensor = allele_tensor.cuda()
        for i in range(self.n_reps):
            for j in range(self.n_condits):
                condit_idx = np.where(
                    (adata.samples.replicate_id == i)
                    & (adata.samples[f"{self.condition_column}_id"] == j)
                )[0]
                assert len(condit_idx) == 1, print(i, j, condit_idx)
                condit_idx = condit_idx.item()
                condit_name = adata.samples.index[condit_idx]
                condit_allele_df = (
                    reindexed_df.loc[:, condit_name]
                    .unstack(level=1, fill_value=0)
                    .astype(int)
                )
                condit_allele_df = condit_allele_df.reindex(
                    adata.guides.index, fill_value=0
                )
                condit_bcmatch_counts = adata.layers["X_bcmatch"][:, condit_idx].astype(
                    int
                )
                # if not (condit_bcmatch_counts >= condit_allele_df.sum(axis=1)).all():
                #     print(
                #         f"Allele counts are larger than total bcmatch counts in rep {i}, {j} by {(condit_bcmatch_counts - condit_allele_df.sum(axis=1)).min()}."
                #     )
                condit_allele_df[0] = condit_bcmatch_counts - condit_allele_df.loc[
                    :, condit_allele_df.columns != 0
                ].sum(axis=1)
                condit_allele_df[0].loc[condit_allele_df[0] < 0] = 0
                condit_allele_df = condit_allele_df.sort_values(
                    "allele_id_for_guide", axis=1, ascending=True
                )
                assert all(condit_allele_df.columns == list(range(self.n_max_alleles)))
                assert condit_allele_df.to_numpy().shape == (
                    self.n_guides,
                    self.n_max_alleles,
                )
                allele_tensor[i, j, :, :] = torch.as_tensor(condit_allele_df.to_numpy())

        try:
            assert (allele_tensor >= 0).all(), allele_tensor[allele_tensor < 0]
        except AssertionError:
            print("Allele tensor doesn't match condit_allele_df")
            return (allele_tensor, reindexed_df)
        return allele_tensor

    def transform_allele_control(self, adata, reindexed_df):
        """
        Transform reindexed allele dataframe reindexed_df of (guide, allele_id_for_guide) -> (per sample count)
        to (n_reps, n_bins, n_guides, n_alleles) tensor.
        """

        allele_tensor = torch.empty(
            (self.n_reps, 1, self.n_guides, self.n_max_alleles),
        )
        if self.device is not None:
            allele_tensor = allele_tensor.cuda()
        for i in range(self.n_reps):
            condit_idx = np.where(adata.samples.replicate_id == i)[0]
            assert len(condit_idx) == 1, print(i, self.control_condition, condit_idx)
            condit_idx = condit_idx.item()
            condit_name = adata.samples.index[condit_idx]
            condit_allele_df = (
                reindexed_df.loc[:, condit_name]
                .unstack(level=1, fill_value=0)
                .astype(int)
            )
            condit_allele_df = condit_allele_df.reindex(
                adata.guides.index, fill_value=0
            )
            condit_bcmatch_counts = adata.layers["X_bcmatch"][:, condit_idx].astype(int)
            # if not (condit_bcmatch_counts >= condit_allele_df.sum(axis=1)).all():
            #     print(
            #         "Allele counts are larger than total bcmatch counts in rep {}, .".format(
            #             i,
            #         )
            #     )
            #     print(
            #         f"min value: {(condit_bcmatch_counts - condit_allele_df.sum(axis=1)).min()}"
            #     )
            condit_allele_df[0] = condit_bcmatch_counts - condit_allele_df.loc[
                :, condit_allele_df.columns != 0
            ].sum(axis=1)
            condit_allele_df[0].loc[condit_allele_df[0] < 0] = 0
            condit_allele_df = condit_allele_df.sort_values(
                "allele_id_for_guide", axis=1, ascending=True
            )
            assert all(condit_allele_df.columns == list(range(self.n_max_alleles)))
            assert condit_allele_df.to_numpy().shape == (
                self.n_guides,
                self.n_max_alleles,
            )
            allele_tensor[i, 0, :, :] = torch.as_tensor(condit_allele_df.to_numpy())
        try:
            assert (allele_tensor >= 0).all(), allele_tensor[allele_tensor < 0]
        except AssertionError:
            print("Negative values in allele_tensor")
            return (allele_tensor, reindexed_df)
        return allele_tensor

    def get_allele_mask(
        self, adata, guide_allele_id_to_allele_df, keep_one_allele: bool = False
    ):
        """
        =====
        Arguments
        keep_one_allele -- If True, if no valid allele exist, still provide mask for one valid allele.
        """
        guide_to_allele_tbl = guide_allele_id_to_allele_df.reset_index()
        n_valid_allele = (
            guide_to_allele_tbl.groupby("guide")["allele_id_for_guide"]
            .max()
            .reindex(adata.guides.index, fill_value=0)
        )
        if keep_one_allele:
            n_valid_allele[n_valid_allele == 0] = 1
        mask = torch.zeros((self.n_guides, self.n_max_alleles))
        for i in range(mask.shape[0]):
            mask[i, 0] = 1
            for j in range(n_valid_allele[i]):
                mask[i, j + 1] = 1
        return mask.bool()


@dataclass
class SortingScreenData(ScreenData):
    upper_bounds: torch.Tensor
    lower_bounds: torch.Tensor

    def __init__(
        self,
        screen: be.ReporterScreen,
        repguide_mask: str = None,
        sample_mask_column: str = None,
        shrink_alpha: bool = False,
        condition_column: str = "sort",
        sample_covariate_column: List[str] = [],
        control_condition: str = "bulk",
        lower_quantile_column: str = "lower_quantile",
        upper_quantile_column: str = "upper_quantile",
        replicate_column: str = "replicate",
        **kwargs,
    ):
        """
        Args
        ---
        lower_quantile_column: column in screen.samples that indicate lower quantile threshold of sorting bin.
        upper_quantile_column: column in screen.samples that indicate upper quantile threshold of sorting bin.
        """
        self._pre_init(lower_quantile_column, upper_quantile_column, replicate_column)
        self.n_bins = self.n_condits
        super().__init__(
            screen=screen,
            repguide_mask=repguide_mask,
            sample_mask_column=sample_mask_column,
            shrink_alpha=shrink_alpha,
            condition_column=condition_column,
            control_condition=control_condition,
            **kwargs,
        )
        self.condition_column = condition_column

    def _pre_init(
        self,
        lower_quantile_column,
        upper_quantile_column,
    ):
        if (self.screen.samples[lower_quantile_column] < 0.0).any() or (
            self.screen.samples[lower_quantile_column] > 1.0
        ).any():
            raise ValueError(
                f"Invalid quantile value({self.screen.samples[lower_quantile_column]}) in screen.samples[{lower_quantile_column}]: check input."
            )
        if (self.screen.samples[upper_quantile_column] < 0.0).any() or (
            self.screen.samples[upper_quantile_column] > 1.0
        ).any():
            raise ValueError(
                f"Invalid quantile value ({self.screen.samples[upper_quantile_column]}) in screen.samples[{upper_quantile_column}]: check input."
            )
        if (
            self.screen.samples[upper_quantile_column]
            - self.screen.samples[lower_quantile_column]
            < 0
        ).any():
            raise ValueError(
                f"Not all screen.samples[{upper_quantile_column}] larger than screen.samples[{upper_quantile_column}]: check input."
            )

        if not (
            self.screen.samples.groupby(
                [upper_quantile_column, lower_quantile_column]
            ).size()
            == len(self.screen.samples[self.replicate_column].unique())
        ).all():
            raise ValueError(
                "Not all replicate share same quantile bin definition. If you have missing bin data, add the sample and add 'mask' column in 'screen.samples' or run `bean-qc` that automatically handles this."
            )
        sorting_bins = self.screen_selected.samples.sort_values(
            [upper_quantile_column, lower_quantile_column]
        )[[upper_quantile_column, lower_quantile_column]].drop_duplicates()
        for j, (idx, row) in enumerate(sorting_bins.iterrows()):
            self.screen_selected.samples.loc[
                (
                    self.screen_selected.samples[upper_quantile_column]
                    == row[upper_quantile_column]
                )
                & (
                    self.screen_selected.samples[lower_quantile_column]
                    == row[lower_quantile_column]
                ),
                f"{self.condition_column}_id",
            ] = j
        self.upper_bounds = torch.as_tensor(sorting_bins[upper_quantile_column].values)
        self.lower_bounds = torch.as_tensor(sorting_bins[lower_quantile_column].values)
        self.screen.samples[f"{self.condition_column}_id"] = -1
        self.screen.samples.loc[
            self.screen_selected.samples.index, f"{self.condition_column}_id"
        ] = self.screen_selected.samples[f"{self.condition_column}_id"]
        self.screen = _assign_rep_ids_and_sort(
            self.screen, self.replicate_column, self.condition_column
        )
        if hasattr(self, "sample_covariates"):
            self.rep_by_cov = torch.as_tensor(
                (
                    self.screen.samples[["_rc"] + self.sample_covariates]
                    .drop_duplicates()
                    .set_index("_rc")
                    .values.astype(int)
                )
            )
        self.screen_selected = _assign_rep_ids_and_sort(
            self.screen_selected, self.replicate_column, self.condition_column
        )
        self.screen_control = _assign_rep_ids_and_sort(
            self.screen_control,
            self.replicate_column,
        )


@dataclass
class SurvivalScreenData(ScreenData):
    def __init__(
        self,
        screen: be.ReporterScreen,
        repguide_mask: str = None,
        sample_mask_column: str = None,
        shrink_alpha: bool = False,
        condition_column: str = "condition",
        control_condition: str = "bulk",
        control_can_be_selected=True,
        time_column: str = "time",
        replicate_column: str = "replicate",
        **kwargs,
    ):
        self._pre_init(condition_column)
        super().__init__(
            screen=screen,
            repguide_mask=repguide_mask,
            sample_mask_column=sample_mask_column,
            shrink_alpha=shrink_alpha,
            condition_column=condition_column,
            control_condition=control_condition,
            control_can_be_selected=control_can_be_selected,
            **kwargs,
        )
        self._post_init()

    def _pre_init(self, time_column: str, condition_column: str):
        self.condition_column = self.time_column = time_column
        try:
            self.screen.samples[time_column] = self.screen.samples[time_column].astype(
                float
            )
        except ValueError as e:
            raise ValueError(
                f"Invalid timepoint value({self.screen.samples[time_column]}) in screen.samples[{time_column}]: check input."
            ) from e

        if not (
            self.screen.samples.groupby(condition_column).size()
            == len(self.screen.samples[self.replicate_column].unique())
        ).all():
            raise ValueError(
                f"Not all replicate share same timepoint definition. If you have missing bin data, add the sample and add 'mask' column in 'screen.samples', or run `bean-qc` that automatically handles this. \n{self.screen.samples}"
            )

    def _post_init(
        self,
    ):
        self.timepoints = torch.as_tensor(
            self.screen_selected.samples[self.time_column].unique()
        )
        control_timepoint = self.screen_control.samples[self.time_column].unique()
        if len(control_timepoint) != 1:
            info(self.screen_control)
            info(self.screen_control.samples)
            info(control_timepoint)
            raise ValueError(
                "All samples with --control-condition should have the same --time-col column in ReporterScreen.samples[time_col]. Check your input ReporterScreen object."
            )
        else:
            self.control_timepoint = control_timepoint[0]
        self.n_timepoints = self.n_condits
        timepoints = self.screen_selected.samples.sort_values(self.time_column)[
            self.time_column
        ].drop_duplicates()
        if timepoints.isnull().any():
            raise ValueError(
                f"NaN values in time points provided in input: {self.screen_selected.samples[self.time_column]}"
            )
        for j, time in enumerate(timepoints):
            self.screen_selected.samples.loc[
                self.screen_selected.samples[self.time_column] == time,
                f"{self.time_column}_id",
            ] = j
        self.screen.samples[f"{self.time_column}_id"] = -1
        self.screen.samples.loc[
            self.screen_selected.samples.index, f"{self.time_column}_id"
        ] = self.screen_selected.samples[f"{self.time_column}_id"]
        self.screen = _assign_rep_ids_and_sort(
            self.screen, self.replicate_column, self.time_column
        )
        if hasattr(self, "sample_covariates"):
            self.rep_by_cov = self.screen.samples.groupby(self.replicate_column)[
                self.sample_covariates
            ].values
        self.screen_selected = _assign_rep_ids_and_sort(
            self.screen_selected,
            self.replicate_column,
            self.time_column,
        )
        self.screen_control = _assign_rep_ids_and_sort(
            self.screen_control,
            self.replicate_column,
        )


@dataclass
class VariantReporterScreenData(VariantScreenData, ReporterScreenData):
    def __init__(
        self,
        screen: be.ReporterScreen,
        repguide_mask: str = None,
        sample_mask_column: str = None,
        pi_popt: Tuple[float] = None,
        impute_pi_popt: bool = False,
        shrink_alpha: bool = False,
        condition_column: str = "condition",
        control_condition: str = "bulk",
        accessibility_col: str = None,
        accessibility_bw_path: str = None,
        use_const_pi: bool = False,
        pi_prior_count: int = 10,
        **kwargs,
    ):
        VariantScreenData.__init__(
            self,
            screen=screen,
            repguide_mask=repguide_mask,
            sample_mask_column=sample_mask_column,
            condition_column=condition_column,
            control_condition=control_condition,
            accessibility_col=accessibility_col,
            accessibility_bw_path=accessibility_bw_path,
            **kwargs,
        )
        ReporterScreenData._post_init(
            self,
            screen,
            use_const_pi,
            impute_pi_popt,
            pi_prior_count,
            shrink_alpha,
            pi_popt,
        )


@dataclass
class VariantSortingScreenData(VariantScreenData, SortingScreenData):
    def __init__(
        self,
        screen,
        *args,
        lower_quantile_column="lower_quantile",
        upper_quantile_column="upper_quantile",
        replicate_column="replicate",
        condition_column="bin",
        target_col="target",
        sample_mask_column="mask",
        shrink_alpha: bool = False,
        use_bcmatch=False,
        **kwargs,
    ):
        ScreenData.__init__(
            self,
            screen,
            *args,
            sample_mask_column=sample_mask_column,
            replicate_column=replicate_column,
            condition_column=condition_column,
            shrink_alpha=shrink_alpha,
            **kwargs,
        )
        SortingScreenData._pre_init(
            self,
            lower_quantile_column,
            upper_quantile_column,
        )
        ScreenData._post_init(self)
        VariantScreenData._post_init(self, target_col)
        if use_bcmatch:
            self.set_bcmatch(
                screen,
            )

    def set_bcmatch(self, screen):
        screen.samples["size_factor_bcmatch"] = self.get_size_factor(
            screen.layers["X_bcmatch"]
        )
        self.screen_selected.samples["size_factor_bcmatch"] = screen.samples.loc[
            self.screen_selected.samples.index, "size_factor_bcmatch"
        ]
        self.screen_control.samples["size_factor_bcmatch"] = screen.samples.loc[
            self.screen_control.samples.index, "size_factor_bcmatch"
        ]
        self.X_bcmatch = self.transform_data(self.screen_selected.layers["X_bcmatch"])
        self.X_bcmatch_masked = self.X_bcmatch * self.sample_mask[:, :, None]
        self.X_bcmatch_control = self.transform_data(
            self.screen_control.layers["X_bcmatch"], 1
        )
        self.X_bcmatch_control_masked = (
            self.X_bcmatch_control * self.bulk_sample_mask[:, None, None]
        )
        self.size_factor_bcmatch = torch.as_tensor(
            self.screen_selected.samples["size_factor_bcmatch"].to_numpy()
        ).reshape(self.n_reps, self.n_condits)
        self.size_factor_bcmatch_control = torch.as_tensor(
            self.screen_control.samples["size_factor_bcmatch"].to_numpy()
        ).reshape(self.n_reps, 1)
        a0_bcmatch = get_pred_alpha0(
            self.X_bcmatch.clone().cpu(),
            self.size_factor_bcmatch.clone().cpu(),
            self.popt,
            self.sample_mask.cpu(),
        )
        self.a0_bcmatch = torch.as_tensor(a0_bcmatch)

    def __getitem__(self, guide_idx):
        ndata = super().__getitem__(guide_idx)
        if hasattr(ndata, "X_bcmatch"):
            ndata.X_bcmatch = ndata.X_bcmatch[:, :, guide_idx]
        if hasattr(ndata, "X_bcmatch_masked"):
            ndata.X_bcmatch_masked = ndata.X_bcmatch_masked[:, :, guide_idx]
        if hasattr(ndata, "X_bcmatch_control"):
            ndata.X_bcmatch_control = ndata.X_bcmatch_control[:, :, guide_idx]
        if hasattr(ndata, "X_bcmatch_control_masked"):
            ndata.X_bcmatch_control_masked = ndata.X_bcmatch_control_masked[
                :, :, guide_idx
            ]
        return ndata


@dataclass
class VariantSortingReporterScreenData(VariantReporterScreenData, SortingScreenData):
    def __init__(
        self,
        screen,
        *args,
        lower_quantile_column="lower_quantile",
        upper_quantile_column="upper_quantile",
        replicate_column="replicate",
        condition_column="bin",
        target_col="target",
        sample_mask_column="mask",
        use_const_pi: bool = False,
        impute_pi_popt: bool = False,
        pi_prior_count: int = 10,
        shrink_alpha: bool = False,
        pi_popt: Tuple[float] = None,
        **kwargs,
    ):
        ScreenData.__init__(
            self,
            screen,
            *args,
            sample_mask_column=sample_mask_column,
            replicate_column=replicate_column,
            condition_column=condition_column,
            shrink_alpha=shrink_alpha,
            **kwargs,
        )
        SortingScreenData._pre_init(
            self,
            lower_quantile_column,
            upper_quantile_column,
        )
        ScreenData._post_init(self)
        VariantScreenData._post_init(self, target_col)
        ReporterScreenData._post_init(
            self,
            screen,
            use_const_pi,
            impute_pi_popt,
            pi_prior_count,
            shrink_alpha,
            pi_popt,
        )


@dataclass
class TilingSortingReporterScreenData(TilingReporterScreenData, SortingScreenData):
    def __init__(
        self,
        screen,
        *args,
        lower_quantile_column="lower_quantile",
        upper_quantile_column="upper_quantile",
        replicate_column="replicate",
        condition_column="bin",
        sample_mask_column="mask",
        use_const_pi: bool = False,
        impute_pi_popt: bool = False,
        pi_prior_count: int = 10,
        shrink_alpha: bool = False,
        pi_popt: Tuple[float] = None,
        allele_df_key: str = None,
        allele_col: str = None,
        control_guide_tag: str = None,
        **kwargs,
    ):
        ScreenData.__init__(
            self,
            screen,
            *args,
            sample_mask_column=sample_mask_column,
            replicate_column=replicate_column,
            condition_column=condition_column,
            shrink_alpha=shrink_alpha,
            **kwargs,
        )
        SortingScreenData._pre_init(
            self,
            lower_quantile_column,
            upper_quantile_column,
        )
        ScreenData._post_init(self)
        TilingReporterScreenData._post_init(
            self,
            allele_df_key=allele_df_key,
            control_guide_tag=control_guide_tag,
        )
        ReporterScreenData._post_init(
            self,
            screen,
            use_const_pi,
            impute_pi_popt,
            pi_prior_count,
            shrink_alpha,
            pi_popt,
        )


@dataclass
class VariantSurvivalScreenData(VariantScreenData, SurvivalScreenData):
    def __init__(
        self,
        screen,
        *args,
        replicate_column="replicate",
        condition_column="condition",
        time_column="time",
        control_can_be_selected=True,
        target_col="target",
        sample_mask_column="mask",
        shrink_alpha: bool = False,
        use_bcmatch=False,
        **kwargs,
    ):
        ScreenData.__init__(
            self,
            screen,
            *args,
            sample_mask_column=sample_mask_column,
            replicate_column=replicate_column,
            condition_column=condition_column,
            time_column=time_column,
            shrink_alpha=shrink_alpha,
            control_can_be_selected=control_can_be_selected,
            **kwargs,
        )
        SurvivalScreenData._pre_init(self, time_column, condition_column)
        ScreenData._post_init(self)
        SurvivalScreenData._post_init(self)
        VariantScreenData._post_init(self, target_col)
        if use_bcmatch:
            self.set_bcmatch(
                screen,
            )

    def __getitem__(self, guide_idx):
        ndata = super().__getitem__(guide_idx)
        if hasattr(ndata, "X_bcmatch"):
            ndata.X_bcmatch = ndata.X_bcmatch[:, :, guide_idx]
        if hasattr(ndata, "X_bcmatch_masked"):
            ndata.X_bcmatch_masked = ndata.X_bcmatch_masked[:, :, guide_idx]
        if hasattr(ndata, "X_bcmatch_control"):
            ndata.X_bcmatch_control = ndata.X_bcmatch_control[:, :, guide_idx]
        if hasattr(ndata, "X_bcmatch_control_masked"):
            ndata.X_bcmatch_control_masked = ndata.X_bcmatch_control_masked[
                :, :, guide_idx
            ]
        return ndata

    def set_bcmatch(self, screen):
        screen.samples["size_factor_bcmatch"] = self.get_size_factor(
            screen.layers["X_bcmatch"]
        )
        self.screen_selected.samples["size_factor_bcmatch"] = screen.samples.loc[
            self.screen_selected.samples.index, "size_factor_bcmatch"
        ]
        self.screen_control.samples["size_factor_bcmatch"] = screen.samples.loc[
            self.screen_control.samples.index, "size_factor_bcmatch"
        ]
        self.X_bcmatch = self.transform_data(self.screen_selected.layers["X_bcmatch"])
        self.X_bcmatch_masked = self.X_bcmatch * self.sample_mask[:, :, None]
        self.X_bcmatch_control = self.transform_data(
            self.screen_control.layers["X_bcmatch"], 1
        )
        self.X_bcmatch_control_masked = (
            self.X_bcmatch_control * self.bulk_sample_mask[:, None, None]
        )
        self.size_factor_bcmatch = torch.as_tensor(
            self.screen_selected.samples["size_factor_bcmatch"].to_numpy()
        ).reshape(self.n_reps, self.n_condits)
        self.size_factor_bcmatch_control = torch.as_tensor(
            self.screen_control.samples["size_factor_bcmatch"].to_numpy()
        ).reshape(self.n_reps, 1)
        a0_bcmatch = get_pred_alpha0(
            self.X_bcmatch.clone().cpu(),
            self.size_factor_bcmatch.clone().cpu(),
            self.popt,
            self.sample_mask.cpu(),
        )
        self.a0_bcmatch = torch.as_tensor(a0_bcmatch)


@dataclass
class VariantSurvivalReporterScreenData(VariantReporterScreenData, SurvivalScreenData):
    def __init__(
        self,
        screen,
        *args,
        replicate_column="replicate",
        condition_column="condition",
        time_column="time",
        control_can_be_selected=True,
        target_col="target",
        sample_mask_column="mask",
        use_const_pi: bool = False,
        impute_pi_popt: bool = False,
        pi_prior_count: int = 10,
        shrink_alpha: bool = False,
        pi_popt: Tuple[float] = None,
        **kwargs,
    ):
        ScreenData.__init__(
            self,
            screen,
            *args,
            sample_mask_column=sample_mask_column,
            replicate_column=replicate_column,
            condition_column=condition_column,
            time_column=time_column,
            shrink_alpha=shrink_alpha,
            control_can_be_selected=control_can_be_selected,
            **kwargs,
        )
        SurvivalScreenData._pre_init(self, time_column, condition_column)
        ScreenData._post_init(self)
        SurvivalScreenData._post_init(self)
        VariantScreenData._post_init(self, target_col)
        ReporterScreenData._post_init(
            self,
            screen,
            use_const_pi,
            impute_pi_popt,
            pi_prior_count,
            shrink_alpha,
            pi_popt,
        )


@dataclass
class TilingSurvivalReporterScreenData(TilingReporterScreenData, SurvivalScreenData):
    def __init__(
        self,
        screen,
        *args,
        replicate_column="replicate",
        condition_column="condition",
        time_column="time",
        control_can_be_selected=True,
        sample_mask_column="mask",
        use_const_pi: bool = False,
        impute_pi_popt: bool = False,
        pi_prior_count: int = 10,
        shrink_alpha: bool = False,
        pi_popt: Tuple[float] = None,
        allele_df_key: str = None,
        allele_col: str = None,
        control_guide_tag: str = None,
        **kwargs,
    ):
        ScreenData.__init__(
            self,
            screen,
            *args,
            sample_mask_column=sample_mask_column,
            replicate_column=replicate_column,
            condition_column=condition_column,
            time_column=time_column,
            shrink_alpha=shrink_alpha,
            control_can_be_selected=control_can_be_selected,
            **kwargs,
        )
        SurvivalScreenData._pre_init(self, time_column, condition_column)
        ScreenData._post_init(self)
        SurvivalScreenData._post_init(self)
        TilingReporterScreenData._post_init(
            self,
            allele_df_key=allele_df_key,
            control_guide_tag=control_guide_tag,
        )
        ReporterScreenData._post_init(
            self,
            screen,
            use_const_pi,
            impute_pi_popt,
            pi_prior_count,
            shrink_alpha,
            pi_popt,
        )


DATACLASS_DICT = {
    "sorting": {
        "Normal": VariantSortingScreenData,
        "MixtureNormal": VariantSortingReporterScreenData,
        "MixtureNormal+Acc": VariantSortingReporterScreenData,
        "MixtureNormalConstPi": VariantSortingScreenData,
        "MultiMixtureNormal": TilingSortingReporterScreenData,
        "MultiMixtureNormal+Acc": TilingSortingReporterScreenData,
    },
    "survival": {
        "Normal": VariantSurvivalScreenData,
        "MixtureNormal": VariantSurvivalReporterScreenData,
        "MixtureNormal+Acc": VariantSurvivalReporterScreenData,
        "MultiMixtureNormal": TilingSurvivalReporterScreenData,
        "MultiMixtureNormal+Acc": TilingSurvivalReporterScreenData,
    },
}
