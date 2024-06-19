import pickle as pkl
import numpy as np
import torch
from bean.model.run import _get_guide_target_info
from bean.model.parser import parse_args
from bean.cli.run import main as get_screendata
from bean.preprocessing.data_class import SortingScreenData


def generate_prior_data_for_disjoint_library_pair(
    command1: str, command2: str, output1_path: str, prior_params_path: str
):
    """Generate prior for a two batches with disjoint guides but with shared variants."""
    with open(output1_path, "rb") as f:
        data = pkl.load(f)
        ndata = data["data"]
    parser = parse_args()
    command1 = command1.split("bean run ")[-1]
    command2 = command2.split("bean run ")[-1]
    args = parser.parse_args(command1.split(" "))
    args2 = parser.parse_args(command2.split(" "))
    ndata2 = get_screendata(args2, return_data=True)
    target_df = _get_guide_target_info(
        ndata.screen, args, cols_include=[args.negctrl_col]
    )
    target_df2 = _get_guide_target_info(
        ndata2.screen, args2, cols_include=[args2.negctrl_col]
    )
    batch1_idx = np.where(
        target_df.index.map(lambda s: s in target_df2.index.tolist())
    )[0]
    batch2_idx = []
    for i in batch1_idx:
        batch2_idx.append(
            np.where(target_df.index.tolist()[i] == target_df2.index)[0].item()
        )
    batch2_idx = np.array(batch2_idx)
    if isinstance(ndata, SortingScreenData):
        mu_loc = torch.zeros((ndata2.n_targets, 1))
        mu_loc[batch2_idx, :] = data["params"]["mu_loc"][batch1_idx, :]
        mu_scale = torch.ones((ndata2.n_targets, 1))
        mu_scale[batch2_idx, :] = data["params"]["mu_scale"][batch1_idx, :]
        sd_loc = torch.zeros((ndata2.n_targets, 1))
        sd_loc[batch2_idx, :] = data["params"]["sd_loc"][batch1_idx, :]
        sd_scale = torch.ones((ndata2.n_targets, 1)) * 0.01
        sd_scale[batch2_idx, :] = data["params"]["sd_scale"][batch1_idx, :]
        prior_params = {
            "mu_loc": mu_loc,
            "mu_scale": mu_scale,
            "sd_loc": sd_loc,
            "sd_scale": sd_scale,
        }
    else:
        mu_loc = torch.zeros((ndata2.n_targets, 1))
        mu_loc[batch2_idx, :] = data["params"]["mu_loc"][batch1_idx, :]
        mu_scale = torch.ones((ndata2.n_targets, 1))
        mu_scale[batch2_idx, :] = data["params"]["mu_scale"][batch1_idx, :]
        prior_params = {
            "mu_loc": mu_loc,
            "mu_scale": mu_scale,
        }
    with open(prior_params_path, "wb") as f:
        pkl.dump(prior_params, f)
    print(
        f"Successfully generated prior parameters at {prior_params_path}. To use this parameter, run:\nbean run {command2+' --prior-params '+prior_params_path}"
    )


def main(args):
    generate_prior_data_for_disjoint_library_pair(
        args.command1, args.command2, args.raw_run_output1, args.output_path
    )
