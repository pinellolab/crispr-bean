import matplotlib.pyplot as plt
from simulate import (
    simulate_reps,
    _get_effect_sizes
    )

def write_mageck_input(
    outfile_name, 
    total_cells_per_rep = 10**5, 
    nreps = 4, 
    nreads_per_sample = 10**6,
    sorting_mode = "topbot"
):
    if sorting_mode == "topbot":
        sorting_bins = ["top", "bot"]
    elif sorting_mode == "bins":
        sorting_bins = ["top", "high", "low", "bot"]
    measures = ["guide", "target_edit", "reporter_edit"]
    res_reps = simulate_reps(total_cells_per_rep, nreps, nreads_per_sample, 
                            sorting_mode = sorting_mode)
    to_merge = [x.set_index(
        ["target_id", "guide_num"]
    )[["{}_{}".format(sb, m) for sb in sorting_bins + ["bulk"] for m in measures]
    ].add_prefix("rep{}_".format(i)) for i, x in enumerate(res_reps)]
    
    count_tbl = reduce(lambda left, right: left.join(right), to_merge)
    count_tbl = count_tbl.reset_index().fillna(0)

    count_tbl.insert(0, "sgRNA", count_tbl.reset_index().target_id.astype(str) + "_" + count_tbl.guide_num.astype(str))
    count_tbl = count_tbl.rename(columns = {"target_id":"gene"})
    count_tbl = count_tbl.drop('guide_num', axis = 1).astype(int)
    count_tbl.to_csv(outfile_name, sep = "\t", index = False)
    
    return(count_tbl)

def run_mageck(sorting_mode = "topbot",
               total_cells_per_rep = 10**5, 
               nreps = 4, 
               nreads_per_sample = 10**6,
               measure = "guide", 
               iteration_per_screen = 10, 
               rerun = True, 
              ):
    """Simulate & save the simulated counts and run mageck on it."""
    assert measure in ["guide", "target_edit", "reporter_edit"]
    if sorting_mode == "topbot":
        sorting_bins = ["top", "bot"]
    elif sorting_mode == "bins":
        sorting_bins = ["top", "high", "low", "bot"]
    mageck_reses = None
    mean_sensitivities = None
    
    for r in range(iteration_per_screen):
        mageck_input_name = "sorted_reads/mageck_input_{}_{}_{}_{}_{}.tsv".format(
            sorting_mode, total_cells_per_rep, nreps, nreads_per_sample, r
        )


        if not os.path.exists(mageck_input_name) or rerun:
            _ = write_mageck_input(
                mageck_input_name, 
                total_cells_per_rep, 
                nreps, 
                nreads_per_sample, 
                sorting_mode = sorting_mode
            )
        
        mageck_res_name_prefix = "mageck_results/{}_{}_{}_{}_{}_{}".format(
            sorting_mode, measure, total_cells_per_rep, 
            nreps, nreads_per_sample, r
        )

        _, guide_info = _get_effect_sizes(total_cells = total_cells_per_rep)
        target_info = guide_info[['target_id', 'effect_size']].drop_duplicates()

        def get_samples_input(sort_bin):
            return(",".join(
            ["rep{}_{}_guide".format(i, sort_bin) for i in range(nreps)]))
        
        mageck_ctrl_samples = get_samples_input("bulk")
        for sorting_bin in sorting_bins:
            mageck_res_name = mageck_res_name_prefix + "_{}.gene_summary.txt".format(
                sorting_bin)
            mageck_treatment_samples = get_samples_input(sorting_bin)
            if not os.path.exists(mageck_res_name) or rerun:
                os.system("mageck test -k "+mageck_input_name+
                          " -t "+mageck_treatment_samples + 
                          " -c "+mageck_ctrl_samples+
                      " -n " + mageck_res_name_prefix + "_{}".format(sorting_bin) + 
                          " --paired")

def get_mageck_rra_sensitivity(sorting_mode = "topbot",
                               total_cells_per_rep = 10**5, 
                               nreps = 4, 
                               nreads_per_sample = 10**6,
                               measure = "guide", 
                               iteration_per_screen = 10, 
                               rerun = False, 
                              ):
    if sorting_mode == "topbot":
        sorting_bins = ["top", "bot"]
    elif sorting_mode == "bins":
        sorting_bins = ["top", "high", "low", "bot"]
        # read data
    sensitivities = []
    for sorting_bin in sorting_bins:
        mageck_res_name = mageck_res_name_prefix + "_{}.gene_summary.txt".format(
            sorting_bin)
        try:
            mageck_res = pd.read_csv(mageck_res_name, sep = "\t")
        except:
            print("ERROR in reading {}".format(mageck_res_name))
            return((None, None))
        mageck_res = target_info.merge(mageck_res, right_on = "id", left_on = "target_id")
        mageck_res["significant"] = (mageck_res["neg|fdr"] < 0.05) | (mageck_res["pos|fdr"] < 0.05)
        mean_sensitivity = mageck_res.groupby("effect_size").agg({"significant":"mean"})
        sensitivities.append(mean_sensitivity)

    mean_sensitivity = pd.concat(sensitivities, axis = 1)

    if r == 0:
        mean_sensitivities = mean_sensitivity
        mageck_reses = mageck_res
        mageck_reses["R"] = r
    else:
        mean_sensitivities = pd.concat((mean_sensitivities, mean_sensitivity), axis = 1)
        mageck_res["R"] = r
        mageck_reses = pd.concat((mageck_reses, mageck_res), axis = 0)


    plt.scatter(mageck_res["neg|rank"], mageck_res.effect_size, color= "grey")
    plt.xlabel("negative selection rank")
    plt.ylabel("effect size")
    plt.title("{}K cells, {}K reads, {} reps".format(
        total_cells_per_rep/1000, 
        nreads_per_sample/1000, 
        nreps))
    plt.scatter(mageck_res.loc[mageck_res.significant, "neg|rank"], 
                mageck_res.loc[mageck_res.significant, "effect_size"], 
                color = "red")
    plt.savefig("plots/ranks_{}_{}_{}_{}_{}.pdf".format(
        sorting_mode, measure, total_cells_per_rep, nreps, nreads_per_sample,))
    
    print("DONE running for {}_{}_{}_{}".format(
        sorting_mode, total_cells_per_rep, nreps, nreads_per_sample))

    return((mean_sensitivities, mageck_reses))      