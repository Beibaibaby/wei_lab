using JLD2

# Path to one of your JLD2 files
file_path = "/gpfs/data/doiron-lab/draco/results_corr/exp_20231114_233353/1/cross_corr_E_E.jld2"

jldopen(file_path, "r") do file
    println(keys(file))
end