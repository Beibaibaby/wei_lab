using Plots
include("LIF_1.0.0_compare.jl")

function main_simulation(d_ee, f_ee, d_ie, f_ie, output_filename)
    A = 20.0
    tau_d = 103.0
    tau_f = 96.0
    dt = 1.0
    T = 3000.0
    first_spike_time = 50.0
    temporal_frequencies = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 20, 25, 30, 40, 50, 60, 75, 80, 100, 110, 130, 150]

    depression_ratios = []
    depression_ratios_fixed_F = []

    for tf in temporal_frequencies
        S_input = generate_spike_train(T, dt, first_spike_time, tf)

        # Simulate with variable F
        _, _, _, Hs = simulate_LIF_neuron(A, d_ee, f_ee, tau_d, tau_f, dt, T, S_input)
        if length(Hs) >= 2
            push!(depression_ratios, Hs[2] / Hs[1])
        end

        # Simulate with fixed F
        _, _, _, Hs_fixed_F = simulate_LIF_neuron_fixed_F(A, d_ie, f_ie, tau_d, tau_f, dt, T, S_input)
        if length(Hs_fixed_F) >= 2
            push!(depression_ratios_fixed_F, Hs_fixed_F[2] / Hs_fixed_F[1])
        end
    end

    # Plot and save results
    p = plot(temporal_frequencies, depression_ratios, xlabel="Temporal Frequency (Hz)", title="d_ee=$d_ee f_ee=$f_ee d_ie=$d_ie f_ie=$f_ie", ylabel="Depression Ratio", label="E->E", marker=:circle, markercolor=:red, linewidth=2, linecolor=:red)
    plot!(p, temporal_frequencies, depression_ratios_fixed_F, label="E->I", marker=:circle, linewidth=2, linecolor=:deepskyblue2, markercolor=:deepskyblue2)
    savefig(p, output_filename)
end

main_simulation(0.1, 0.1, 0.1, 0.55, "depression_ratio_vs_temporal_frequency_compare.png")