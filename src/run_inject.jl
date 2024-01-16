using Statistics
using Plots  
using Dates  # for generating timestamps
using JSON3  # Use JSON3 package for JSON handling
using Random


struct NetworkParameters
    K::Int
    Ne::Int
    Ni::Int
    T::Int
    sqrtK::Float64
    taue::Int
    taui::Int
    jie::Float64
    jee::Float64
    jei::Float64
    jii::Float64
end

doplot = true

include("network.jl")
include("sim_inject.jl")
K = 800
Ne = 4000 #Ne is the Number of E cells
Ni = 1000 #Ni is the Number of I cells
T =  1000 #T is the simulation time (ms)
sqrtK = sqrt(K) 
taue = 10 #membrane time constant for exc. neurons (ms)
taui = 15 #membrane time constant for inh. neurons (ms)

jie = 4. / (taui*sqrtK)        #strength from E to I
jee = 10. / (taue*sqrtK)        #strength from E to E 

jei = -16. * 1.2 /(taue*sqrtK)      #strength from I to E
jii = -16. / (taui*sqrtK)      #strength from I to I 

#times, ns, Ne, Ncells, T = sim(K,Ne,Ni,T,jie,jei,jii,jee)



#params = NetworkParameters(800, 4000, 1000, 1500, sqrt(800), 10, 15, 0.2 / (15 * sqrt(800)), 0.4 / (10 * sqrt(800)), -0.8 * 1.2 / (10 * sqrt(800)), -0.8 / (15 * sqrt(800)))
params = NetworkParameters(K, Ne, Ni, T, sqrtK, taue, taui, jie, jee, jei, jii)
#times, ns, Ne, Ncells, T = sim(params.K, params.Ne, params.Ni, params.T, params.jie, params.jei, params.jii, params.jee)
times, ns, Ne, Ncells, T, v_history, E_input, I_input, weights = sim_old()

EPSP_EPSP_pool=v_history[1:100,end-9999:end]
TT_pool=v_history[501:600,end-9999:end]
IPSP_IPSP_pool=v_history[1001:1100,end-9999:end]

IPSP_IPSP_pool= IPSP_IPSP_pool .- 1
TT_pool=TT_pool .- 1
EPSP_EPSP_pool=EPSP_EPSP_pool .+ 0.2

# Initialize accumulators
EPSP_EPSP_accumulator = zeros(100, 10000)
TT_accumulator = zeros(100, 10000)
IPSP_IPSP_accumulator = zeros(100, 10000)

for _ in 1:50
    # Assuming you've defined the NetworkParameters somewhere above
    # params = NetworkParameters(K, Ne, Ni, T, sqrtK, taue, taui, jie, jee, jei, jii)
    times, ns, Ne, Ncells, T, v_history, E_input, I_input, weights = sim_old()
    
    EPSP_EPSP_pool = v_history[1:100, end-9999:end]
    TT_pool = v_history[501:600, end-9999:end]
    IPSP_IPSP_pool = v_history[1001:1100, end-9999:end]
    
    IPSP_IPSP_pool = IPSP_IPSP_pool .- 2
    TT_pool = TT_pool .- 1
    EPSP_EPSP_pool = EPSP_EPSP_pool .+ 0.2
    
    # Add to accumulator
    EPSP_EPSP_accumulator .+= EPSP_EPSP_pool
    TT_accumulator .+= TT_pool
    IPSP_IPSP_accumulator .+= IPSP_IPSP_pool
end
# Compute the average
EPSP_EPSP_avg = EPSP_EPSP_accumulator ./ 50
TT_avg = TT_accumulator ./ 50
IPSP_IPSP_avg = IPSP_IPSP_accumulator ./ 50

IPSP_IPSP_pool= EPSP_EPSP_avg
TT_pool=TT_avg
EPSP_EPSP_pool=IPSP_IPSP_avg

# Run simulation 10 times



#indices = find_positive_weights_indices(weights, 4000)
#println("indices ", indices)

println("mean excitatory firing rate: ", mean(1000 * ns[1:params.Ne] / params.T), " Hz")
println("mean inhibitory firing rate: ", mean(1000 * ns[(params.Ne+1):Ncells] / params.T), " Hz")

function compute_sliding_rate(spiketimes, window_size, step_size, T)
    n_neurons, _ = size(spiketimes)
    n_steps = floor(Int, (T - window_size) / step_size) + 1
    rates = zeros(n_steps)

    for i = 1:n_steps
        t_start = (i-1) * step_size + 1  # Ensure the start time is non-zero
        t_end = t_start + window_size - 1  # Adjust end time based on start
        
        # Check for out-of-bounds
        if t_end > T
            break
        end

        n_spikes = sum([sum((t_start .<= spiketimes[j, :]) .& (spiketimes[j, :] .< t_end)) for j=1:n_neurons])
        rates[i] = n_spikes / (n_neurons * window_size) * 1000  # rate in Hz
    end
    #println(rates)
    return rates
end

if doplot 

        # Generate a timestamp for unique filenames
        timestamp = Dates.now()
        timestamp_str = Dates.format(timestamp, "yyyy-mm-dd_HH-MM-SS")

    println("creating plot")

    # Define plot size: (width, height)
    plot_size = (600, 400) 

 # Parameters for sliding window
   window_size = 100  # in ms
   step_size = 5     # in ms

   e_rate = compute_sliding_rate(times[1:Ne, :], window_size, step_size, T)
   i_rate = compute_sliding_rate(times[(Ne+1):Ncells, :], window_size, step_size, T)

# Compute the time values based on window_size and step_size
   n_steps = length(e_rate)  # or length(i_rate), assuming they have the same length
   time_values = [i * step_size + (window_size / 2) for i in 1:n_steps]

   p2 = plot(time_values, e_rate, xlabel="Time (ms)", ylabel="Firing rate (Hz)",
          label="Excitatory", lw=2, linecolor=:red, size=plot_size)
plot!(time_values, i_rate, label="Inhibitory", lw=2, linecolor=:deepskyblue2)

   fig_filename = "./figs/$timestamp_str.png"

# Save the figure with the timestamped filename
   savefig(p2, fig_filename)

println("Figure saved as $fig_filename")
    


    # Define filenames for the figure and JSON file
    
json_filename = "./figs/$timestamp_str.json"

param_dict = Dict(
    "K" => params.K,
    "Ne" => params.Ne,
    "Ni" => params.Ni,
    "T" => params.T,
    "sqrtK" => params.sqrtK,
    "taue" => params.taue,
    "taui" => params.taui,
    "jie" => params.jie,
    "jee" => params.jee,
    "jei" => params.jei,
    "jii" => params.jii
)

   JSON3.open(json_filename, "w") do io
    JSON3.print(io, param_dict)
end

println("Parameters saved as $json_filename")

avg_correlation_E_I=compute_correlation(E_input, I_input)
avg_correlation_E_E=compute_correlation(E_input, E_input)
avg_correlation_I_I=compute_correlation(I_input, I_input)
cross_corr_E_E=compute_cross_correlation(E_input[:, end-999:end], E_input[:, end-999:end])
cross_corr_I_I=compute_cross_correlation(I_input[:, end-999:end], I_input[:, end-999:end])
cross_corr_E_I=compute_cross_correlation(E_input[:, end-999:end], I_input[:, end-999:end])
cross_corr_I_E=compute_cross_correlation(I_input[:, end-999:end], E_input[:, end-999:end])


println("avg correlation(E-I): ", avg_correlation_E_I)
println("avg correlation(E-E): ", avg_correlation_E_E)
println("avg correlation(I-I): ", avg_correlation_I_I)

cross_EPSP_EPSP=compute_cross_correlation(EPSP_EPSP_pool,EPSP_EPSP_pool)
cross_IPSP_IPSP=compute_cross_correlation(IPSP_IPSP_pool,IPSP_IPSP_pool)
cross_EPSP_IPSP=compute_cross_correlation(EPSP_EPSP_pool,IPSP_IPSP_pool)
cross_IPSP_EPSP=compute_cross_correlation(IPSP_IPSP_pool,EPSP_EPSP_pool)
cross_T_T=compute_cross_correlation(TT_pool,TT_pool)

#println(cross_EPSP_EPSP)
#println("cross E-E: ", cross_corr_E_E)


function plot_correlations(cross_corr_E_E, cross_corr_I_I, cross_corr_E_I,cross_corr_I_E)
    # Assuming the dictionaries have the same keys
    sorted_keys = sort(collect(keys(cross_corr_E_E)))
    
    sorted_values_E_E = [cross_corr_E_E[k] for k in sorted_keys]
    sorted_values_I_I = [cross_corr_I_I[k] for k in sorted_keys]
    sorted_values_E_I = [cross_corr_E_I[k] for k in sorted_keys]
    sorted_values_I_E = [cross_corr_I_E[k] for k in sorted_keys]

    plot(sorted_keys, sorted_values_E_E, label="E-E", linewidth=2, marker=:circle, color=:blue, legend=:topright)
    plot!(sorted_keys, sorted_values_I_I, label="I-I", linewidth=2, marker=:circle, color=:red)
    plot!(sorted_keys, sorted_values_E_I, label="E-I", linewidth=2, marker=:circle, color=:yellow)
    plot!(sorted_keys, sorted_values_I_E, label="I-E", linewidth=2, marker=:circle, color=:orange)
    
    xlabel!("Tau")
    ylabel!("Correlation")
    title!("Cross-correlation")
    savefig("cross_correlation_plot_PSP_new.png")
    display(plot)
end


function plot_correlations_mem(cross_corr_E_E, cross_corr_I_I, cross_corr_T_T, cross_corr_E_I, cross_corr_I_E)
    # Assuming the dictionaries have the same keys
    sorted_keys = sort(collect(keys(cross_corr_E_E)))
    
    sorted_values_E_E = [cross_corr_E_E[k] for k in sorted_keys]
    sorted_values_I_I = [cross_corr_I_I[k] for k in sorted_keys]
    sorted_values_T_T = [cross_corr_T_T[k] for k in sorted_keys]
    sorted_values_E_I = [cross_corr_E_I[k] for k in sorted_keys]
    sorted_values_I_E = [cross_corr_I_E[k] for k in sorted_keys]

    plot(sorted_keys, sorted_values_E_E, label="E-E", linewidth=2, marker=:circle, color=:blue, legend=:topright)
    plot!(sorted_keys, sorted_values_I_I, label="I-I", linewidth=2, marker=:circle, color=:red)
    plot!(sorted_keys, sorted_values_E_I, label="E-I", linewidth=2, marker=:circle, color=:yellow)
    plot!(sorted_keys, sorted_values_I_E, label="I-E", linewidth=2, marker=:circle, color=:orange)
    plot!(sorted_keys, sorted_values_T_T, label="T-T", linewidth=2, marker=:circle, color=:green)

    xlabel!("Tau")
    ylabel!("Correlation")
    title!("Cross-correlation")
    savefig("cross_correlation_plot_PSP_new_high.png")
    display(plot)
end

#plot_correlations(cross_corr_E_E, cross_corr_I_I, cross_corr_E_I,cross_corr_I_E)
plot_correlations_mem(cross_EPSP_EPSP, cross_IPSP_IPSP, cross_T_T, cross_EPSP_IPSP, cross_IPSP_EPSP)


function plot_cells(v_history, cells)
    p = plot(layout=(3,1), size=(600,800))

    for (index, cell) in enumerate(cells)
        plot!(p[index], v_history[cell, :], label="Cell $cell", xlabel="Time", ylabel="Voltage", legend=:topright)
    end
    savefig("v_history.png")
    display(p)
end

# Example usage: Plotting cells 1, 2, and 3
plot_cells(v_history, [1, 505, 1005])

println("finish")


end