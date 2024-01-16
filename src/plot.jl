using JLD2  # For loading the data
using Plots         # For plotting
using ProgressMeter # For the progress bar
using Measures

# Load the data
@load "times.jld2" times

# Define your parameters
params = Dict(
    :Ne => 4000,  # The number of excitatory cells, replace with actual value
    :Ncells => 5000,  # Total number of cells, replace with actual value
    :T => 3000
)

function compute_sliding_rate(spiketimes, window_size, step_size, T)
    println("starting computing firing rate")
    n_neurons, _ = size(spiketimes)
    n_steps = floor(Int, (T - window_size) / step_size) + 1
    rates = zeros(n_steps)

    p = Progress(n_steps, 2) # Initialize the progress bar

    for i = 1:n_steps
        next!(p) # Update the progress bar

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


# Parameters for sliding window
window_size = 100  # in ms
step_size = 20     # in ms

# Compute the firing rates
e_rate = compute_sliding_rate(times[1:params[:Ne], :], window_size, step_size, params[:T])
i_rate = compute_sliding_rate(times[(params[:Ne]+1):params[:Ncells], :], window_size, step_size, params[:T])

# Compute the time values
n_steps = length(e_rate)  # Assuming e_rate and i_rate have the same length
time_values = [(i-1) * step_size + window_size for i in 1:n_steps]

# Plot settings
plot_size = (1000, 600)
plot_margin = 5mm

# Plot the firing rates
p2 = plot(time_values, e_rate, xlabel="Time (ms)", ylabel="Firing rate (Hz)",
    label="Excitatory", lw=2, linecolor=:red, size=plot_size,
    title="Firing rate", ylim=(0,5))

plot!(p2, time_values, i_rate, label="Inhibitory", lw=2, linecolor=:deepskyblue2,
    left_margin=plot_margin)

# Save the plot
savefig(p2, "firing_rate_plot.png")  # Save as PDF