using Plots

function simulate_LIF_neuron(A, d, f, tau_d, tau_f, dt, T, S_input)
    τ_m = 10.0       
    V_thresh = -50.0  
    V_reset = -75.0   
    V_rest = -75.0    
    R_m = 10.0        

    time = 0:dt:T      

    V = V_rest                
    Vs = Float64[]           
    spike_times = []         
    D = 1.0                  
    F = 1.0                  
    Hs = Float64[] 

    prev_V = V 
    spike_count = 0  

    for (idx, t) in enumerate(time)
        W = A * D * F
        dV = (-(V - V_rest) + R_m * W * S_input[idx]) / τ_m
        V += dV * dt

        if S_input[idx] == 1
            spike_count += 1
            if spike_count == 1 || spike_count == 5
                # Record the change in potential due to the input spike
                push!(Hs, V - prev_V)
            end
            D *= d  
            F += f  
        end

        dD = (1 - D) / tau_d
        D += dD * dt
    
        dF = (1 - F) / tau_f
        F += dF * dt
        
        prev_V = V
        push!(Vs, V)
    end
    
    return time, Vs, spike_times, Hs
end


function plot_LIF_simulation(time, Vs, spike_times, S_input, filename)
    p1 = plot(time, Vs, lw=2, label="Membrane Potential", xlabel="Time (ms)", ylabel="Membrane Potential (mV)", ylims=(-80, -40), legend=:topright, legendfontsize=8, title="LIF Neuron Spiking Activity", linecolor=:blue)
    scatter!(spike_times, fill(-50, length(spike_times)), markershape=:circle, color="red", label="Spikes", ms=4)
    
    p2 = plot(time, S_input, label="Spike Train", xlabel="Time (ms)", ylabel="Amplitude", color="purple", legend=:topright, ylims=(-0.1, 1.1), legendfontsize=8, linewidth=2)
    title!("Input Spike Train")

    p = plot(p1, p2, layout=(2,1), link=:x, size=(800, 600))

    savefig(p, filename)
end


function generate_spike_train(T, dt, initial_spike_time, tf)
    # T: total simulation time
    # dt: time step
    # initial_spike_time: time for the first spike
    # tf: temporal frequency, indicating how often a spike appears in the train
    
    # Calculate the number of time steps
    num_steps = convert(Int, T/dt) + 1
    
    # Initialize the spike train with all zeros
    S_input = zeros(Int, num_steps)
    
    # Set the initial spike
    S_input[convert(Int, initial_spike_time/dt)] = 1
    
    # Calculate the time interval between spikes based on the temporal frequency
    spike_interval = round(Int, 1000/tf)
    
    # Generate following spikes at evenly spaced intervals
    for i in 2:5  # since the first spike is already set, we generate the next 4 spikes
        next_spike_time = initial_spike_time + (i-1)*spike_interval
        if next_spike_time <= T  # only set spike if it is within the total time T
            S_input[convert(Int, next_spike_time/dt)] = 1
        end
    end
    
    return S_input
end

# Example usage
#A = 20.0        # fixed parameter A
#d = 0.23       # depression fraction upon a spike
#f = 0.083         # facilitation increment upon a spike
#tau_d = 52.0     # time constant for D to recover to 1 (ms)
#tau_f = 96.0     # time constant for F to recover to 1 (ms)
#dt = 1.0        # time step (ms)
#T = 3000.0       # total time to simulate (ms)
#Old parameters above

# Example usage
A = 20.0        # fixed parameter A
d = 0.15       # depression fraction upon a spike
f = 0.92         # facilitation increment upon a spike
tau_d = 103.0     # time constant for D to recover to 1 (ms)
tau_f = 96.0     # time constant for F to recover to 1 (ms)
dt = 1.0        # time step (ms)
T = 3000.0       # total time to simulate (ms)

# Generate spike train
first_spike_time = 50.0  # ms
temporal_frequency = 2.0  # Hz
S_input = generate_spike_train(T, dt, first_spike_time, temporal_frequency)

# Run simulation
time, Vs, spike_times, Hs = simulate_LIF_neuron(A, d, f, tau_d, tau_f, dt, T, S_input)

# Calculate depression ratio H_5/H_1
depression_ratio = Hs[2] / Hs[1]
println("Depression ratio H_5/H_1: ", depression_ratio)


# Plot results and save to file
output_filename = "./figs/working.png"
plot_LIF_simulation(time, Vs, spike_times, S_input, output_filename)

using Plots

# Example temporal_frequencies
temporal_frequencies = [2,3,4,5,6,7,8,9,10,11,12, 13, 15, 20, 25, 30,35, 40,45,50,60,70,80]
#temporal_frequencies = [2,3,4,5,6,7,8,9,10,11,12, 13, 15, 20, 25, 30,35, 40]
# Initialize an empty array to store results
depression_ratios = []

# Loop over all temporal_frequencies
for tf in temporal_frequencies
    global S_input, time, Vs, spike_times, Hs, depression_ratio
    # Example: Generate spike train using the previously defined function
    S_input = generate_spike_train(T, dt, first_spike_time, tf)

    # Example: Run simulation and store depression ratio
    time, Vs, spike_times, Hs = simulate_LIF_neuron(A, d, f, tau_d, tau_f, dt, T, S_input)
    
    # Ensure that we have recorded the required Hs values to calculate the depression ratio
    if length(Hs) >= 2
        depression_ratio = Hs[2] / Hs[1]
        push!(depression_ratios, depression_ratio)
    else
        println("Warning: Not enough Hs values recorded for tf = $tf. Skipping.")
        continue  # skip to the next iteration of the loop
    end
end

# Plot results
plot(temporal_frequencies, depression_ratios, xlabel="Temporal Frequency (Hz)", ylabel="Depression Ratio", label="Depression Ratio vs Temporal Frequency", marker=:circle, linewidth=2, linecolor=:blue)

# Save figure
savefig("depression_ratio_vs_temporal_frequency.png")