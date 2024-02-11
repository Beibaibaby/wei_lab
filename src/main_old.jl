using Plots
using Measures
using Optim

"""
Simulate the dynamics of a Leaky Integrate-and-Fire (LIF) neuron with synaptic depression, facilitation, and dual exponential synaptic input model.

# Arguments
- `A`: Amplitude of synaptic input.
- `d_1`, `d_2`: Depression factors.
- `f`: Facilitation increment.
- `tau_d_1`, `tau_d_2`: Time constants for recovery of depression.
- `tau_f`: Time constant for recovery of facilitation.
- `dt`: Time step for the simulation.
- `T`: Total simulation time.
- `S_input`: Binary array representing the presence of spikes over time.
- `taurise`, `taudecay`: Time constants for the rise and decay of the synaptic input.

# Returns
Tuple of time array, membrane potentials, spike times, H values, D values, and F values.
"""
taurise = 1
taudecay = 3


function simulate_LIF_neuron(A, d_1, d_2, f, tau_d_1, tau_d_2, tau_f, dt, T, S_input, taurise, taudecay)
    # Neuron parameters
    τ_m = 10.0        # Membrane time constant
    V_thresh = -50.0  # Voltage threshold for spike
    V_reset = -75.0   # Reset voltage post-spike
    V_rest = -75.0    # Resting membrane potential
    R_m = 10.0        # Membrane resistance

    # Time vector
    time = 0:dt:T

    # Initialize variables
    V = V_rest                # Membrane potential
    Vs = Float64[]            # Array to store membrane potential over time
    spike_times = []          # Spike times
    D_1 = 1.0                 # Depression factor 1
    D_2 = 1.0                 # Depression factor 2
    F = 1.0                   # Facilitation factor
    Hs = Float64[]            # Array to store H values
    Ds = Float64[]            # Array to store D values
    Fs = Float64[]            # Array to store F values
    prev_V = V                # Previous voltage
    spike_count = 0           # Spike counter

    # Synaptic input variables
    xrise = zeros(length(time))
    xdecay = zeros(length(time))

    is_rising = false  # Flag to track whether synInput is rising
    prev_synInput = 0.0  # Variable to store previous synaptic input

    for (idx, t) in enumerate(time)
        # Update synaptic input variables based on input spikes
        if S_input[idx] == 1
            xrise[idx] += 1
            xdecay[idx] += 1
        end

        # Calculate the synaptic input using the dual exponential model
        if idx > 1
            xrise[idx] += -dt * xrise[idx - 1] / taurise
            xdecay[idx] += -dt * xdecay[idx - 1] / taudecay
        end


        synInput = (xdecay[idx] - xrise[idx]) / (taudecay - taurise)

        # Compute membrane potential
        W = A * D_1 * D_2 * F
        dV = (-(V - V_rest) + R_m * W * synInput) / τ_m
        V += dV * dt

        # Update depression and facilitation factors
        if S_input[idx] == 1
            D_1 *= d_1
            D_2 *= d_2
            F += f  
        end

        # Check for rising edge of synInput
        if synInput > prev_synInput
            is_rising = true
        elseif is_rising && synInput < prev_synInput
            # synInput starts to drop, indicating a peak
            spike_count += 1
            if spike_count == 1 || spike_count == 2
                push!(Hs, V - V_rest) # Change in potential due to the input spike
            end
            is_rising = false
        end

        prev_synInput = synInput


        dD_1 = (1 - D_1) / tau_d_1
        D_1 += dD_1 * dt

        dD_2 = (1 - D_2) / tau_d_2
        D_2 += dD_2 * dt

        dF = (1 - F) / tau_f
        F += dF * dt

        # Store computed values
        prev_V = V
        push!(Vs, V)
        push!(Ds, D_1 * D_2)
        push!(Fs, F)
    end
    
    return time, Vs, spike_times, Hs, Ds, Fs, xrise, xdecay
end




function plot_LIF_simulation(time, Vs, spike_times, S_input, Ds, Fs, filename)
    p1 = plot(time[1:200], Vs[1:200], lw=2, label="Membrane Potential", ylabel="Membrane Potential (mV)", ylims=(-80, -50), legend=:topright, legendfontsize=8, linecolor=:orange,grid=false,dpi=400,left_margin=15mm,linewidth=10)
    #scatter!(spike_times, fill(-50, length(spike_times)), markershape=:circle, color="red", label="Spikes", ms=4)#mediumpurple3
    
    #p2 = plot(time, S_input, label="Spike Train", xlabel="Time (ms)", ylabel="Amplitude", color="purple", legend=:topright, ylims=(-0.1, 1.1), legendfontsize=8, linewidth=2)
    #title!("Input Spike Train")

    p3 = plot(time[1:200], Ds[1:200], lw=2, label="Depression Factor D", xlabel="Time (ms)", ylabel="Value", color=:black, legend=:topright, legendfontsize=8, linewidth=10,grid=false,dpi=400,left_margin=15mm,ylim=(0,5))
    plot!(p3, time[1:200], Fs[1:200], lw=2, label="Facilitation Factor F", color=:darkgreen, legend=:topright, legendfontsize=8, linewidth=10)
    #plot!(p3, time, Fs .* Ds, lw=2, label="F times D", color=:black, legend=:topright, legendfontsize=8, linewidth=2)
    #title!("D and F Over Time")

    p = plot(p1, p3, layout=(2,1), link=:x, size=(1000, 550),dpi=400)  # Adjusted the layout and size

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

#

# Example usage
A = 20.0        # fixed parameter A
d_1 = 0.24       # depression fraction upon a spike
d_2 = 0.24
f = 0.00        # facilitation increment upon a spike
tau_d_1 = 103     # time constant for D to recover to 1 (ms)
tau_d_2 = 103
tau_f = 96.0     # time constant for F to recover to 1 (ms)
dt = 0.1       # time step (ms)
T = 2000.0       # total time to simulate (ms)

# Generate spike train
first_spike_time = 50.0  # ms
#temporal_frequency = 60.0  # Hz
#S_input = generate_spike_train(T, dt, first_spike_time, temporal_frequency)

# Run simulation
#time, Vs, spike_times, Hs, Ds, Fs = simulate_LIF_neuron(A, d, f, tau_d, tau_f, dt, T, S_input) 

# Calculate depression ratio H_5/H_1
#depression_ratio = Hs[2] / Hs[1]
#println("Depression ratio H_5/H_1: ", depression_ratio)


# Plot results and save to file
#output_filename = "./tf=60_F.png"
#plot_LIF_simulation(time, Vs, spike_times, S_input, Ds, Fs, output_filename)  # Added Ds and Fs


function generate_two_spike_train(T, dt, initial_spike_time, interval)
    # Initialize the spike train with all zeros
    num_steps = convert(Int, T/dt) + 1
    S_input = zeros(Int, num_steps)
    
    # Set the first spike
    S_input[convert(Int, initial_spike_time/dt)] = 1
    
    # Set the second spike based on the interval
    second_spike_time = initial_spike_time + interval
    if second_spike_time <= T
        S_input[convert(Int, second_spike_time/dt)] = 1
    end
    
    return S_input
end



#stimulate LIF with two small d, write a function 

# Observed data
observed_intervals = [50, 100, 200, 400, 800, 1600]
observed_pprs = [0.717, 0.728, 0.800, 0.835,  0.899, 1.010]

# Updated cost function to accept a vector of all parameters
function cost_function(params)
    A, d_1,d_2, f, tau_d_1,tau_d_2, tau_f = params  # Unpack all parameters

    simulated_pprs = Float64[]

    for interval in observed_intervals
        local S_input = generate_two_spike_train(T, dt, first_spike_time, interval)
        time, Vs, spike_times, Hs, Ds, Fs,xrise, xdecay = simulate_LIF_neuron(A, d_1, d_2, f, tau_d_1,tau_d_2, tau_f, dt, T, S_input,taurise, taudecay)
        ppr = Hs[2] / Hs[1]  # Compute PPR
        push!(simulated_pprs, ppr)
    end

    # Calculate the sum of squared differences (mean squared error)
    return sum((simulated_pprs - observed_pprs).^2) / length(observed_pprs)
end



# Set initial guesses for all parameters
initial_params = [20.0, d_1,d_2, 0.000000001, 103.0,103.0, 96.0]  # A, d_1, d_2 f, tau_d, tau_f

# Define lower and upper bounds for each parameter
# Assuming all parameters are positive and have an upper limit for the example
lower_bounds = [0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0]
upper_bounds = [100.0, 1.0, 1.0,1.0, 500.0,500.0, 500.0]

# Define the box constraints
box = Fminbox(BFGS())

# Run the optimization with constraints
result = optimize(cost_function, lower_bounds, upper_bounds, initial_params, box, Optim.Options(g_tol = 1e-6))

# Extract optimized values
optimized_params = Optim.minimizer(result)
#rewriting the optimized parameters
optimized_A, optimized_d_1,optimized_d_2, optimized_f, optimized_tau_d_1,optimized_tau_d_2, optimized_tau_f = optimized_params
# Print the optimized values
println("Optimized A: $optimized_A, d_1: $optimized_d_1,d_2: $optimized_d_2, f: $optimized_f, tau_d_1: $optimized_tau_d_1,tau_d_2: $optimized_tau_d_2, tau_f: $optimized_tau_f")

# Use optimized




# Simulate PPRs using optimized parameters but change the agrumants to the function

function simulate_pprs(A, d_1,d_2, f, tau_d_1,tau_d_2, tau_f,intervals)
    simulated_pprs = Float64[]
    for interval in intervals
        S_input = generate_two_spike_train(T, dt, first_spike_time, interval)
        time, Vs, spike_times, Hs, Ds, Fs,xrise, xdecay = simulate_LIF_neuron(A, d_1, d_2, f, tau_d_1,tau_d_2, tau_f, dt, T, S_input,taurise, taudecay)
        print(Hs[2])
        print(Hs[1])
        ppr = Hs[2] / Hs[1]
        
        push!(simulated_pprs, ppr)
    end
    return simulated_pprs
end

simulated_pprs = simulate_pprs(optimized_A, optimized_d_1,optimized_d_2, optimized_f, optimized_tau_d_1,optimized_tau_d_2, optimized_tau_f, observed_intervals)

# Sort the observed data
sorted_observed_indices = sortperm(observed_intervals)
sorted_observed_intervals = observed_intervals[sorted_observed_indices]
sorted_observed_pprs = observed_pprs[sorted_observed_indices]

# Sort the simulated data
sorted_simulated_indices = sortperm(observed_intervals)
sorted_simulated_intervals = observed_intervals[sorted_simulated_indices]
sorted_simulated_pprs = simulated_pprs[sorted_simulated_indices]

# Plotting

# Prepare the title with optimized parameters

title_str = "A= $(round(optimized_A, digits=2)), d_1= $(round(optimized_d_1, digits=2)), d_2= $(round(optimized_d_2, digits=2)), f= $(round(optimized_f, digits=2)), tau_d_1= $(round(optimized_tau_d_1, digits=2)), tau_d_2= $(round(optimized_tau_d_2, digits=2)), tau_f= $(round(optimized_tau_f, digits=2))"

# Plotting with the optimized parameters in the title
# change the frontsize of the title to become smaller (title=title_str)
p = plot(title=title_str, legend=:bottomright, xlabel="Interval (ms)", ylabel="Paired Pulse Ratio",
         legendfontsize=8, grid=false, dpi=400, left_margin=15mm, linewidth=10, titlefontsize=7)

# Plot observed data as a curve with dots
plot!(p, sorted_observed_intervals, sorted_observed_pprs, label="Observed Data", color=:blue, line=:solid)
scatter!(sorted_observed_intervals, sorted_observed_pprs, color=:blue, markersize=4, label=false)

# Plot simulated data as a curve with dots
plot!(p, sorted_simulated_intervals, sorted_simulated_pprs, label="Simulated Data", color=:red, line=:solid)
scatter!(sorted_simulated_intervals, sorted_simulated_pprs, color=:red, markersize=4, label=false)

# Save and display the plot
savefig(p, "PPR_Comparison_Plot.png")
#display(p)





# Generate a spike train with 5 spikes
first_spike_time = 50.0  # ms
temporal_frequency = 10.0  # Hz (adjust as needed for 5 spikes within T)
S_input = generate_spike_train(T, dt, first_spike_time, temporal_frequency)

# Run the simulationoptimized_A, optimized_d_1,optimized_d_2, optimized_f, optimized_tau_d_1,optimized_tau_d_2, optimized_tau_f, observed_intervals

#simulate_LIF_neuron(A, d_1, d_2, f, tau_d_1, tau_d_2, tau_f, dt, T, S_input, taurise, taudecay)

time, Vs, spike_times, Hs, Ds, Fs, xrise, xdecay = simulate_LIF_neuron(optimized_A, optimized_d_1,optimized_d_2, optimized_f, optimized_tau_d_1,optimized_tau_d_2, optimized_tau_f, dt, T, S_input, taurise, taudecay)

# Calculate synaptic input over time
synInputs = [(xdecay - xrise) / (taudecay - taurise) for (xrise, xdecay) in zip(xrise, xdecay)]

# Plot membrane potential and synaptic input
# Plot membrane potential
p1 = plot(time, Vs, label="Membrane Potential (V)", color=:blue, xlabel="Time (ms)", ylabel="Membrane Potential (mV)", left_margin=20mm)

# Plot synaptic input
p2 = plot(time, synInputs, label="Synaptic Input", color=:red, xlabel="Time (ms)", ylabel="Synaptic Input", left_margin=20mm)

# Plot H values
p3 = plot(1:length(Hs), Hs, label="H Values", color=:green, xlabel="Spike Number", ylabel="H Value", left_margin=20mm, legend=:topright)

# Combine all plots into a single layout
combined_plot = plot(p1, p2, p3, layout=(3,1), legend=:topright)

savefig(combined_plot, "Neuron_Simulation.png")