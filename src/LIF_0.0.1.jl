using Plots

function simulate_LIF_neuron(A, d, f, tau_d, tau_f, dt, T, S_input)
    # A is fixed parameter A
    # d is depression fraction upon a spike
    #f is facilitation increment upon a spike
    # tau_d is time constant for D to recover to 1 (ms)
    # tau_f is time constant for F to recover to 1 (ms)
    # dt is time step (ms)
    # T is total time to simulate (ms)

    # LIF neuron parameters
    τ_m = 10.0       
    V_thresh = -50.0  
    V_reset = -75.0   
    V_rest = -75.0    
    R_m = 10.0         
    
    # Simulation parameters
    time = 0:dt:T      
    
    # Initialize variables
    V = V_rest               
    Vs = Float64[]           
    spike_times = []         
    D = 1.0                  
    F = 1.0                  
    
    # Ensure S_input is the right length
    if length(S_input) != length(time)
        error("Length of S_input must be equal to length of time vector")
    end

    # Main simulation loop
    for (idx, t) in enumerate(time)
        # Compute the synaptic weight
        W = A * D * F
        
        # LIF equation with synaptic weight
        dV = (-(V - V_rest) + R_m * W * S_input[idx]) / τ_m
        V += dV * dt
    
        # Check for spike
        if V >= V_thresh
            push!(spike_times, t)
            V = V_reset
        end
        
        # Update D and F on spike occurrence in input spike train
        if S_input[idx] == 1
            D *= d  
            F += f  
        end
    
        # Recover D and F to 1 with exponential dynamics
        dD = (1 - D) / tau_d
        D += dD * dt
    
        dF = (1 - F) / tau_f
        F += dF * dt
        
        push!(Vs, V)
    end
    
    return time, Vs, spike_times
end


# Plotting function
function plot_LIF_simulation(time, Vs, spike_times, S_input, filename)
    
    p1 = plot(time, Vs, lw=2, label="Membrane Potential", xlabel="Time (ms)", ylabel="Membrane Potential (mV)", ylims=(-80, -40), legend=:topright, legendfontsize=8, title="LIF Neuron Spiking Activity", linecolor=:blue)
    scatter!(spike_times, fill(-50, length(spike_times)), markershape=:circle, color="red", label="Spikes", ms=4)
    
    p2 = plot(time, S_input, label="Spike Train", xlabel="Time (ms)", ylabel="Amplitude", color="purple", legend=:topright, ylims=(-0.1, 1.1), legendfontsize=8, linewidth=2)
    title!("Input Spike Train")
    
    # Combine the plots vertically
    p = plot(p1, p2, layout=(2,1), link=:x, size=(800, 600))
    
    # Save figure to file
    savefig(p, filename)
end

# Example usage
A = 20.0        # fixed parameter A
d = 0.5         # depression fraction upon a spike
f = 0.2         # facilitation increment upon a spike
tau_d = 5.0     # time constant for D to recover to 1 (ms)
tau_f = 5.0    # time constant for F to recover to 1 (ms)
dt = 1.0        # time step (ms)
T = 600.0       # total time to simulate (ms)

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
    spike_interval = convert(Int, 1000/tf)
    
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
A = 20.0        # fixed parameter A
d = 0.5         # depression fraction upon a spike
f = 0.2         # facilitation increment upon a spike
tau_d = 5.0     # time constant for D to recover to 1 (ms)
tau_f = 5.0     # time constant for F to recover to 1 (ms)
dt = 1.0        # time step (ms)
T = 600.0       # total time to simulate (ms)

# Generate spike train
initial_spike_time = 50  # time for the first spike (ms)
tf = 10.0  # temporal frequency (Hz)
S_input = generate_spike_train(T, dt, initial_spike_time, tf)

# Run simulation
time, Vs, spike_times = simulate_LIF_neuron(A, d, f, tau_d, tau_f, dt, T, S_input)

# Plot results and save to file
output_filename = "./figs/working.png"
plot_LIF_simulation(time, Vs, spike_times, S_input, output_filename)
