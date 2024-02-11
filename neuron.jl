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
            println("Beep! Input spike")
            println(t)
            println(idx)
            xrise[idx] += 1
            xdecay[idx] += 1
        end

        # Calculate the synaptic input using the dual exponential model
        if idx >1
            xrise[idx] += -dt * xrise[idx-1] / (taurise) +xrise[idx-1]
            xdecay[idx] += -dt * xdecay[idx-1] / (taudecay)+ xdecay[idx-1]
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
        
        if idx ==500 || idx == 501 || idx == 502 || idx == 503 || idx == 504 || idx == 499
           
            println(t)
            println(idx)
            println("Synaptic Input: ", synInput)
            println("change of xrise", dt * xrise[idx - 1] / (taurise) )
            println("change of xdecay", dt * xdecay[idx - 1] / (taudecay))
            println(xrise[idx])
            println(xdecay[idx])
        end


        # Check for rising edge of synInput
        if synInput > prev_synInput #check if input is rising
            is_rising = true
            #println("Rising") 
            #println(idx)
        elseif is_rising && synInput <= prev_synInput && synInput>0
            # synInput starts to drop, indicating a peak
            println(idx)

            spike_count += 1
            #if spike_count == 1 || spike_count == 2
            push!(Hs, synInput) # Change in potential due to the input spike
            #end
            
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

# Example usage
A = 20.0        # fixed parameter A
d_1 = 0.20       # depression fraction upon a spike
d_2 = 0.20
f = 0.00        # facilitation increment upon a spike
tau_d_1 = 103     # time constant for D to recover to 1 (ms)
tau_d_2 = 103
tau_f = 96.0     # time constant for F to recover to 1 (ms)
dt = 0.1       # time step (ms)
T = 1000.0       # total time to simulate (ms)
taurise=0.95
taudecay=1.0

# Generate a spike train with 5 spikes
first_spike_time = 50.0  # ms
temporal_frequency = 10.0  # Hz (adjust as needed for 5 spikes within T)
S_input = generate_spike_train(T, dt, first_spike_time, temporal_frequency)

# Run the simulationoptimized_A, optimized_d_1,optimized_d_2, optimized_f, optimized_tau_d_1,optimized_tau_d_2, optimized_tau_f, observed_intervals

#simulate_LIF_neuron(A, d_1, d_2, f, tau_d_1, tau_d_2, tau_f, dt, T, S_input, taurise, taudecay)

time, Vs, spike_times, Hs, Ds, Fs, xrise, xdecay = simulate_LIF_neuron(A, d_1,d_2, f, d_1,d_2, tau_f, dt, T, S_input, taurise, taudecay)

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



p4 = plot(time[450:600], synInputs[450:600], label="Synaptic Input", color=:red, xlabel="Time (ms)", ylabel="Synaptic Input", left_margin=20mm)

p5 = plot(time[450:600], xdecay[450:600], label="xdecay", color=:red, xlabel="Time (ms)", ylabel="xdecay", left_margin=20mm)

p6 = plot(time[450:600], xrise[450:600], label="xrise", color=:red, xlabel="Time (ms)", ylabel="xrise", left_margin=20mm)

plot_input = plot(p4, p5, p6, layout=(3,1), legend=:topright)

savefig(plot_input, "input_zoomed.png")
