using NPZ
using Plots
using Measures
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
            println(D_1)
            println(D_2)
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
            push!(Hs, W * synInput) # Change in potential due to the input spike
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


# Define the file paths
file_paths = [
    "./spikes/18515004.abf_ch0_sweep0.npy_spikes.npy",
    "./spikes/18515005.abf_ch0_sweep0.npy_spikes.npy",
    "./spikes/18515008.abf_ch0_sweep0.npy_spikes.npy",
    "./spikes/18515010.abf_ch0_sweep0.npy_spikes.npy"
]

# Simulation parameters (these are examples, adjust as needed)
A = 54.39
d_1 = 1.0
d_2 = 0.45
f = 1.0
tau_d_1 = 498.94
tau_d_2 = 473.2
tau_f = 242.49
dt = 0.1
T = 5000.0
taurise = 0.95
taudecay = 1.0

# Load spike data
function load_spike_data(file_path)
    spike_data = NPZ.npzread(file_path)
    return spike_data
end

# Prepare S_input for simulation
function prepare_S_input(spike_data, T, dt)
    time_steps = Int(T / dt) + 1
    S_input = zeros(Int, time_steps)
    for spike_time in findall(x -> x == 1, spike_data)
        if spike_time <= time_steps
            S_input[spike_time] = 1
        end
    end
    return S_input
end

# Function to run simulation and save plot for a given file path
# Corrected function to run simulation and save plot for a given file path
function simulate_and_plot(file_path)
    # Load and prepare spike input
    spike_data = load_spike_data(file_path)
    S_input = prepare_S_input(spike_data, T, dt)
    #println(S_input)
    # Run simulation
    time, Vs, spike_times, Hs, Ds, Fs, xrise, xdecay = simulate_LIF_neuron(A, d_1, d_2, f, tau_d_1, tau_d_2, tau_f, dt, T, S_input, taurise, taudecay)
    # Calculate synaptic inputs from xrise and xdecay values
    synInputs = [(xdecay - xrise) / (taudecay - taurise) for (xrise, xdecay) in zip(xrise, xdecay)]

    # Calculate Current using the synaptic inputs, depression (D), and facilitation (F) factors
    Current = 54.39 .* Ds .* Fs .* synInputs

# Plotting the Current over time
  p = plot(time, Current, label="Current", xlabel="Time (ms)", ylabel="Current (pA)", legend=:topright,size=(3000,500),left_margin=20mm,bottom_margin=15mm)

    # Save plot to file
    output_dir = "./spikes"
    if !isdir(output_dir)
        mkdir(output_dir)
    end
    base_name = splitext(basename(file_path))[1]  # Extract base name without extension
    save_path = joinpath(output_dir, base_name * "_simulation_plot.png")
    savefig(p, save_path)  # Make sure to use the plot object 'p' here
    println("Plot saved to $save_path")
end


# Iterate over file paths, simulating and plotting for each
for file_path in file_paths
    simulate_and_plot(file_path)
end
