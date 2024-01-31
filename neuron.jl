using Plots
include("main.jl")

# Set the simulation parameters
A = 20.0        # Amplitude of synaptic input
d_1 = 0.24      # Depression factor 1
d_2 = 0.24      # Depression factor 2
f = 0.00        # Facilitation increment
tau_d_1 = 103.0 # Time constant for recovery of depression 1
tau_d_2 = 103.0 # Time constant for recovery of depression 2
tau_f = 96.0    # Time constant for recovery of facilitation
dt = 0.01        # Time step (ms)
T = 1000.0      # Total simulation time (ms)
taurise = 1.0   # Time constant for synaptic rise
taudecay = 3.0  # Time constant for synaptic decay

# Generate a spike train with 5 spikes
first_spike_time = 50.0  # ms
temporal_frequency = 10.0  # Hz (adjust as needed for 5 spikes within T)
S_input = generate_spike_train(T, dt, first_spike_time, temporal_frequency)

# Run the simulation
time, Vs, spike_times, Hs, Ds, Fs, xrise, xdecay = simulate_LIF_neuron(A, d_1, d_2, f, tau_d_1, tau_d_2, tau_f, dt, T, S_input, taurise, taudecay)

# Calculate synaptic input over time
synInputs = [(xdecay - xrise) / (taudecay - taurise) for (xrise, xdecay) in zip(xrise, xdecay)]

# Plot membrane potential and synaptic input

p1 = plot(time, Vs, label="Membrane Potential (V)", color=:blue, xlabel="Time (ms)", ylabel="Membrane Potential (mV)", left_margin=20mm)

# Plot synaptic input
p2 = plot(time, synInputs, label="Synaptic Input", color=:red, xlabel="Time (ms)", ylabel="Synaptic Input", left_margin=20mm)

plot(p1, p2, layout=(2,1), legend=:topright)
savefig("Neuron_Simulation.png")