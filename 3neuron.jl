using Plots
using Measures


# Modified neuron simulation function for network context
function simulate_LIF_neuron_network(A, d_1, d_2, f, tau_d_1, tau_d_2, tau_f, dt, T, taurise, taudecay, external_input, inhibitory_input)
    τ_m = 10.0        # Membrane time constant
    V_thresh = -50.0  # Voltage threshold for spike
    V_reset = -75.0   # Reset voltage post-spike
    V_rest = -75.0    # Resting membrane potential
    R_m = 10.0        # Membrane resistance
    Inhibitory_weight = -0.5 # Weight for inhibitory effect

    # Time vector
    time = 0:dt:T
    n = length(time)

    # Initialize variables
    V = V_rest                # Membrane potential
    Vs = -75.0*ones(n)             # Array to store membrane potential over time
    spikes = zeros(Int, n)    # Spike times as binary array
    D_1, D_2, F = 1.0, 1.0, 1.0 # Synaptic variables
    xrise, xdecay = zeros(n), zeros(n)

    for idx in 2:n
        # Update synaptic input variables based on external and inhibitory inputs
        xrise[idx] += external_input[idx] - inhibitory_input[idx] * Inhibitory_weight
        xdecay[idx] += external_input[idx] - inhibitory_input[idx] * Inhibitory_weight

        xrise[idx] += -dt * xrise[idx-1] / taurise + xrise[idx-1]
        xdecay[idx] += -dt * xdecay[idx-1] / taudecay + xdecay[idx-1]

        synInput = (xdecay[idx] - xrise[idx]) / (taudecay - taurise)

        # Compute membrane potential
        W = A * D_1 * D_2 * F
        dV = (-(V - V_rest) + R_m * W * synInput) / τ_m
        V += dV * dt

        # Spike occurred
        if V >= V_thresh
            V = V_reset
            spikes[idx] = 1
        end

        # Update synaptic variables
        if external_input[idx] == 1
            D_1 *= d_1
            D_2 *= d_2
            F += f
        end

        # Update depression and facilitation factors
        D_1 += (1 - D_1) / tau_d_1 * dt
        D_2 += (1 - D_2) / tau_d_2 * dt
        F += (1 - F) / tau_f * dt

        Vs[idx] = V
    end
    
    return time, Vs, spikes
end

# Initialize parameters
A, d_1, d_2, f = 50.0, 0.9, 0.9, 0.1
tau_d_1, tau_d_2, tau_f = 100.0, 100.0, 100.0
dt, T = 0.1, 1000.0
taurise, taudecay = 1.0, 5.0

# External spiking inputs for each neuron
external_input_1 = [rand() < 0.05 ? 1 : 0 for _ in 0:dt:T]
external_input_2 = [rand() < 0.05 ? 1 : 0 for _ in 0:dt:T]
external_input_3 = [rand() < 0.05 ? 1 : 0 for _ in 0:dt:T]

# Simulate the network
time, Vs_1, spikes_1 = simulate_LIF_neuron_network(A, d_1, d_2, f, tau_d_1, tau_d_2, tau_f, dt, T, taurise, taudecay, external_input_1, zeros(length(external_input_1)))
_, Vs_2, _ = simulate_LIF_neuron_network(A, d_1, d_2, f, tau_d_1, tau_d_2, tau_f, dt, T, taurise, taudecay, external_input_2, spikes_1)
_, Vs_3, _ = simulate_LIF_neuron_network(A, d_1, d_2, f, tau_d_1, tau_d_2, tau_f, dt, T, taurise, taudecay, external_input_3, spikes_1)

# Initialize a new plot with 3 subplots for each neuron's membrane potential
p1 = plot(time, Vs_1, label="Neuron 1 Membrane Potential", xlabel="Time (ms)", ylabel="MP (mV)",left_margin = 10mm, bottom_margin = 10mm)
p2 = plot(time, Vs_2, label="Neuron 2 Membrane Potential", xlabel="Time (ms)", ylabel="MP (mV)",left_margin = 10mm, bottom_margin = 10mm)
p3 = plot(time, Vs_3, label="Neuron 3 Membrane Potential", xlabel="Time (ms)", ylabel="MP (mV)",left_margin = 10mm, bottom_margin = 10mm)

# Combine the plots into a single figure with subplots arranged vertically
combined_plot = plot(p1, p2, p3, layout = (3, 1), legend = false)

# Save the figure
savefig(combined_plot, "3neuron.png")
