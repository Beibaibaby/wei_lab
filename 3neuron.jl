using Plots
using Measures

# Revised neuron simulation function for network context
function simulate_LIF_neuron_network_modified(A, d_1, d_2, f, tau_d_1, tau_d_2, tau_f, dt, T, tauerise, tauedecay, tauirise, tauidecay, external_input, inhibitory_input)
    τ_m = 10.0        # Membrane time constant
    V_thresh = -50.0  # Voltage threshold for spike
    V_reset = -75.0   # Reset voltage post-spike
    V_rest = -75.0    # Resting membrane potential
    R_m = 10.0        # Membrane resistance

    # Time vector
    time = 0:dt:T
    n = length(time)

    # Initialize variables
    V = V_rest                # Membrane potential
    Vs = -75.0*ones(n)        # Array to store membrane potential over time
    spikes = zeros(Int, n)    # Spike times as binary array
    D_1, D_2, F = 1.0, 1.0, 1.0 # Synaptic variables for inhibitory input
    xerise, xedecay = zeros(n), zeros(n)  # Excitatory synaptic variables
    xirise, xidecay = zeros(n), zeros(n)  # Inhibitory synaptic variables

    for idx in 2:n
        # Update excitatory synaptic input variables
        xerise[idx] = xerise[idx-1] -dt * xerise[idx-1] / tauerise + external_input[idx]
        xedecay[idx] = xedecay[idx-1]-dt * xedecay[idx-1] / tauedecay + external_input[idx]

        # Update inhibitory synaptic input variables
        xirise[idx] = xirise[idx-1] -dt * xirise[idx-1] / tauirise + inhibitory_input[idx]
        xidecay[idx] = xidecay[idx-1] -dt * xidecay[idx-1] / tauidecay + inhibitory_input[idx]

        # Calculate separate synaptic inputs for excitatory and inhibitory
        synInput_E = (xedecay[idx] - xerise[idx]) / (tauedecay - tauerise)
        synInput_I = (xidecay[idx] - xirise[idx]) / (tauidecay - tauirise)
        
        # Combine excitatory and inhibitory inputs
        synInput = A  * synInput_E - 3*A * D_1 * D_2 * F * synInput_I

        # Compute membrane potential
        dV = (-(V - V_rest) + R_m * synInput) / τ_m
        V += dV * dt

        # Spike occurred
        if V >= V_thresh
            V = V_reset
            spikes[idx] = 1
            # Update synaptic variables if the spike is due to inhibitory input
            D_1 *= d_1
            D_2 *= d_2
            F += f
        end

        # Update depression and facilitation factors for inhibitory inputs
        D_1 += (1 - D_1) / tau_d_1 * dt
        D_2 += (1 - D_2) / tau_d_2 * dt
        F += (1 - F) / tau_f * dt

        Vs[idx] = V
    end
    
    return time, Vs, spikes
end

# Initialize parameters
A, d_1, d_2, f = 13.0, 0.25, 0.99, 0.1
tau_d_1, tau_d_2, tau_f = 300.0, 300.0, 300.0


tauerise, tauedecay = 1.0, 5.0
tauirise, tauidecay = 2.0, 10.0

# External (excitatory) spiking inputs for each neuron
dt, T = 0.1, 4000.0
p1 = 0.0150
p2 = 0.0300

t_s = 500.0
# Defining the periods for each external input
period_1_end = 2000.0
period_2_end = period_1_end + t_s
period_3_end = period_2_end + t_s

# External input generation with two periods
external_input_1 = [rand() < (t <= period_1_end ? p1 : p2) for t in 0:dt:T]
external_input_2 = [rand() < (t <= period_2_end ? p1 : p2) for t in 0:dt:T]
external_input_3 = [rand() < (t <= period_3_end ? p1 : p2) for t in 0:dt:T]


# Simulate the network
time, Vs_1, spikes_1 = simulate_LIF_neuron_network_modified(A, d_1, d_2, f, tau_d_1, tau_d_2, tau_f, dt, T, tauerise, tauedecay, tauirise, tauidecay, external_input_1, zeros(length(external_input_1)))
_, Vs_2, spikes_2 = simulate_LIF_neuron_network_modified(A, d_1, d_2, f, tau_d_1, tau_d_2, tau_f, dt, T, tauerise, tauedecay, tauirise, tauidecay, external_input_2, spikes_1)
_, Vs_3, spikes_3 = simulate_LIF_neuron_network_modified(A, 0.001, 0.01, f, 1000, 1000, tau_f, dt, T, tauerise, tauedecay, tauirise, tauidecay, external_input_3, spikes_2)

# Function to find spike times for plotting
function spike_times(time, spikes)
    return time[findall(x -> x == 1, spikes)]
end
V_thresh = -50.0

# Plotting
p1 = plot(time, Vs_1, label="Neuron 1", xlabel="Time (ms)", ylabel=" (mV)", left_margin=15mm, bottom_margin=10mm)
scatter!(spike_times(time, spikes_1), V_thresh*ones(size(spike_times(time, spikes_1))), color=:red, markersize=1, label="Spikes 1")

p2 = plot(time, Vs_2, label="Neuron 2", xlabel="Time (ms)", ylabel=" (mV)", left_margin=15mm, bottom_margin=10mm)
scatter!(spike_times(time, spikes_2), V_thresh*ones(size(spike_times(time, spikes_2))), color=:red, markersize=1, label="Spikes 2")

p3 = plot(time, Vs_3, label="Neuron 3", xlabel="Time (ms)", ylabel=" (mV)", left_margin=15mm, bottom_margin=10mm)
scatter!(spike_times(time, spikes_3), V_thresh*ones(size(spike_times(time, spikes_3))), color=:red, markersize=1, label="Spikes 3") # Assuming spikes_3 exists

vspan!(p1, [2000, 2500], color = :green, alpha = 0.1, label = "")

vspan!(p1, [2500, 3000], color = :orange, alpha = 0.1, label = "")

# Apply similar blocks to p2 and p3 as needed
vspan!(p2, [2000, 2500], color = :green, alpha = 0.1, label = "")
vspan!(p2, [2500, 3000], color = :orange, alpha = 0.1, label = "")

vspan!(p3, [2000, 2500], color = :green, alpha = 0.1, label = "")
vspan!(p3, [2500, 3000], color = :orange, alpha = 0.1, label = "")


# Combine the plots into a single figure with subplots arranged vertically
combined_plot = plot(p1, p2, p3, layout = (3, 1), legend = false,dpi=300)

# Save the figure
savefig(combined_plot, "3neuron.png")

#####################Cut here for the next part#####################

time, Vs_4, spikes_4 = simulate_LIF_neuron_network_modified(A, d_1, d_2, f, tau_d_1, tau_d_2, tau_f, dt, T, tauerise, tauedecay, tauirise, tauidecay, external_input_1, zeros(length(external_input_1)))
_, Vs_5, spikes_5 = simulate_LIF_neuron_network_modified(A, d_1, d_2, f, tau_d_1, tau_d_2, tau_f, dt, T, tauerise, tauedecay, tauirise, tauidecay, external_input_2, zeros(length(external_input_1)))
_, Vs_6, spikes_6 = simulate_LIF_neuron_network_modified(A, 0.001, 0.01, f, 1000, 1000, tau_f, dt, T, tauerise, tauedecay, tauirise, tauidecay, external_input_3, spikes_5)


# Plotting
p4 = plot(time, Vs_4, label="Neuron 1", xlabel="Time (ms)", ylabel=" (mV)", left_margin=15mm, bottom_margin=10mm)
scatter!(spike_times(time, spikes_4), V_thresh*ones(size(spike_times(time, spikes_4))), color=:red, markersize=1, label="Spikes 1")

p5 = plot(time, Vs_5, label="Neuron 2", xlabel="Time (ms)", ylabel=" (mV)", left_margin=15mm, bottom_margin=10mm)
scatter!(spike_times(time, spikes_5), V_thresh*ones(size(spike_times(time, spikes_5))), color=:red, markersize=1, label="Spikes 2")

p6 = plot(time, Vs_6, label="Neuron 3", xlabel="Time (ms)", ylabel=" (mV)", left_margin=15mm, bottom_margin=10mm)
scatter!(spike_times(time, spikes_6), V_thresh*ones(size(spike_times(time, spikes_6))), color=:red, markersize=1, label="Spikes 3") # Assuming spikes_3 exists

vspan!(p4, [2000, 2500], color = :green, alpha = 0.1, label = "")
vspan!(p4, [2500, 3000], color = :orange, alpha = 0.1, label = "")

vspan!(p5, [2000, 2500], color = :green, alpha = 0.1, label = "")
vspan!(p5, [2500, 3000], color = :orange, alpha = 0.1, label = "")

vspan!(p6, [2000, 2500], color = :green, alpha = 0.1, label = "")
vspan!(p6, [2500, 3000], color = :orange, alpha = 0.1, label = "")



# Combine the plots into a single figure with subplots arranged vertically
combined_plot_2 = plot(p4, p5, p6, layout = (3, 1), legend = false,dpi=300)

# Save the figure
savefig(combined_plot_2, "3neuron_cko.png")
