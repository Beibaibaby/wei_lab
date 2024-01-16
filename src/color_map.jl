using Plots

# Include your existing functions here: simulate_LIF_neuron_fixed_F, generate_spike_train
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

function find_theta_for_df_pair(A, d, f, tau_d, tau_f, dt, T, initial_spike_time, df_threshold=0.6)
    temporal_frequency = 1.0
    step = 1  # Increment step for temporal frequency
    max_tf = 175.0  # Maximum temporal frequency to try
    depression_ratio = 1.0

    while temporal_frequency <= max_tf
        S_input = generate_spike_train(T, dt, initial_spike_time, temporal_frequency)
        time, Vs, spike_times, Hs = simulate_LIF_neuron(A, d, f, tau_d, tau_f, dt, T, S_input)
        
        if length(Hs) >= 2
            depression_ratio = Hs[2] / Hs[1]
            if depression_ratio <= df_threshold
                return temporal_frequency
            end
        end
        
        temporal_frequency += step
    end
    
    return NaN  # Return NaN if no suitable temporal frequency is found
end

function generate_theta_colormap(d_values, f_values, A, tau_d, tau_f, dt, T, initial_spike_time, filename)
    theta_matrix = zeros(length(f_values), length(d_values))

    for (i, d) in enumerate(d_values)
        for (j, f) in enumerate(f_values)
            theta = find_theta_for_df_pair(A, d, f, tau_d, tau_f, dt, T, initial_spike_time)
            theta_matrix[j, i] = theta
        end
    end

    heatmap_plot = heatmap(d_values, f_values, theta_matrix, xlabel="d", ylabel="f", title=" ", color=:viridis,grid=false,dpi=300,size=(300,300))
    savefig(heatmap_plot, filename)  # Save the plot to a file
end

# Usage example

# Define parameters
A = 20.0
tau_d = 103.0
tau_f = 96.0
dt = 1.0
T = 3000.0
initial_spike_time = 50.0

# Define d and f values
# Generate d and f values
d_values = 0.1:0.005:0.5  # d from 0.1 to 1.0, incrementing by 0.05
f_values = 0.1:0.01:1.0   # f from 0.1 to 2.0, incrementing by 0.1
# Generate colormap
generate_theta_colormap(d_values, f_values, A, tau_d, tau_f, dt, T, initial_spike_time,"color_map.png")
