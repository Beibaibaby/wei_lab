using Plots
using Measures
using Optim

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

    Ds = Float64[]
    Fs = Float64[] 

    prev_V = V 
    spike_count = 0  

    for (idx, t) in enumerate(time)
        # Store the D and F at each timestep
        push!(Ds, D)
        push!(Fs, F)

        W = A * D * F
        dV = (-(V - V_rest) + R_m * W * S_input[idx]) / τ_m
        V += dV * dt

        if S_input[idx] == 1
            spike_count += 1
            if spike_count == 1 || spike_count == 2
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
    
    return time, Vs, spike_times, Hs, Ds, Fs  # Added Ds and Fs to the return values
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

# Example usage
A = 20.0        # fixed parameter A
d = 0.24       # depression fraction upon a spike
f = 0.00        # facilitation increment upon a spike
tau_d = 103.0     # time constant for D to recover to 1 (ms)
tau_f = 96.0     # time constant for F to recover to 1 (ms)
dt = 1.0        # time step (ms)
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




# Updated intervals
intervals = [50, 100, 200, 400, 800, 1600] 
pprs = Float64[]  # Store paired pulse ratios

# ... [Rest of your code] ...

for interval in intervals
    local S_input = generate_two_spike_train(T, dt, first_spike_time, interval)
    time, Vs, spike_times, Hs, Ds, Fs = simulate_LIF_neuron(A, d, f, tau_d, tau_f, dt, T, S_input)
    ppr = Hs[2] / Hs[1]
    push!(pprs, ppr)
end





# Create the plot
ppr_plot = plot(intervals, pprs, title="Paired Pulse Ratio vs. Interval", xlabel="Interval (ms)", ylabel="Paired Pulse Ratio", legend=false, linewidth=2, marker=:circle)

# Save the plot to a file
savefig(ppr_plot, "Paired_Pulse_Ratio_Plot.png")  # Saves the plot as a PNG file

lower_bounds = [0.0, 0.0] # Lower bounds for d and f, respectively
upper_bounds = [1.0, 5.0] # Upper bounds for d and f, you can adjust these as necessary

# Define the box constraints
box = Fminbox(BFGS())


# Observed data
observed_intervals = [800, 100, 400, 200, 1600, 50]
observed_pprs = [0.78506405666351, 0.706324711130814, 0.768659152138408, 0.70210158333243, 1.06176397890907,0.602094787420307]



# Cost function
function cost_function(params)
    d, f = params
    simulated_pprs = Float64[]

    for interval in observed_intervals
        local S_input = generate_two_spike_train(T, dt, first_spike_time, interval)
        time, Vs, spike_times, Hs, Ds, Fs = simulate_LIF_neuron(A, d, f, tau_d, tau_f, dt, T, S_input)
        ppr = Hs[2] / Hs[1]
        push!(simulated_pprs, ppr)
    end

    # Calculate the difference (e.g., mean squared error) between observed and simulated PPRs
    return sum((simulated_pprs - observed_pprs).^2) / length(observed_pprs)
end

# Optimization
initial_guess = [0.5, 0.1]  # Initial guess for d and f
result = optimize(cost_function, lower_bounds, upper_bounds, initial_guess, box, Optim.Options(g_tol = 1e-6))
optimized_d, optimized_f = Optim.minimizer(result)

println("Optimized d: $optimized_d, Optimized f: $optimized_f")




# Simulate PPRs using optimized parameters
function simulate_pprs(d, f, intervals)
    simulated_pprs = Float64[]
    for interval in intervals
        S_input = generate_two_spike_train(T, dt, first_spike_time, interval)
        time, Vs, spike_times, Hs, Ds, Fs = simulate_LIF_neuron(A, d, f, tau_d, tau_f, dt, T, S_input)
        ppr = Hs[2] / Hs[1]
        push!(simulated_pprs, ppr)
    end
    return simulated_pprs
end

simulated_pprs = simulate_pprs(optimized_d, optimized_f, observed_intervals)



# Plotting
p = plot(legend=:topright, xlabel="Interval (ms)", ylabel="Paired Pulse Ratio")

# Plot observed data
scatter!(p, observed_intervals, observed_pprs, label="Observed Data", color=:blue, marker=:circle)

# Plot simulated data
scatter!(p, observed_intervals, simulated_pprs, label="Simulated Data", color=:red, marker=:square)

# Save and display the plot
savefig(p, "PPR_Comparison_Plot.png")
