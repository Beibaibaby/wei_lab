#this file is part of litwin-kumar_doiron_cluster_2012
#Copyright (C) 2014 Ashok Litwin-Kumar
#see README for more information
using Statistics
using ProgressMeter
using Plots

include("./src/LIF_1.0.0_compare.jl")

function sim_dynamic(Ne,Ni,T,taue,taui,pei,pie,pii,pee,K,stimstr_para,Nstim,jie_para,
    jei_para,jii_para,jee_para, stim_duration,stim_start_time,ie_sign,ee_sign,corr_flag,
    add_noise,sigma_noise,scale_noise,d_ee,f_ee,d_ie,f_ie,
    stim_duration_2,stim_start_2,stimstr_2,c_noise)
    println("Setting up parameters")
    #corr_flag=false
    # Network parameters

    Ncells = Ne+Ni
    sqrtK = sqrt(K)
    if add_noise

       print("noise level is ")
       println(sigma_noise)

       print("noise lamba is ")
       println(scale_noise)
    end

    # Synaptic weights
    jie = jie_para / (taui * sqrtK)
    jei = jei_para / (taue * sqrtK)
    jii = jii_para / (taui * sqrtK)
    jee = jee_para / (taue * sqrtK)
   

    # Stimulation
    stimstr = stimstr_para/taue
    stimstr_2 = stimstr_2/taue
    print("stim is ")
    println(stimstr)
    print("stim 2 is ")
    println(stimstr_2)

    stimstart = stim_start_time
    stimend = stim_start_time+stim_duration

    stim_end_2 = stim_duration_2 + stim_start_2
    
    print("stim starting time is ")
    println(stim_start_time)

    print("stim ending time is ")
    println(stimend)

    print("stim 2 starting time is ")
    println(stim_start_2)

    print("stim 2 ending time is ")
    println(stim_end_2)


    #d = 0.15       # depression fraction upon a spike
    #f = 0.92  (orignal)      # facilitation increment upon a spike
    #f = 0.5
    #tau_d = 103.0     # time constant for D to recover to 1 (ms)
    #tau_f = 96.0     # time constant for F to recover to 1 (ms)
    
    tau_d = 1030 # time constant for D to recover to 1 (step;not ms)
    tau_f = 960 # time constant for D to recover to 1 (step;not ms)

    #todo: might need adjust the time scale

    # Constants and thresholds
    muemin = 1.1
    muemax = 1.2
    muimin = 1
    muimax = 1.05

    vre = 0.0
    threshe = 1
    threshi = 1

    dt = 0.1
    refrac = 5

    tauerise = 1
    taurise = 1
    tauedecay = 3
    tauidecay = 2

    maxrate = 500  # Maximum average firing rate (Hz) was 100 orignally

    # Initialize parameters
    mu = zeros(Ncells)
    thresh = zeros(Ncells)
    tau = zeros(Ncells)

    mu[1:Ne] .= (muemax - muemin) .* rand(Ne) .+ muemin
    mu[(Ne + 1):Ncells] .= (muimax - muimin) .* rand(Ni) .+ muimin

    thresh[1:Ne] .= threshe
    thresh[(Ne + 1):Ncells] .= threshi

    tau[1:Ne] .= taue
    tau[(Ne + 1):Ncells] .= taui

    weights = zeros(Ncells, Ncells)
    weights_D_ee = ones(Ne) #the sending D (E->E)
    weights_F_ee = ones(Ne) #the sending F (E->E)

    weights_D_ie = ones(Ne) #the sending F (E->I)
    weights_F_ie = ones(Ne) #the sending F (E->I)

    # Here we only need one decay/facilitation factor for one given neuron i, the factors from i to j are all the same

    # Random connections
    weights[1:Ne, 1:Ne] .= jee .* (rand(Ne, Ne) .< pee)
    weights[1:Ne, (1 + Ne):Ncells] .= jei .* (rand(Ne, Ni) .< pei)
    weights[(1 + Ne):Ncells, 1:Ne] .= jie .* (rand(Ni, Ne) .< pie)
    weights[(1 + Ne):Ncells, (1 + Ne):Ncells] .= jii .* (rand(Ni, Ni) .< pii)
    
    for ci = 1:Ncells
        weights[ci, ci] = 0
    end

    maxTimes = round(Int, maxrate * T / 1000)
    times = zeros(Ncells, maxTimes)
    ns = zeros(Int, Ncells)
    Nsteps = round(Int, T / dt)


    weights_D_ee_track = zeros(Nsteps)  
    weights_F_ee_track = zeros(Nsteps)

    weights_D_ie_track = zeros(Nsteps)  
    weights_F_ie_track = zeros(Nsteps)

    weights_IE_mean_history = zeros(Nsteps)
    weights_EE_mean_history = zeros(Nsteps)

    v_history = zeros(Ncells, Nsteps)  # Nsteps because we're saving at each time step, not just spikes
    E_input=zeros(Ncells, Nsteps) 
    I_input=zeros(Ncells, Nsteps)
    forwardInputsE = zeros(Ncells)
    forwardInputsI = zeros(Ncells)
    forwardInputsEPrev = zeros(Ncells)
    forwardInputsIPrev = zeros(Ncells)

    xerise = zeros(Ncells)
    xedecay = zeros(Ncells)
    xirise = zeros(Ncells)
    xidecay = zeros(Ncells)

    v = rand(Ncells)
    #v[1:3] .=0
    #v = vre*ones(Ncells)
    lastSpike = -100 * ones(Ncells)
    
    
    println("Starting simulation")
    pl = Progress(Nsteps, 5)
    
    corr_pairs=100
    weights_copy = copy(weights)

    for ti = 1:Nsteps
        t = dt * ti
        forwardInputsE[:] .= 0
        forwardInputsI[:] .= 0
        
        weights_D_ee .+= (1 .- weights_D_ee) ./ tau_d
        weights_F_ee .+= (1 .- weights_F_ee) ./ tau_f

        weights_D_ie .+= (1 .- weights_D_ie) ./ tau_d
        weights_F_ie .+= (1 .- weights_F_ie) ./ tau_f

        weights_D_ee_track[ti] = mean(weights_D_ee)
        weights_F_ee_track[ti] = mean(weights_F_ee)
        weights_D_ie_track[ti] = mean(weights_D_ie)
        weights_F_ie_track[ti] = mean(weights_F_ie)
        
        if ie_sign == true
           weights[(1 + Ne):Ncells, 1:Ne] = weights_copy[(1 + Ne):Ncells, 1:Ne] .* (weights_D_ie .* weights_F_ie)'
        end 

        if ee_sign == true
            weights[1:Ne, 1:Ne] = weights_copy[1:Ne, 1:Ne] .* (weights_D_ee .* weights_F_ee)'
        end 

        W_sub_view_ie = @view weights[(1 + Ne):Ncells, 1:Ne]
        weights_IE_mean_history[ti] = mean(W_sub_view_ie)
        W_sub_view_ee = @view weights[1:Ne, 1:Ne]
        weights_EE_mean_history[ti] = mean(W_sub_view_ee)
        
        if add_noise
            
           gaussian_noise_global = sqrt(c_noise) * scale_noise * randn() * sigma_noise * dt
           
        end

        for ci = 1:Ncells
            
            xerise[ci] += -dt * xerise[ci] / tauerise + forwardInputsEPrev[ci]
            xedecay[ci] += -dt * xedecay[ci] / tauedecay + forwardInputsEPrev[ci]
            xirise[ci] += -dt * xirise[ci] / taurise + forwardInputsIPrev[ci]
            xidecay[ci] += -dt * xidecay[ci] / tauidecay + forwardInputsIPrev[ci]
   
            #synInput = (xedecay[ci] - xerise[ci]) / (tauedecay - tauerise) + (xidecay[ci] - xirise[ci]) / (tauidecay - taurise)
            synInput_E = (xedecay[ci] - xerise[ci]) / (tauedecay - tauerise)
            synInput_I = (xidecay[ci] - xirise[ci]) / (tauidecay - taurise)
            E_input[ci, ti] = synInput_E
            I_input[ci, ti] = synInput_I 
            synInput = synInput_E+synInput_I # (xedecay[ci] - xerise[ci]) / (tauedecay - tauerise) + (xidecay[ci] - xirise[ci]) / (tauidecay - taurise)
            
            if corr_flag == true
              if ci<corr_pairs
                synInput=synInput_E+(-0.2+mu[ci])/tau[ci]
                #synInput = synInput_E
              end

              if ci >= 500 && ci <600
                synInput=synInput+(1-mu[ci])/tau[ci]
              end

              if ci >= 1001 && ci <1100
                synInput=synInput_I+(2-mu[ci])/tau[ci]
              end
            else
               synInput=synInput_E+synInput_I
            end

            if (ci < Nstim) && (t > stimstart) && (t < stimend)
                synInput += stimstr
            end
            

            if (ci < Nstim) && (t > stim_start_2) && (t < stim_end_2)
                synInput += stimstr_2
            end
            

            if add_noise

            
               gaussian_noise_local = sqrt(1-c_noise) * randn() * sigma_noise * dt
           

               synInput += gaussian_noise_global + gaussian_noise_local
            end
           
            if t > (lastSpike[ci] + refrac)
                v[ci] += dt * ((1 / tau[ci]) * (mu[ci] - v[ci]) + synInput)     
                v_history[ci, ti] = v[ci]
                if v[ci] > thresh[ci]
                    if corr_flag == true
                        if (ci>=corr_pairs && ci<500 )||(ci>=600 && ci<=1000)  || ci>=1100
                           v[ci] = vre
                        end
                    else
                        v[ci] = vre
                    end

                    lastSpike[ci] = t
                    ns[ci] += 1
                    if ns[ci] <= maxTimes
                        times[ci, ns[ci]] = t
                    end
                    
                    #test if remove
                   if ci < Ne
                    weights_D_ee[ci] = d_ee * weights_D_ee[ci]
                    weights_F_ee[ci] = f_ee + weights_F_ee[ci]
                    weights_D_ie[ci] = d_ie * weights_D_ie[ci]
                    weights_F_ie[ci] = f_ie + weights_F_ie[ci]
                   end                             

                    for j = 1:Ncells
                        if weights[j, ci] > 0  # E synapse
                            forwardInputsE[j] += weights[j, ci]
                        elseif weights[j, ci] < 0  # I synapse
                            forwardInputsI[j] += weights[j, ci]
                        end
                    end
                end
            else
                v_history[ci, ti] = v[ci]#
            end
        end

        forwardInputsEPrev .= forwardInputsE
        forwardInputsIPrev .= forwardInputsI
        next!(pl)
    end

    print("\r")
    println(maximum(ns))
    println(maximum(ns))
    if maximum(ns)<=maxTimes
         times = times[:, 1:maximum(ns)]
    else
       println("no over max")
    end
    
    return times, ns, Ne, Ncells, T, v_history, E_input, I_input, weights, weights_D_ee_track, weights_F_ee_track , weights_IE_mean_history, weights_EE_mean_history, weights_D_ie_track, weights_F_ie_track
end



function compute_correlation(E_input::Matrix, I_input::Matrix, num_pairs::Int=100)
    # Ensure input matrices have the same dimensions
    if size(E_input) != size(I_input)
        error("E_input and I_input must have the same dimensions.")
    end

    Ncells, Nsteps = size(E_input)

    # Ensure the number of pairs is less than or equal to Ncells
    if num_pairs > Ncells
        error("The number of pairs cannot exceed the number of cells.")
    end

    # Randomly select indices for both E_input and I_input
    random_e_indices = rand(1:Ncells, num_pairs)
    random_i_indices = rand(1:Ncells, num_pairs)

    total_correlation = 0.0

    for i = 1:num_pairs
        e_data = E_input[random_e_indices[i], :]
        i_data = I_input[random_i_indices[i], :]

        # Compute Pearson correlation and add to the total
        total_correlation += cor(e_data, i_data)
    end

    # Compute average correlation
    avg_correlation = total_correlation / num_pairs

    return avg_correlation
end

function compute_cross_correlation(E_input::Matrix, I_input::Matrix, tau_range::UnitRange{Int}=-15:15, num_pairs::Int=50)
    # Ensure input matrices have the same dimensions
    if size(E_input) != size(I_input)
        error("E_input and I_input must have the same dimensions.")
    end

    Ncells, Nsteps = size(E_input)
    
    num_pairs=round(Int,Ncells)


    # Ensure the number of pairs is less than or equal to Ncells
    if num_pairs > Ncells * Ncells
        error("The number of pairs cannot exceed the square of the number of cells.")
    end
    

    # Randomly select indices for both E_input and I_input
    random_e_indices = rand(1:Ncells, num_pairs)
    random_i_indices = rand(1:Ncells, num_pairs)

    avg_correlations = Dict{Int, Float64}()

    for tau in tau_range
        total_correlation = 0.0
        for i = 1:num_pairs
            e_data = E_input[random_e_indices[i], :]
            i_data = I_input[random_i_indices[i], :]
            
            # Shift data according to tau
            if tau > 0
                e_data = e_data[1:end-tau]
                i_data = i_data[tau+1:end]
            elseif tau < 0
                e_data = e_data[-tau+1:end]
                i_data = i_data[1:end+tau]
            end
            
            # Compute Pearson correlation for this tau and add to the total
            total_correlation += cor(e_data, i_data)
        end

        # Compute average correlation for this tau
        avg_correlations[tau] = total_correlation / num_pairs
    end

    return avg_correlations
end

    

function run_simulation_and_tune_params(K, Ne, Ni, T, jie, jee, jei, jii)
 params = NetworkParameters(K, Ne, Ni, T, sqrt(K), 10, 15, jie, jee, jei, jii)

 times, ns, Ne, Ncells, T = sim(params.K, params.Ne, params.Ni, params.T, params.jie, params.jei, params.jii, params.jee)

 excitatory_rate = mean(1000 * ns[1:params.Ne] / params.T)
 inhibitory_rate = mean(1000 * ns[(params.Ne + 1):Ncells] / params.T)

 return times, excitatory_rate, inhibitory_rate
end

function tune_parameters()
 doplot = true

 successful_adjustments = 0
 required_successful_adjustments = 50
 # Initial parameter values
 K = 800
 Ne = 4000
 Ni = 1000
 T = 1500
 sqrtK = sqrt(K)
 taue = 10
 taui = 15

 jie = 0.2 / (taui * sqrtK)
 jee = 0.4 / (taue * sqrtK)
 jei = -0.8 * 1.2 / (taue * sqrtK)
 jii = -0.8 / (taui * sqrtK)

 adjustment_step = 0.001  # Adjust step size based on sensitivity
 random_factor = 0.001   # Introduce randomness
 while successful_adjustments < required_successful_adjustments
 while true
    println("Round+1")
    times, excitatory_rate, inhibitory_rate = run_simulation_and_tune_params(K, Ne, Ni, T, jie, jee, jei, jii)

    if excitatory_rate <= 5 && inhibitory_rate <= 5 && inhibitory_rate < excitatory_rate && inhibitory_rate >=0.5
        println("Parameters tuned successfully:")
        println("jie = $jie, jee = $jee, jei = $jei, jii = $jii")
        println("Mean excitatory firing rate: $excitatory_rate Hz")
        println("Mean inhibitory firing rate: $inhibitory_rate Hz")
        successful_adjustments=successful_adjustments+1
        if doplot
            timestamp = Dates.now()
            timestamp_str = Dates.format(timestamp, "yyyy-mm-dd_HH-MM-SS")

            plot_size = (600, 400)

            window_size = 100
            step_size = 5

            e_rate = compute_sliding_rate(times[1:Ne, :], window_size, step_size, T)
            i_rate = compute_sliding_rate(times[(Ne + 1):Ncells, :], window_size, step_size, T)

            n_steps = length(e_rate)
            time_values = [i * step_size + (window_size / 2) for i in 1:n_steps]

            p2 = plot(time_values, e_rate, xlabel="Time (ms)", ylabel="Firing rate (Hz)",
                      label="Excitatory", lw=2, linecolor=:red, size=plot_size)
            plot!(time_values, i_rate, label="Inhibitory", lw=2, linecolor=:skyblue)

            fig_filename = "./figs/ss_$timestamp_str.png"
            savefig(p2, fig_filename)
            println("Figure saved as $fig_filename")

            json_filename = "./figs/ss_$timestamp_str.json"
            param_dict = Dict(
                "K" => K,
                "Ne" => Ne,
                "Ni" => Ni,
                "T" => T,
                "sqrtK" => sqrtK,
                "taue" => taue,
                "taui" => taui,
                "jie" => jie,
                "jee" => jee,
                "jei" => jei,
                "jii" => jii
            )
            JSON3.open(json_filename, "w") do io
                JSON3.print(io, param_dict)
            end
            println("Parameters saved as $json_filename")
        end
    

    else
        # Randomly adjust parameters with noise
        jie += adjustment_step * (-1 + 2*random_factor * randn())
        jee += adjustment_step * (-1 + 2*random_factor * randn())
        jei -= adjustment_step * (-1 + 2*random_factor * randn())
        jii -= adjustment_step * (-1 + 2*random_factor * randn())

        jie = max(0.0001, jie)
        jee = max(0.0001, jee)

        # Ensure jei and jii are less than 0
        jei = min(-0.0001, jei)
        jii = min(-0.0001, jii)
    end
 end
 end
end

#tune_parameters()



function compute_sliding_rate(spiketimes, window_size, step_size, T)
    println("starting computing firing rate")
    n_neurons, _ = size(spiketimes)
    n_steps = floor(Int, (T - window_size) / step_size) + 1
    rates = zeros(n_steps)

    p = Progress(n_steps, 2) # Initialize the progress bar

    for i = 1:n_steps
        next!(p) # Update the progress bar

        t_start = (i-1) * step_size + 1  # Ensure the start time is non-zero
        t_end = t_start + window_size - 1  # Adjust end time based on start
        
        # Check for out-of-bounds
        if t_end > T
            break
        end

        n_spikes = sum([sum((t_start .<= spiketimes[j, :]) .& (spiketimes[j, :] .< t_end)) for j=1:n_neurons])
        rates[i] = n_spikes / (n_neurons * window_size) * 1000  # rate in Hz
    end
    #println(rates)
    return rates
end




function plot_correlations(cross_corr_E_E, cross_corr_I_I, cross_corr_E_I, cross_corr_I_E, cross_corr_C_C, output_path)
    # Assuming the dictionaries have the same keys
    sorted_keys = sort(collect(keys(cross_corr_E_E)))
    
    sorted_values_E_E = [cross_corr_E_E[k] for k in sorted_keys]
    sorted_values_I_I = [cross_corr_I_I[k] for k in sorted_keys]
    sorted_values_E_I = [cross_corr_E_I[k] for k in sorted_keys]
    sorted_values_I_E = [cross_corr_I_E[k] for k in sorted_keys]
    sorted_values_C_C = [cross_corr_C_C[k] for k in sorted_keys]
    plot_margin = 25mm
    p_size = (1200, 800)

    plot(sorted_keys, sorted_values_E_E, label="E-E", linewidth=2, marker=:circle, color=:blue, legend=:topright,size=p_size,left_margin=plot_margin)
    plot!(sorted_keys, sorted_values_I_I, label="I-I", linewidth=2, marker=:circle, color=:red)
    plot!(sorted_keys, sorted_values_E_I, label="E-I", linewidth=2, marker=:circle, color=:yellow)
    plot!(sorted_keys, sorted_values_I_E, label="I-E", linewidth=2, marker=:circle, color=:orange)
    plot!(sorted_keys, sorted_values_C_C, label="C-C", linewidth=2, marker=:circle, color=:green)

    xlabel!("Tau")
    ylabel!("Correlation")
    title!("Cross-correlation")
    savefig(output_path)
end


function plot_correlations_mem(cross_corr_E_E, cross_corr_I_I, cross_corr_T_T, cross_corr_E_I, cross_corr_I_E, output_file)
    # Assuming the dictionaries have the same keys
    sorted_keys = sort(collect(keys(cross_corr_E_E)))
    
    sorted_values_E_E = [cross_corr_E_E[k] for k in sorted_keys]
    sorted_values_I_I = [cross_corr_I_I[k] for k in sorted_keys]
    sorted_values_T_T = [cross_corr_T_T[k] for k in sorted_keys]
    sorted_values_E_I = [cross_corr_E_I[k] for k in sorted_keys]
    sorted_values_I_E = [cross_corr_I_E[k] for k in sorted_keys]

    # Create the plot
    plot(sorted_keys, sorted_values_E_E, label="E-E", linewidth=2, marker=:circle, color=:blue, legend=:topright)
    plot!(sorted_keys, sorted_values_I_I, label="I-I", linewidth=2, marker=:circle, color=:red)
    plot!(sorted_keys, sorted_values_E_I, label="E-I", linewidth=2, marker=:circle, color=:yellow)
    plot!(sorted_keys, sorted_values_I_E, label="I-E", linewidth=2, marker=:circle, color=:orange)
    plot!(sorted_keys, sorted_values_T_T, label="T-T", linewidth=2, marker=:circle, color=:green)

    # Set labels and title
    xlabel!("Tau")
    ylabel!("Correlation")
    title!("Cross-correlation")

    # Save the figure to the specified output file
    savefig(output_file)
end


function plot_cells(v_history, cells)
    p = plot(layout=(3,1), size=(600,800))

    for (index, cell) in enumerate(cells)
        plot!(p[index], v_history[cell, :], label="Cell $cell", xlabel="Time", ylabel="Voltage", legend=:topright)
    end
    savefig("v_history.png")
end

function plot_curve(d_ee, f_ee, d_ie, f_ie, output_filename)
    A = 20.0
    tau_d = 103.0
    tau_f = 96.0
    dt = 1.0
    T = 3000.0
    first_spike_time = 50.0
    temporal_frequencies = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 20, 25, 30, 40, 50, 60, 75, 80, 100, 110, 130, 150]

    depression_ratios = []
    depression_ratios_fixed_F = []

    for tf in temporal_frequencies
        S_input = generate_spike_train(T, dt, first_spike_time, tf)

        # Simulate with variable F
        _, _, _, Hs = simulate_LIF_neuron(A, d_ee, f_ee, tau_d, tau_f, dt, T, S_input)
        if length(Hs) >= 2
            push!(depression_ratios, Hs[2] / Hs[1])
        end

        # Simulate with fixed F
        _, _, _, Hs_fixed_F = simulate_LIF_neuron_fixed_F(A, d_ie, f_ie, tau_d, tau_f, dt, T, S_input)
        if length(Hs_fixed_F) >= 2
            push!(depression_ratios_fixed_F, Hs_fixed_F[2] / Hs_fixed_F[1])
        end
    end

    # Plot and save results
    p = plot(temporal_frequencies, depression_ratios, xlabel="Temporal Frequency (Hz)", title="d_ee=$d_ee f_ee=$f_ee d_ie=$d_ie f_ie=$f_ie", ylabel="Depression Ratio", label="E->E", marker=:circle, markercolor=:red, linewidth=2, linecolor=:red)
    plot!(p, temporal_frequencies, depression_ratios_fixed_F, label="E->I", marker=:circle, linewidth=2, linecolor=:deepskyblue2, markercolor=:deepskyblue2)
    savefig(p, output_filename)
end
