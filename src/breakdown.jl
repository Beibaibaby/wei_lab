using Statistics
using Plots  
using Dates  # for generating timestamps
using JSON3  # Use JSON3 package for JSON handling
using Random
using Profile
using ProfileView

include("sim_inject_2.0.0.jl")

#plot or not
doplot = true
do_v_corr = false
do_repeat_v_corr=false

#Setting the parameters
###################################################

struct NetworkParameters
    Ncells::Int
    Ne::Int
    Ni::Int
    T::Int
    taue::Int
    taui::Int
    pei::Float64
    pie::Float64
    pii::Float64
    pee::Float64
    K::Int
    jie_para::Float64
    jei_para::Float64
    jii_para::Float64
    jee_para::Float64
    Nstim::Int
    stimstr_para::Float64
end

# Now, using the provided values, create an instance of the struct:
params = NetworkParameters(
    5000,  # Ncells
    4000,  # Ne
    1000,  # Ni
    100,  # T
    15,    # taue
    10,    # taui
    0.5,   # pei
    0.5,   # pie
    0.5,   # pii
    0.2,   # pee
    800,   # K
    4.0,   # jie
    -18.0 * 1.2,  # jei
    -16.0,  # jii
    10.0,   # jee
    400,    # Nstim
    0       # stimstr
)


#store it
#run the stimulus
#times, ns, Ne, Ncells, T, v_history, E_input, I_input, weights = sim_old()
times, ns, Ne, Ncells, T, v_history, E_input, I_input, weights,weights_D_history, weights_F_history=sim_dynamic_EI(params.Ne,params.Ni,params.T,params.taue,params.taui,params.pei,params.pie,params.pii,params.pee,params.K,params.stimstr_para,params.Nstim,params.jie_para,params.jei_para,params.jii_para,params.jee_para)

