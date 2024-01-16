
function sim(K,Ne,Ni)
	println("setting up parameters")
	 #K is the average number connections per neuron
     #Ne is the Number of E cells
     #Ni is the Number of I cells
	Ncells=Ne+Ni
	T = 3000 #simulation time (ms)

	taue = 15 #membrane time constant for exc. neurons (ms)
	taui = 10 #membrane time constant for inh. neurons (ms)

	#connection probabilities
	
	p = K/Ncells
	pei = p
	pie = p
	pii = p
    pee = p

	sqrtK = sqrt(K) 

	#Nepop = 80 #Number of neurons in a cluster(block) of E neuron

	jie = 4. / (taui*sqrtK)        #strength from I to E
	jei = -16. * 1.2 /(taue*sqrtK) #strength from E to I
	jii = -16. / (taui*sqrtK)      #strength from I to I 

	#set up connection probabilities within and without blocks
	#ratioejee = 1.9 
	jeeout = 10. / (taue*sqrtK)
	#jeein = ratioejee*10. / (taue*sqrtK)

	ratiopee = 2.5
	#pee = K/(Nepop*(ratiopee-1) + Ne)
	#peein = ratiopee*pee

	#stimulation
	Nstim = 400 #number of neurons to stimulate (indices 1 to Nstim will be stimulated)
	stimstr = .07/taue
	stimstart = T-500
	stimend = T 

	#constant bias to each neuron type
	muemin = 1.1
	muemax = 1.2
	muimin = 1
	muimax = 1.05

	vre = 0. #reset voltage

	threshe = 1 #threshold for exc. neurons
	threshi = 1

	dt = .1 #simulation timestep (ms)
	refrac = 5 #refractory period (ms)

	#synaptic time constants (ms)
	tauerise = 1
	tauedecay = 3
	tauirise = 1
	tauidecay = 2

	maxrate = 100 #(Hz) maximum average firing rate.  if the average firing rate across the simulation for any neuron exceeds this value, some of that neuron's spikes will not be saved

	mu = zeros(Ncells)
	thresh = zeros(Ncells)
	tau = zeros(Ncells)

	mu[1:Ne] = (muemax-muemin)*rand(Ne) .+ muemin
	mu[(Ne+1):Ncells] = (muimax-muimin)*rand(Ni) .+ muimin

	thresh[1:Ne] .= threshe
	thresh[(1+Ne):Ncells] .= threshi

	tau[1:Ne] .= taue
	tau[(1+Ne):Ncells] .= taui

	#Npop = round(Int,Ne/Nepop)

	weights = zeros(Ncells,Ncells)

	#random connections
	weights[1:Ne,1:Ne] = jeeout*(rand(Ne,Ne) .< pee)
	weights[1:Ne,(1+Ne):Ncells] = jei*(rand(Ne,Ni) .< pei)
	weights[(1+Ne):Ncells,1:Ne] = jie*(rand(Ni,Ne) .< pie)
	weights[(1+Ne):Ncells,(1+Ne):Ncells] = jii*(rand(Ni,Ni) .< pii)

	#connections within cluster
	#for pii = 1:Npop
	#	ipopstart = 1 + Nepop*(pii-1)
	#	ipopend = pii*Nepop
    #
	#	weights[ipopstart:ipopend,ipopstart:ipopend] = jeein*(rand(Nepop,Nepop) .< peein)
	#end

	for ci = 1:Ncells
		weights[ci,ci] = 0
	end

	#maxTimes = round(Int,maxrate*T/1000)
	times = zeros(Ncells,T)
	ns = zeros(Int,Ncells)

	forwardInputsE = zeros(Ncells) #summed weight of incoming E spikes
	forwardInputsI = zeros(Ncells)
	forwardInputsEPrev = zeros(Ncells) #as above, for previous timestep
	forwardInputsIPrev = zeros(Ncells)

	xerise = zeros(Ncells) #auxiliary variables for E/I currents (difference of exponentials)
	xedecay = zeros(Ncells)
	xirise = zeros(Ncells)
	xidecay = zeros(Ncells)

	v = rand(Ncells) #membrane voltage 

	lastSpike = -100*ones(Ncells) #time of last spike

	Nsteps = round(Int,T/dt)
		
	println("starting simulation")

	#begin main simulation loop
	for ti = 1:Nsteps
		if mod(ti,Nsteps/100) == 1  #print percent complete
			print("\r",round(Int,100*ti/Nsteps))
		end
		t = dt*ti
		forwardInputsE[:] .= 0
		forwardInputsI[:] .= 0
		
		for ci = 1:Ncells
			xerise[ci] += -dt*xerise[ci]/tauerise + forwardInputsEPrev[ci]
			xedecay[ci] += -dt*xedecay[ci]/tauedecay + forwardInputsEPrev[ci]
			xirise[ci] += -dt*xirise[ci]/tauirise + forwardInputsIPrev[ci]
			xidecay[ci] += -dt*xidecay[ci]/tauidecay + forwardInputsIPrev[ci]

			synInput = (xedecay[ci] - xerise[ci])/(tauedecay - tauerise) + (xidecay[ci] - xirise[ci])/(tauidecay - tauirise)

			if (ci < Nstim) && (t > stimstart) && (t < stimend) 
				synInput += stimstr;
			end

			if t > (lastSpike[ci] + refrac)  #not in refractory period
				v[ci] += dt*((1/tau[ci])*(mu[ci]-v[ci]) + synInput)

				if v[ci] > thresh[ci]  #spike occurred
					v[ci] = vre
					lastSpike[ci] = t
					ns[ci] = ns[ci]+1
					if ns[ci] <= T
						times[ci,ns[ci]] = t
					end

					for j = 1:Ncells
						if weights[j,ci] > 0  #E synapse
							forwardInputsE[j] += weights[j,ci]
						elseif weights[j,ci] < 0  #I synapse
							forwardInputsI[j] += weights[j,ci]
						end
					end #end loop over synaptic projections
				end #end if(spike occurred)
			end #end if(not refractory)
		end #end loop over neurons

		forwardInputsEPrev = copy(forwardInputsE)
		forwardInputsIPrev = copy(forwardInputsI)

	end #end loop over time
	print("\r")
	#println("Size of times: ", size(times))
	#println("Maximum ns: ", maximum(ns))
	#println(ns)
	# if maxTimes>maximum(ns)+1
	#  times = times[:,1:maximum(ns)]
	# end

	return times,ns,Ne,Ncells,T
end
