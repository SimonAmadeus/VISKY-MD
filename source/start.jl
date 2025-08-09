# Initialize and start the simulation
# read control file
# call the run_MD function

include("./VISKY.jl")

using .VISKY

function parse_control_file(filename)
    params = Dict{Symbol, Any}()
    open(filename) do file
        for line in eachline(file)
            line = strip(split(line, '#')[1])  # Strip comments from line
            if !isempty(line)
                key, value = split(line, r"\s+", limit=2)
                params[Symbol(key)] = strip(value)
            end
        end
    end
    return params
end

function main()
    control_params = parse_control_file("control.txt")

    sys = sys_init()

    for (key, value) in control_params

        if key == :timestep
            sys.∆t = parse(Float64, strip(value))
        
        elseif key == :temperature
            sys.T = parse(Float64, strip(value))

        elseif key == :non_bonded_interactions

            interaction_type, args = split(value, r"\s+", limit = 2)
            
            if interaction_type == "original_hPF_Interactions"
                kappa, x, y, z = map(x -> parse(Float64, strip(x)), split(args, r"\s+"))
                sys.non_bonded_interactions = original_hPF_Interactions(kappa, [x, y, z])
            elseif interaction_type == "spectral_hPF_Interactions"
                kappa, sigma, x, y, z = map(x -> parse(Float64, strip(x)), split(args, r"\s+"))
                sys.non_bonded_interactions = spectral_hPF_Interactions(kappa, sigma, [x, y, z])
            
            elseif interaction_type == "LJ_12_6"
                ε, σ = map(x -> parse(Float64, strip(x)), split(args, r"\s+"))
                sys.non_bonded_interactions = LJ_12_6(ε, σ)
            
            elseif interaction_type == "DPD"
                a, γ = map(x -> parse(Float64, strip(x)), split(args, r"\s+"))
                sys.non_bonded_interactions = DPD(a, γ)
            end
        
        elseif key == :bonded_interactions
            
            interaction_type, args = split(value, r"\s+", limit = 2)

            if interaction_type == "Harmonic_Bonded"
                k, r_0 = map(x -> parse(Float64, strip(x)), split(args, r"\s+"))
                sys.bonded_interactions = Harmonic_Bonded(k, r_0)
            end

        elseif key == :angle_interactions
            
            interaction_type, args = split(value, r"\s+", limit = 2)

            if interaction_type == "Harmonic_Angles"
                k, θ_0 = map(x -> parse(Float64, strip(x)), split(args, r"\s+"))
                sys.angle_interactions = Harmonic_Angles(k, θ_0)
            end

        elseif key == :thermostat

            thermostat_type, arg = split(strip(value), r"\s+", limit = 2)
            if thermostat_type == "Berendsen_Thermostat"
                sys.thermostat = Berendsen_Thermostat(parse(Float64, arg))
            elseif thermostat_type == "Andersen_Thermostat"
                sys.thermostat = Andersen_Thermostat(parse(Float64, arg))
            elseif thermostat_type == "Lowe_Andersen_Thermostat"
                arg1, arg2 = split(arg, r"\s+")
                sys.thermostat = Lowe_Andersen_Thermostat(parse(Float64, arg1), parse(Float64, arg2))
            elseif thermostat_type == "Langevin_Thermostat"
                sys.thermostat = Langevin_Thermostat(parse(Float64, arg))    
            elseif thermostat_type == "Nose_Hoover_Thermostat"
                # Set Q initially to 1
                sys.thermostat = Nose_Hoover_Thermostat(1.0, parse(Float64, arg))    
            end

        elseif key == :collisions
            _, args = split(strip(value), r"\s+", limit = 2)
            x, y, z, dt, theta = map(x -> parse(Float64, strip(x)), split(args, r"\s+"))
            sys.collisions = MPCD([x, y, z], Int64(dt), theta)
        
        elseif key == :shear
            args = split(strip(value), r"\s+")
            mode, n, n_slabs, interval = args
            sys.shear = Shear_CTRL(mode, parse(Int64, n), parse(Int64, n_slabs), parse(Int64, interval))

        elseif key == :slip_springs
            n, n_md, n_mc, k, r_0, restart = split(strip(value), r"\s+")
            sys.t_bonds = T_Bonds_CTRL(parse(Int64, n), parse(Int64, n_md), parse(Int64, n_mc), parse(Float64, k), parse(Float64, r_0), restart)
               
        else
            # Assign the values to the corresponding fields of the sys structure
            field = Symbol(key)
            if field ∈ fieldnames(System)
                field_type = typeof(getfield(sys, field))
                if field_type == String
                    setfield!(sys, field, strip(value))
                else
                    setfield!(sys, field, parse(field_type, value))
                end
            end
        end
    end

    # Initialize logger
    logger = logger_init()
    logger.PotentialEnergy = true
    logger.KineticEnergy = true
    logger.Momentum = true
    sys.logger = logger
    sys.trajectorydump = XYZTrajectoryDump()
    #sys.trajectorydump = LammpsTrajectoryDump()

    sys.restart = true # Write restart file (restart.data)

    @time run_MD(sys)
    #@profile run_MD(sys)
    #@time Profile.@profile run_MD(sys)
    #Profile.print()
    #ProfileView.view()
    #ProfileVega.view()
end

main()