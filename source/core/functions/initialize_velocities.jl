function initialize_velocities!(args::System, particle_vec::Vector{Particle}, velocity_vec)
    if args.init_velocity == "Z"
        for i = 1:length(particle_vec)
            velocity_vec[i] = zeros(3)
        end
        println("Initial velocities are set to zero.")
    elseif args.init_velocity == "G"
        for i = 1:length(particle_vec)
            d = Normal(0.0, sqrt(args.T / particle_vec[i].mass))
            velocity_vec[i] = rand(d,3)
        end
        # Remove any net momentum:
        total_momentum = [0.0, 0.0, 0.0]
        for i = 1:length(particle_vec)
            total_momentum += particle_vec[i].mass * velocity_vec[i]
        end
        avg_momentum = total_momentum / length(particle_vec)
        for i = 1:length(particle_vec)
            velocity_vec[i] -= avg_momentum / particle_vec[i].mass
        end
        println("Initial velocities are set by a Gaussian distribution.")
    elseif args.init_velocity == "R"
        println("Initial velocities are read from restart file.") 
    else
        println("Please set a mode for the initial velocity!")
    end
end