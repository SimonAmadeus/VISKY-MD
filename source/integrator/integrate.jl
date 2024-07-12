# Velocity Verlet intgratior.
function update_positions!(args::System, particle_vec::Vector{Particle}, velocity_vec, force_vec)
    # Update positions.
    #@inbounds @simd for i = 1:length(particle_vec)
    for i = 1:length(particle_vec)
        particle_vec[i].pos .+= velocity_vec[i] .* args.∆t
        # Check for infinity
        if any(isinf, particle_vec[i].pos)
            println("Particle $i position is Inf.")
            break
        end
    end 
    # Update images.
    #@inbounds @simd for i = 1:length(particle_vec)
    for i = 1:length(particle_vec)
        if any(abs.(particle_vec[i].pos) .> 1e3) 
            @warn "Numerical instability detected in particle's position."
            println(particle_vec[i])
            continue  # Skip updating this particle, or do some corrective action.
        end
        particle_vec[i].image .+= floor.((particle_vec[i].pos .+ args.box ./ 2) ./ args.box)
        particle_vec[i].pos   .-= args.box .* floor.((particle_vec[i].pos .+ args.box ./ 2) ./ args.box)
        # Restrict positions.
        particle_vec[i].pos   .= clamp.(particle_vec[i].pos, - args.box / 2, args.box / 2)
    end
end

function update_velocities!(args::System, particle_vec::Vector{Particle}, velocity_vec, force_vec)
    # Update velocities.
    @inbounds @simd for i = 1:length(particle_vec)
        velocity_vec[i] .+= 0.5 .* (force_vec[i] ./ particle_vec[i].mass) .* args.∆t
    end 
end
