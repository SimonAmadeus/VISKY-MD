function calc_totalmomentum(particle_vec::Vector{Particle}, velocity_vec)
    P_tot_x::Float64 = 0
    P_tot_y::Float64 = 0
    P_tot_z::Float64 = 0
    P_tot::Float64 = 0
    for i = 1:length(velocity_vec)
        P_tot_x += velocity_vec[i][1] * particle_vec[i].mass
        P_tot_y += velocity_vec[i][2] * particle_vec[i].mass
        P_tot_z += velocity_vec[i][3] * particle_vec[i].mass
    end
    P_tot = P_tot_x + P_tot_y + P_tot_z
    return P_tot
end

function calc_momentum_mag(particle_vec::Vector{Particle}, velocity_vec)
    P_tot_x::Float64 = 0
    P_tot_y::Float64 = 0
    P_tot_z::Float64 = 0
    for i = 1:length(velocity_vec)
        P_tot_x += velocity_vec[i][1] * particle_vec[i].mass
        P_tot_y += velocity_vec[i][2] * particle_vec[i].mass
        P_tot_z += velocity_vec[i][3] * particle_vec[i].mass
    end
    P_mag = sqrt(P_tot_x^2 + P_tot_y^2 + P_tot_z^2)
    return P_mag
end

function velocity_rescale!(particle_vec::Vector{Particle}, velocity_vec)
    # Calculate total momentum.
    P_tot_x::Float64 = 0
    P_tot_y::Float64 = 0
    P_tot_z::Float64 = 0
    total_mass::Float64 = 0
    for i = 1:length(velocity_vec)
        P_tot_x += velocity_vec[i][1] * particle_vec[i].mass
        P_tot_y += velocity_vec[i][2] * particle_vec[i].mass
        P_tot_z += velocity_vec[i][3] * particle_vec[i].mass
        total_mass += particle_vec[i].mass
    end

    # Calculate average velocity.
    avg_velocity_x = P_tot_x / total_mass
    avg_velocity_y = P_tot_y / total_mass
    avg_velocity_z = P_tot_z / total_mass

    # Subtract average velocity from each particle's velocity.
    for i = 1:length(velocity_vec)
        velocity_vec[i][1] -= avg_velocity_x
        velocity_vec[i][2] -= avg_velocity_y
        velocity_vec[i][3] -= avg_velocity_z
    end
end

function momentum_drift!(∆t, momentum_drift, P_mag, P_mag_prev)
    current_drift = (P_mag - P_mag_prev) / ∆t
#=    
    println("current: ", current_drift)
    println("P: ", P_tot)
    println("prev: ", P_tot_prev)
=#
    append!(momentum_drift, current_drift)
end