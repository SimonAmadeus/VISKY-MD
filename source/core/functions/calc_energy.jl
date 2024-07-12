function bond_energy(args::System, i::Int64, j::Int64, particle_vec::Vector{Particle})
    pos_1 = particle_vec[i].pos .+ particle_vec[i].image .* args.box
    pos_2 = particle_vec[j].pos .+ particle_vec[j].image .* args.box
    ∆d = pos_2 .- pos_1
    r = norm(∆d)
    energy = 0.5 * args.t_bonds.k * (r - args.t_bonds.r_0)^2
    return energy
end

function calc_kinetic_energy(particle_vec::Vector{Particle}, velocity_vec, E_kin::Vector{Float64})
    for i in 1:length(particle_vec)
        E_kin[1] += 0.5 * particle_vec[i].mass * velocity_vec[i][1]^2
        E_kin[2] += 0.5 * particle_vec[i].mass * velocity_vec[i][2]^2
        E_kin[3] += 0.5 * particle_vec[i].mass * velocity_vec[i][3]^2
    end
    return E_kin
end

function calc_total_kinetic_energy(particle_vec::Vector{Particle}, velocity_vec)
    K_tot::Float64 = 0
    for i = 1:length(velocity_vec)
        m = particle_vec[i].mass
        v = velocity_vec[i]
        K_tot += 0.5 * m * dot(v, v)
    end
    return K_tot
end