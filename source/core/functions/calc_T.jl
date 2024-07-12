function calc_inst_temp(particle_vec::Vector{Particle}, K_tot::Float64)
    dim::Int64 = 3
    T::Float64 = 2 * K_tot / (length(particle_vec) * dim)
    return T
end
