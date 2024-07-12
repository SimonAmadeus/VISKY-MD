function apply_non_bonded_interactions!(args::System, ::Any, particle_vec::Vector{Particle}, c_l::Int64, force_vec, energy_vec::Vector{Float64}, stress::Vector{Float64}, mesh::Mesh, ::No_Non_Bonded) 
    return energy, stress
end
