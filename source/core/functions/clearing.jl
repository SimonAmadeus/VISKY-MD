#function clear_force!(force_vec::Vector{Force})
#    for i = 1:length(force_vec)
#        force_vec[i].force[1] = 0.0
#        force_vec[i].force[2] = 0.0
#        force_vec[i].force[3] = 0.0
#    end
#end

function clear_force!(force_vec::Vector{Vector{Float64}})
    for f in force_vec
        fill!(f, 0.0)
    end
end

function clear_energy!(energy_vec::Vector{Float64})
    for i = 1:length(energy_vec)
        energy_vec[i] = 0.0
    end
end

function clear_stress_tensor!(stress)
    for i in 1:length(stress)
        stress[i] = 0.0
    end
end

function clear_kin_energy_tensor!(E_kin)
    for i in 1:length(E_kin)
        E_kin[i] = 0.0
    end
end
