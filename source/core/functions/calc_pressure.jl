function calc_pressure(stress::Vector{Float64}, E_kin::Vector{Float64}, box::Vector{Float64})
    virial = stress[1] + stress[2] + stress[3]
    V = box[1] * box[2] * box[3]
    curr_p = 2 / V * (E_kin[1] + E_kin[2] + E_kin[3] + virial) / 3
    return curr_p
end