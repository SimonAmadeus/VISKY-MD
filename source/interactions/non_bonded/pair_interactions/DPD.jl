export DPD

struct DPD <: Non_Bonded_Interactions
    a::Float64 # Repulstion parameter.
    γ::Float64 # Friction parametter.
end

function apply_non_bonded_interactions!(args::System, ::System_Type, particle_vec::Vector{Particle}, velocity_vec, neighbor_list::Any, c_l::Int64, force_vec, energy::Float64, stress::Vector{Float64}, mesh::Mesh, DPD::DPD)
    σ = sqrt(2 * args.T * DPD.γ)
    invsqrt_dt = 1 / sqrt(args.∆t)
    #dtfac = 1 / sqrt(args.∆t)
    #for i = 2:length(particle_vec)
        #for j = 1:i-1
    for i = 1:length(neighbor_list)
        for k = 1:length(neighbor_list[i])
            j = neighbor_list[i][k]
            if j > i
                ∆r = particle_vec[i].pos .- particle_vec[j].pos
                ∆r .-= args.box .* floor.((∆r .+ args.box / 2) ./ args.box) # PBCs
                r = norm(∆r)
                r_unit = ∆r / r
                if r < args.cutoff
                    rinv = 1 / r
                    ∆v = velocity_vec[i] .- velocity_vec[j]
                    ω = 1 - r / args.cutoff 
                    f_c = DPD.a * ω # Conservative force.
                    f_d = - DPD.γ * ω^2 * dot(r_unit, ∆v) # Dissipative force.
                    f_r = σ * ω * randn() * invsqrt_dt # Random force.
                    #f_r = (2 * rand() - 1) * σ * ω * r_unit # Random forcse.

                    f = (f_c + f_d + f_r) * r_unit

                    force_vec[i][1] += f[1] 
                    force_vec[j][1] -= f[1]
                    force_vec[i][2] += f[2]
                    force_vec[j][2] -= f[2]
                    force_vec[i][3] += f[3]
                    force_vec[j][3] -= f[3]
                    energy += 0.5 * DPD.a * ω^2
                    # xx, yy, zz, xy, xz, yz
                    stress[1] += f[1] * ∆r[1]
                    stress[2] += f[2] * ∆r[2]
                    stress[3] += f[3] * ∆r[3]
                    stress[4] += f[1] * ∆r[2]
                    stress[5] += f[1] * ∆r[3]
                    stress[6] += f[2] * ∆r[3]
                end    
            end
        end
    end
    return energy, stress

end