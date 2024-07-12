export LJ_12_6

struct LJ_12_6 <: Non_Bonded_Interactions
    ε::Float64
    σ::Float64
end

# Lennard-Jones interactions for monoatomic systems.
function apply_non_bonded_interactions!(args::System, ::No_Mol, particle_vec::Vector{Particle}, velocity_vec, neighbor_list::Any, c_l::Int64, force_vec, energy::Float64, stress::Vector{Float64}, mesh::Mesh, vertex::Array{Int64, 3}, cell_vertices::Array{Int64, 2}, lj126::LJ_12_6)

    #for i = 2:length(particle_vec)
    #    for j = 1:i-1
    for i = 1:length(neighbor_list)
        for k = 1:length(neighbor_list[i])
            j = neighbor_list[i][k]
            if j > i
                ∆r = particle_vec[i].pos .- particle_vec[j].pos
                ∆r .-= args.box .* floor.((∆r .+ args.box / 2) ./ args.box) # PBCs
                r = norm(∆r)
                r6i = 0
                if r < args.cutoff
                    r6i = 1.0 / r^6
                    f = 48 * lj126.ε * (lj126.σ * r6i^2 - 0.5 * lj126.σ * r6i)
                    f_x = f * ∆r[1] / r
                    f_y = f * ∆r[2] / r
                    f_z = f * ∆r[3] / r
                    force_vec[i][1] += f_x
                    force_vec[j][1] -= f_x
                    force_vec[i][2] += f_y
                    force_vec[j][2] -= f_y
                    force_vec[i][3] += f_z
                    force_vec[j][3] -= f_z
                    energy += 4 * lj126.ε * (lj126.σ * r6i^2 - lj126.σ * r6i)
                    # xx, yy, zz, xy, xz, yz
                    stress[1] += f_x * ∆r[1]
                    stress[2] += f_y * ∆r[2]
                    stress[3] += f_z * ∆r[3]
                    stress[4] += f_x * ∆r[2]
                    stress[5] += f_x * ∆r[3]
                    stress[6] += f_y * ∆r[3]
                end
            end
        end
    end
    return energy, stress
end

# 1-2 and 1-3 exclusion for bonded atoms.
function apply_non_bonded_interactions!(args::System, ::Molecular, particle_vec::Vector{Particle}, velocity_vec, neighbor_list::Any, c_l::Int64, force_vec, energy::Float64, stress::Vector{Float64}, mesh::Mesh, lj126::LJ_12_6)

    #for i = 2:length(particle_vec)
    #    for j = 1:i-1
    for i = 1:length(neighbor_list)
        for k = 1:length(neighbor_list[i])
            j = neighbor_list[i][k]
            if j > i
                # Check if 2 beads are neighbors in the same polymer chain:
                if abs(i - j) > 2 || ( mod(i, c_l) == 1 && (j - i == -1 || j - i == - 2)) || ( mod(j, c_l) == 1 && (i - j == -1 || i - j == - 2)) || ( mod(i, c_l) == 2 && j - i == - 2 ) || ( mod(j, c_l) == 2 && i - j == - 2 ) || ( mod(i, c_l) == 0 && (j - i == 1 || j - i == 2)) || ( mod(j, c_l) == 0 && (i - j == 1 || i - j == 2)) || ( mod(i, c_l) == c_l - 1 && j - i == 2) || ( mod(j, c_l) == c_l - 1 && i - j == 2)
                    ∆r = particle_vec[i].pos .- particle_vec[j].pos
                    ∆r .-= args.box .* floor.((∆r .+ args.box / 2) ./ args.box) # PBCs
                    r = norm(∆r)
                    r6i = 0
                    if r < args.cutoff
                        r6i = 1.0 / r^6
                        f = 48 * lj126.ε * (lj126.σ * r6i^2 - 0.5 * lj126.σ * r6i)
                        f_x = f * ∆r[1] / r
                        f_y = f * ∆r[2] / r
                        f_z = f * ∆r[3] / r
                        force_vec[i][1] += f_x
                        force_vec[j][1] -= f_x
                        force_vec[i][2] += f_y
                        force_vec[j][2] -= f_y
                        force_vec[i][3] += f_z
                        force_vec[j][3] -= f_z
                        energy += 4 * lj126.ε * (lj126.σ * r6i^2 - lj126.σ * r6i)
                        # xx, yy, zz, xy, xz, yz
                        stress[1] += f_x * ∆r[1]
                        stress[2] += f_y * ∆r[2]
                        stress[3] += f_z * ∆r[3]
                        stress[4] += f_x * ∆r[2]
                        stress[5] += f_x * ∆r[3]
                        stress[6] += f_y * ∆r[3]
                    end
                end
            end
        end
    end
    return energy, stress
end
