export Harmonic_Bonded, FENE_Bonded

struct Harmonic_Bonded <: Bonded_Interactions 
    κ::Float64
    r0::Float64
end

struct FENE_Bonded <: Bonded_Interactions
    K::Float64 #    (energy/distance^2)
    R0::Float64 #   (distance)
    ϵ::Float64 #    (energy)
    σ::Float64 #    (distance)
end 

function apply_bonded_interactions!(args::System, particle_vec::Vector{Particle}, force_vec, energy::Float64, stress::Vector{Float64}, bond_vec::Vector{Bond}, ::No_Bonds)
    return energy, stress                                                                                                                               
end

function apply_bonded_interactions!(args::System, particle_vec::Vector{Particle}, force_vec, energy::Float64, stress::Vector{Float64}, bond_vec::Vector{Bond}, harmonic_bond::Harmonic_Bonded) 
    
    # 0.5×k×(r-r0)^2
    for i = 1:length(bond_vec)
        pos_1 = particle_vec[bond_vec[i].p_1].pos .+ particle_vec[bond_vec[i].p_1].image .* args.box
        pos_2 = particle_vec[bond_vec[i].p_2].pos .+ particle_vec[bond_vec[i].p_2].image .* args.box
        ∆d = pos_2 .- pos_1
        #∆d .-= args.box .* floor.((∆d + args.box / 2)./ args.box) # PBCs.
        r = norm(∆d)
        f = - harmonic_bond.κ * (r - harmonic_bond.r0) 
        force_divr = f / r
        force_vec[bond_vec[i].p_1].force .-= ∆d .* force_divr 
        force_vec[bond_vec[i].p_2].force .+= ∆d .* force_divr
        energy += 0.5 * harmonic_bond.κ * (r - harmonic_bond.r0)^2
        # xx, yy, zz, xy, xz, yz
        stress[1] += force_vec[bond_vec[i].p_1].force[1] * ∆d[1] 
        stress[2] += force_vec[bond_vec[i].p_1].force[2] * ∆d[2] 
        stress[3] += force_vec[bond_vec[i].p_1].force[3] * ∆d[3] 
        stress[4] += force_vec[bond_vec[i].p_1].force[1] * ∆d[2] 
        stress[5] += force_vec[bond_vec[i].p_1].force[1] * ∆d[3] 
        stress[6] += force_vec[bond_vec[i].p_1].force[2] * ∆d[3] 
    end

    return energy, stress

end

function apply_bonded_interactions!(args::System, particle_vec::Vector{Particle}, force_vec, energy::Float64, stress::Vector{Float64}, bond_vec::Vector{Bond}, fenebond::FENE_Bonded) 

    for i = 1:length(bond_vec)

        pos_1 = particle_vec[bond_vec[i].p_2].pos .+ particle_vec[bond_vec[i].p_2].image .* args.box
        pos_2 = particle_vec[bond_vec[i].p_1].pos .+ particle_vec[bond_vec[i].p_1].image .* args.box
        ∆d = pos_2 .- pos_1

        r = norm(∆d)
        rsq = r^2

        r0sq = fenebond.R0^2
        rlogarg = 1 - rsq / r0sq

        fbond = - fenebond.K / rlogarg

        if rsq < fenebond.σ^2
            sr2 = fenebond.σ^2 / rsq
            sr6 = sr2*sr2*sr2
            force_divr += 48*fenebond.ϵ * sr6 * (sr6 - 0.5) / rsq
        end

        E = -0.5*fenebond.K*fenebond.R0^2*log(1-(r/fenebond.R0)^2)+4*fenebond.ϵ*((fenebond.σ/r)^12-(fenebond.σ/r)^6)+fenebond.ϵ

        energy += E

        # calculate forces
        force_vec[bond_vec[i].p_1+1].force .+= delta_d .* force_divr
        force_vec[bond_vec[i].p_2+1].force .+= -delta_d .* force_divr

    end

    return energy, stress

end
