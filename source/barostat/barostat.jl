export Berendsen_Barostat

struct Berendsen_Barostat <: Barostat
    τ_p::Float64
end

function apply_barostat!(args::System, particle_vec::Vector{Particle}, curr_p::Float64, ::No_Barostat) end

function apply_barostat!(args::System, particle_vec::Vector{Particle}, curr_p::Float64, barostat::Berendsen_Barostat)
    
    µ = 1 - args.∆t * args.κ / (3 * barostat.τ_p) * (args.p - curr_p)
    args.box *= µ
    for i in 1:length(particle_vec)
        particle_vec[i].pos *= µ
        args.box *= µ
    end

end
