export Berendsen_Thermostat, Andersen_Thermostat, Lowe_Andersen_Thermostat, Langevin_Thermostat, Nose_Hoover_Thermostat

struct Berendsen_Thermostat <: Thermostat 
    τ::Float64
end

struct Andersen_Thermostat <: Thermostat 
    σ::Float64
end

struct Lowe_Andersen_Thermostat <: Thermostat
    Γ::Float64
    r_c::Float64
end

struct Langevin_Thermostat <: Thermostat
    γ::Float64
end

mutable struct Nose_Hoover_Thermostat <: Thermostat
    #ζ::Float64
    Q::Float64
    τ::Float64
end

# no thermostat.
function apply_thermostat!(args::System, particle_vec::Vector{Particle}, velocity_vec, neihbor_list::Any, ::No_Thermostat) end

# Berendsen
function apply_thermostat!(args::System, particle_vec::Vector{Particle}, velocity_vec, neihbor_list::Any, Thermostat::Berendsen_Thermostat)
    # Determine scaling factor.
    K_tot = calc_total_kinetic_energy(particle_vec, velocity_vec)
    curr_T = calc_inst_temp(particle_vec, K_tot)
    scalT  = sqrt(1.0 + args.∆t * (args.T / curr_T - 1.0) / Thermostat.τ)
    # Update velocities.
    if curr_T != 0
        for i = 1:length(particle_vec)
            velocity_vec[i][1] *= scalT
            velocity_vec[i][2] *= scalT
            velocity_vec[i][3] *= scalT
        end
    end
end

# Andersen
function apply_thermostat!(args::System, particle_vec::Vector{Particle}, velocity_vec, neihbor_list::Any, Thermostat::Andersen_Thermostat)
    # Update velocities.
    for i = 1:length(particle_vec)
        rn = rand(1)
        if rn[1] < Thermostat.σ * args.∆t
            d = Normal(0.0, sqrt(args.T / particle_vec[i].mass))
            velocity_vec[i] = rand(d, 3)
        end
    end

end

# Lowe-Andersen
function apply_thermostat!(args::System, particle_vec::Vector{Particle}, velocity_vec, neighbor_list::Any, Thermostat::Lowe_Andersen_Thermostat)
    # Update velocities in a pairwise fashion.
    for i = 1:length(neighbor_list)
        for j = 1:length(neighbor_list[i])
            rn = rand(1)
            if rn[1] < Thermostat.Γ * args.∆t
                d = Normal(0.0, sqrt(2 * args.T / particle_vec[i].mass))
                ∆d = particle_vec[i].pos .- particle_vec[neighbor_list[i][j]].pos
                ∆d .-= args.box .* floor.((∆d .+ args.box / 2) ./ args.box) # PBCs
                r = norm(∆d)
                if r < Thermostat.r_c
                    r_unit = ∆d / r
                    ∆v = velocity_vec[i] .- velocity_vec[neighbor_list[i][j]] 
                    v_random = rand(d, 3)
                    ∆v_relative = v_random - ∆v
                    ∆ij = 0.5 * r_unit * dot(∆v_relative, r_unit)
                    velocity_vec[i] += ∆ij
                    velocity_vec[neighbor_list[i][j]] -= ∆ij
                end
            end    
        end
    end
end

# Langevin
function apply_thermostat!(args::System, particle_vec::Vector{Particle}, velocity_vec, neighbor_list::Any, thermostat::Langevin_Thermostat)
#=    γ = thermostat.γ
    kBT = args.T
    dt = args.∆t

    for i = 1:length(particle_vec)
        m = particle_vec[i].mass
        σ = sqrt(2 * γ * kBT * dt / m)
        f_r = σ * randn(3)  # Random force for thermal fluctuations

        # First half of velocity Verlet integration for Langevin
        velocity_vec[i].velocity += dt * (-γ * velocity_vec[i].velocity + f_r) / m
    end
=#
end

function update_velocities_langevin!(args, particle_vec, velocity_vec, force_vec, thermostat::Langevin_Thermostat)
    γ = thermostat.γ
    kBT = args.T
    dt = args.∆t

    for i in 1:length(particle_vec)
        m = particle_vec[i].mass
        σ = sqrt(2 * γ * kBT / dt)
        #σ = sqrt(γ * kBT / dt)
        f_r = σ * randn(3)

        velocity_vec[i] .+= dt / m * (force_vec[i] - γ * velocity_vec[i] + f_r)
        #velocity_vec[i] .+= 0.5 * dt / m * (force_vec[i] - γ * velocity_vec[i] + f_r)
    end
end

# Nose-Hoover
function apply_thermostat!(args::System, particle_vec::Vector{Particle}, velocity_vec, neighbor_list::Any, thermostat::Nose_Hoover_Thermostat)
    ∆t = args.∆t
    N = length(particle_vec)
    T_0 = args.T
    # total kinetic energy
    K_tot = calc_total_kinetic_energy(particle_vec, velocity_vec)
    # current temperature
    T_current = calc_inst_temp(particle_vec, K_tot)
    # thermostating
    #Q = 400
    #Q_dot = ((K_tot / 3) / (N * T_0) - 1) / thermostat.τ^2
    #ζ_dot = (T_current / T_0 - 1) / (thermostat.τ^2 * Q)
    #thermostat.Q += ζ_dot * ∆t
    #η = exp(-∆t * thermostat.Q)
    #thermostat.Q += ∆t * (curr_T / args.T - 1) / thermostat.Q
    #thermostat.Q += ∆t * Q_dot
    #println("Q: ", thermostat.Q)
    η = sqrt(1 + (∆t / thermostat.Q) * (T_0 / T_current - 1))

    for i in 1:length(velocity_vec)
        #velocity_vec[i].velocity *= exp(- ∆t * 0.5 * thermostat.ζ)
        velocity_vec[i] *= η
    end
    Q_dot = ((K_tot / 3) / (N * T_0) - 1) / thermostat.τ^2
    thermostat.Q += ∆t * Q_dot
end