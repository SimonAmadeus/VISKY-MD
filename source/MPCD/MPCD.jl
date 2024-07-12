export MPCD

mutable struct MPCD <: Collisions
    #a::Float64 # Cell size.
    n_cells::Vector{Int64} # number of cells.
    δt::Int64 # Time between collisions.
    θ::Float64 # Angle in degrees.
end

function apply_collisions!(args::System, particle_vec::Vector{Particle}, velocity_vec, collisions::MPCD)
    
    len_cells = args.box ./ collisions.n_cells
    nx = collisions.n_cells[1]
    ny = collisions.n_cells[2]
    nz = collisions.n_cells[3]
    #nx = Int64(div(args.box[1], collisions.a))
    #ny = Int64(div(args.box[2], collisions.a))
    #nz = Int64(div(args.box[3], collisions.a))
    p_counter = zeros(Int64, nx, ny, nz)
    v_cm = zeros(Float64, nx, ny, nz, 3)
    E_kin_cell = zeros(Float64, nx, ny, nz)
    random_axis = Array{Array{Float64, 1}}(undef, nx, ny, nz)
    #random_axis = zeros(Float64, nx, ny, nz, 3)
    step_size = 1e-6
    shift = zeros(Float64, 3)
    #shift = rand(((- collisions.a / 2) : step_size : (collisions.a / 2)), 3)
    shift[1] = rand(((- len_cells[1] / 2) : step_size : (len_cells[1] / 2)))
    shift[2] = rand(((- len_cells[2] / 2) : step_size : (len_cells[2] / 2)))
    shift[3] = rand(((- len_cells[3] / 2) : step_size : (len_cells[3] / 2)))
    x = similar(particle_vec, Float64) 
    y = similar(particle_vec, Float64)
    z = similar(particle_vec, Float64)
    for i in eachindex(particle_vec)
        x[i] = mod(particle_vec[i].pos[1] + args.box[1] / 2 + shift[1], args.box[1])
        y[i] = mod(particle_vec[i].pos[2] + args.box[2] / 2 + shift[2], args.box[2])
        z[i] = mod(particle_vec[i].pos[3] + args.box[3] / 2 + shift[3], args.box[3])
    end
    for i in eachindex(particle_vec)
        #ix = Int64(div(x[i] + args.box[1] / 2, collisions.a)) + 1
        #iy = Int64(div(y[i] + args.box[2] / 2, collisions.a)) + 1
        #iz = Int64(div(z[i] + args.box[3] / 2, collisions.a)) + 1
        ix = Int64(div(x[i], len_cells[1])) + 1
        iy = Int64(div(y[i], len_cells[2])) + 1
        iz = Int64(div(z[i], len_cells[3])) + 1

        p_counter[ix, iy, iz] += 1

        # Accumulate velocity of particles in cell.
        v_cm[ix, iy, iz, 1] += velocity_vec[i][1]
        v_cm[ix, iy, iz, 2] += velocity_vec[i][2]
        v_cm[ix, iy, iz, 3] += velocity_vec[i][3]
        E_kin_cell[ix, iy, iz] += 0.5 * particle_vec[i].mass * dot(velocity_vec[i], velocity_vec[i])
    end

    # Compute center of mass velocity for each cell (average).
    for i in 1:nx
        for j in 1:ny
            for k in 1:nz
                if p_counter[i, j, k] > 0
                    v_cm[i, j, k, 1] /= p_counter[i, j, k]
                    v_cm[i, j, k, 2] /= p_counter[i, j, k]
                    v_cm[i, j, k, 3] /= p_counter[i, j, k]
                end
            end
        end
    end

    # Create random axis for each cell.
    for i in 1:nx
        for j in 1:ny
            for k in 1:nz
                random_axis[i, j, k] = normalize(Float64.(rand(3)))
            end
        end
    end

    #θ = rand(Uniform(30, 180))
    #θ = rand(Normal(130, 10))
    #θ = generate_sine_angle()
    θ = collisions.θ

    # velocity update in collision step
    for i in eachindex(particle_vec)
        ix = Int64(div(x[i], len_cells[1])) + 1
        iy = Int64(div(y[i], len_cells[2])) + 1
        iz = Int64(div(z[i], len_cells[3])) + 1

        # update particle velocity
        v_cm_cell = [v_cm[ix, iy, iz, 1], v_cm[ix, iy, iz, 2], v_cm[ix, iy, iz, 3]]
        v_rel = velocity_vec[i] .- v_cm_cell
        R = generate_rotation_matrix(random_axis[ix, iy, iz], θ)

        velocity_vec[i][1] = v_cm_cell[1] + R[1, 1] * v_rel[1] + R[1, 2] * v_rel[2] + R[1, 3] * v_rel[3]
        velocity_vec[i][2] = v_cm_cell[2] + R[2, 1] * v_rel[1] + R[2, 2] * v_rel[2] + R[2, 3] * v_rel[3]
        velocity_vec[i][3] = v_cm_cell[3] + R[3, 1] * v_rel[1] + R[3, 2] * v_rel[2] + R[3, 3] * v_rel[3]
    end

    # Thermostat:
    #apply_MBS_thermostat!(args, particle_vec, shift, velocity_vec, v_cm, E_kin_cell, collisions)

end

function apply_MBS_thermostat!(args::System, particle_vec::Vector{Particle}, shift::Vector{Float64}, velocity_vec, v_cm::Array{Float64, 4}, E_kin_cell::Array{Float64}, collisions::MPCD)

    total_velocity = [0.0, 0.0, 0.0]
    len_cells = args.box ./ collisions.n_cells
    
    @inbounds @simd for i = 1:length(particle_vec)
        # position shift  
        x = mod(particle_vec[i].pos[1] + args.box[1] / 2 + shift[1], args.box[1])
        y = mod(particle_vec[i].pos[2] + args.box[2] / 2 + shift[2], args.box[2])
        z = mod(particle_vec[i].pos[3] + args.box[3] / 2 + shift[3], args.box[3])
    
        ix = Int64(div(x, len_cells[1])) + 1
        iy = Int64(div(y, len_cells[2])) + 1
        iz = Int64(div(z, len_cells[3])) + 1

        v_cm_cell = [v_cm[ix, iy, iz, 1], v_cm[ix, iy, iz, 2], v_cm[ix, iy, iz, 3]]
        v_rel = velocity_vec[i] .- v_cm_cell
        
        # create a normal distribution with mean 0 and standard deviation sqrt(k_B*T/m)
        d = Normal(0.0, sqrt(args.T / particle_vec[i].mass))

        # draw velocities from Maxwell-Boltzmann distribution and calculate kinetic energy
        v_sample = [rand(d) for _ in 1:3]
        E_k_target = 0.5 * particle_vec[i].mass * sum(v_sample .^ 2)

        ξ = sqrt(E_k_target / E_kin_cell[ix, iy, iz])

        # update particle velocity according to MBS thermostat
        velocity_vec[i] = v_cm_cell .+ ξ * v_rel 
        total_velocity .+= velocity_vec[i]
    end
    
    v_cm_total = total_velocity ./ length(particle_vec)

    @inbounds @simd for i = 1:length(particle_vec)
        velocity_vec[i] .-= v_cm_total
    end
end

# generate a rotation matrix around a random axis with a given angle 
# Rodrigues' rotation formula
function generate_rotation_matrix(axis::Vector{Float64}, angle::Float64)

    u = normalize(axis)
    cos_angle = cos(angle)
    sin_angle = sin(angle)
    one_minus_cos_angle = 1 - cos_angle

    R = Matrix{Float64}(LinearAlgebra.I, 3, 3)
    R[1, 1] = cos_angle + u[1]^2 * one_minus_cos_angle
    R[1, 2] = u[1] * u[2] * one_minus_cos_angle - u[3] * sin_angle
    R[1, 3] = u[1] * u[3] * one_minus_cos_angle + u[2] * sin_angle
    R[2, 1] = u[2] * u[1] * one_minus_cos_angle + u[3] * sin_angle
    R[2, 2] = cos_angle + u[2]^2 * one_minus_cos_angle
    R[2, 3] = u[2] * u[3] * one_minus_cos_angle - u[1] * sin_angle
    R[3, 1] = u[3] * u[1] * one_minus_cos_angle - u[2] * sin_angle
    R[3, 2] = u[3] * u[2] * one_minus_cos_angle + u[1] * sin_angle
    R[3, 3] = cos_angle + u[3]^2 * one_minus_cos_angle

    return R
end

# function to generate a random angle with sine distribution
function generate_sine_angle()
    u = rand()  # generate a uniform random number [0, 1]
    angle = acos(1 - 2 * u)  # apply the inverse cumulative distribution function
    return rad2deg(angle)  # convert to degrees
end
