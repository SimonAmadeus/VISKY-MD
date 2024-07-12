export Force_Sum

# Determine properties fieldwise, like force per cell, e.g.

struct Force_Sum <: Field_Sum
    a::Float64
    freq::Int64
    massage::AbstractString
end

function initialize_force_cells(args::System, a::Float64)
    nx = Int64(div(args.box[1], a))
    ny = Int64(div(args.box[2], a))
    nz = Int64(div(args.box[3], a))

    force_sum_x = zeros(Float64, nx, ny, nz)
    force_sum_y = zeros(Float64, nx, ny, nz)
    force_sum_z = zeros(Float64, nx, ny, nz)
    return force_sum_x, force_sum_y, force_sum_z
end

function sum_force!(args::System, particle_vec::Vector{Particle}, force_vec, force_sum_x::Array{Float64, 3}, force_sum_y::Array{Float64, 3}, force_sum_z::Array{Float64, 3}, a::Float64)
    for i in 1:length(particle_vec)
        ix = Int64(div(particle_vec[i].pos[1] + args.box[1] / 2, a)) + 1
        iy = Int64(div(particle_vec[i].pos[2] + args.box[2] / 2, a)) + 1
        iz = Int64(div(particle_vec[i].pos[3] + args.box[3] / 2, a)) + 1

        force_sum_x[ix, iy, iz] += force_vec[i][1]
        force_sum_y[ix, iy, iz] += force_vec[i][2]
        force_sum_z[ix, iy, iz] += force_vec[i][3] 
    end
end

function clear_sums!(force_sum_x::Array{Float64, 3}, force_sum_y::Array{Float64, 3}, force_sum_z::Array{Float64, 3})
    for i in 1:size(force_sum_x, 1)
        for j in 1:size(force_sum_x, 2)
            for k in 1:size(force_sum_x, 3)
                force_sum_x[i,j,k] = 0
                force_sum_y[i,j,k] = 0
                force_sum_z[i,j,k] = 0
            end
        end
    end
end

function force_massaging!(args::System, particle_vec::Vector{Particle}, force_vec, force_sum_x::Array{Float64, 3}, force_sum_y::Array{Float64, 3}, force_sum_z::Array{Float64, 3}, a::Float64)
    nx = Int64(div(args.box[1], a))
    ny = Int64(div(args.box[2], a))
    nz = Int64(div(args.box[3], a))

    p_counter = zeros(Int64, nx, ny, nz)

    for i in 1:length(particle_vec)
        ix = Int64(div(particle_vec[i].pos[1] + args.box[1] / 2, a)) + 1
        iy = Int64(div(particle_vec[i].pos[2] + args.box[2] / 2, a)) + 1
        iz = Int64(div(particle_vec[i].pos[3] + args.box[3] / 2, a)) + 1

        p_counter[ix, iy, iz] += 1

        force_sum_x[ix, iy, iz] += force_vec[i][1]
        force_sum_y[ix, iy, iz] += force_vec[i][2]
        force_sum_z[ix, iy, iz] += force_vec[i][3] 
    end

    for i in 1:length(force_vec)
        ix = Int64(div(particle_vec[i].pos[1] + args.box[1] / 2, a)) + 1
        iy = Int64(div(particle_vec[i].pos[2] + args.box[2] / 2, a)) + 1
        iz = Int64(div(particle_vec[i].pos[3] + args.box[3] / 2, a)) + 1

        force_vec[i][1] -= force_sum_x[ix, iy, iz] / p_counter[ix, iy, iz]
        force_vec[i][2] -= force_sum_y[ix, iy, iz] / p_counter[ix, iy, iz]
        force_vec[i][3] -= force_sum_z[ix, iy, iz] / p_counter[ix, iy, iz]
    end
end

function print_force_info(step_i::Int64, force_sum_x::Array{Float64, 3}, force_sum_y::Array{Float64, 3}, force_sum_z::Array{Float64, 3})
    println("Step $(step_i)")
    println(force_sum_x)
    println(std(force_sum_x))
    println(mean(force_sum_x))
    println(sum(force_sum_x))
    println(force_sum_y)
    println(std(force_sum_y))
    println(mean(force_sum_x))
    println(sum(force_sum_y))
    println(force_sum_z)
    println(std(force_sum_z))
    println(mean(force_sum_x))
    println(sum(force_sum_z))
    println("Average: ")
    av_std_sum += (std(force_sum_x) + std(force_sum_x) + std(force_sum_z)) / (3 * (args.n_steps / args.sum.freq))
    av_sum += (sum(force_sum_x) + sum(force_sum_x) + sum(force_sum_z)) / (3 * (args.n_steps / args.sum.freq))
    println(av_std_sum)
    println(av_sum)
end