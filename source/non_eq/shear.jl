export Shear_CTRL

struct Shear_CTRL <: Shear
    mode::AbstractString
    n_exc::Int64
    n_slabs::Int64
    freq::Int64
end

function apply_shear!(args::System, particle_vec::Vector{Particle}, velocity_vec, ::No_Shear) end

# Exchange of velocities takes place
function apply_shear!(args::System, particle_vec::Vector{Particle}, velocity_vec, shear::Shear_CTRL)
    
    i_M = floor(shear.n_slabs / 2) + 1 # Number of middle slab.
    d_slab = args.box[3] / shear.n_slabs # Thickness of slabs in z coordinate.

    i_list_0 = [] #::Vector{Int64} # List of velocity indices in slab 0.
    i_list_M = [] #::Vector{Int64} # List of velocity indices in slab M.
    vel_0 = [] #::Vector{Float64} # List of velocities in slab 0.
    vel_M = [] #::Vector{Float64} # List of velocities in slab M.
    
    # Sort particle indices into slabs.
    for i in 1:length(particle_vec)

        z = particle_vec[i].pos[3] + args.box[3] / 2
        
        if z <= d_slab
            append!(i_list_0, i)
        elseif z > (i_M - 1) * d_slab && z <= i_M * d_slab
            append!(i_list_M, i)       
        end
    
    end
    
    if isempty(i_list_0) == false && isempty(i_list_M) == false

        # Exchange velocity_vec.
        for i in i_list_0
            append!(vel_0, velocity_vec[i][1])
        end
        for i in i_list_M
            append!(vel_M, velocity_vec[i][1])
        end

        v_0 = findmin(vel_0)
        v_M = findmax(vel_M)
        i_0 = v_0[2]
        i_M = v_M[2]

        velocity_vec[i_list_0[i_0]][1] = v_M[1]
        velocity_vec[i_list_M[i_M]][1] = v_0[1]

        args.j += v_M[1] * particle_vec[i_list_0[i_0]].mass -  v_0[1] * particle_vec[i_list_M[i_M]].mass

    end

    i_list_0 = []
    i_list_M = []
    vel_0 = []
    vel_M = []

end

# create bins for the velocity profile
function initialize_shear!(args::System, shear::Shear_CTRL)
    d_slab = args.box[3] / shear.n_slabs
    z_bins = [] # Storage for v_x.

    for i in 1:shear.n_slabs
        push!(z_bins, [])
    end
    return d_slab, z_bins
end

function sample_velocity_profile!(args::System, particle_vec::Vector{Particle}, velocity_vec, d_slab::Float64, z_bins::Vector{}, shear::Shear_CTRL)    
# Sample velocity profile in a single simulation step
    for i in 1:length(particle_vec)
        z = particle_vec[i].pos[3] + args.box[3] / 2
        i_bin = floor(Int64, z / d_slab) + 1
        push!(z_bins[i_bin], velocity_vec[i][1])
    end
end

# output
function apply_shear_log(file0, file1, step::Int64, args::System, d_slab, z_bins, shear::Shear_CTRL)

    string_out = "z v_x\n" # Header.
    write(file0, string_out)

    for i in 1:shear.n_slabs
        string_out = ""
        if isempty(z_bins[i]) == false
            string_out = string(i * d_slab) * " " * string(mean(z_bins[i]))
        else
            string_out = string(i * d_slab) * " 0"    
        end
        string_out *= "\n"
        write(file0, string_out)
    end

    args.j = args.j / (2 * step * args.∆t * args.box[1] * args.box[2])

    string_out = "Total momentum exchange j: $(args.j)\n \n With j = p_x / (2 * t * A)" # Header.
    write(file1, string_out)
end
