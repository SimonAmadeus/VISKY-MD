export System, Particle, Mesh, sys_init

# Defining data structures.
abstract type System_Type end

abstract type Thermostat end
abstract type Barostat end

abstract type Bonded_Interactions end
abstract type Non_Bonded_Interactions end
abstract type Angle_Interactions end

abstract type NL end # Neighbor list.

abstract type Field_Sum end

abstract type TrajectoryDump end
abstract type Logger end

# Special functions.
abstract type Shear end # Shear flow.
abstract type T_Bonds end # Slip-springs.
abstract type Collisions end # Collision dynamics.

# Structs and empty structs.
struct Molecular <: System_Type end

struct No_Mol <: System_Type end

struct No_Thermostat <: Thermostat end
struct No_Barostat <: Barostat end
    
struct No_Bonds <: Bonded_Interactions end
struct No_Non_Bonded <: Non_Bonded_Interactions end
struct No_Angles <: Angle_Interactions end

struct No_NL <: NL end

struct No_Sum <: Field_Sum end
    
struct No_Logger <: Logger end
struct No_TrajectoryDump <: TrajectoryDump end

struct No_Shear <: Shear end
struct No_T_Bonds <: T_Bonds end
struct No_Collisions <: Collisions end

# particle information, such as positions, velocities, forces, etc.
mutable struct Particle
    i_mol::Int64
    spec::Int64
    mass::Float64
    charge::Float64
    pos::Vector{Float64}
    image::Vector{Int64}
end
        
struct Bond
    type::Int64
    p_1::Int64
    p_2::Int64
end

mutable struct T_Bond # Temporary bond (Slip-Spring).
    type::Int64
    p_1::Int64
    p_2::Int64
end

struct Angle
    type::Int64
    p_1::Int64
    p_2::Int64
    p_3::Int64
end

mutable struct Neighbors
   list::Any 
end

# mesh used for hPF
struct Mesh
    edge_size::Array{Float64,1}
    NC::Array{Int64,1}
    NC_total::Int64
    NV_total::Int64
    particle_grid::Array{Float64,1}
    gradient_grid::Array{Float64,2}
    boxsize::Array{Float64,1}
    use_staggered::Bool
end

# System structure. Specify control information.
mutable struct System
    system_type::System_Type

    cutoff::Float64
    velocity_rescale::Int64

    # Special features:
    shear::Shear # Shear flow control.
    shear_sample_period::Int64
    j::Float64 # Exchanged momentum. Summation in shear.jl.

    t_bonds::T_Bonds # Slip-springs.

    collisions::Collisions # Collision dynamics.

    bonded_interactions::Bonded_Interactions
    non_bonded_interactions::Non_Bonded_Interactions
    angle_interactions::Angle_Interactions

    nl::NL # Neighbor list.

    sum::Field_Sum

    thermostat::Thermostat
    barostat::Barostat

    init_velocity::AbstractString
    
    ∆t::Float64 # Timestep.
    T::Float64
    p::Float64
    κ::Float64 # Compressibilitty of the system.

    box::Vector{Float64}
    
    n_steps::Int64
    traj_period::Int64
    log_period::Int64

    data_file::AbstractString
    traj_file::AbstractString
    log_file::AbstractString
    shear_file::AbstractString
    momentum_file::AbstractString
    restart_file::AbstractString
    
    trajectorydump::TrajectoryDump
    logger::Logger
    restart::Bool

    first_step::Bool
end

# Initialize system. 
function sys_init()
    system_type = No_Mol()

    cutoff = 2.5
    velocity_rescale = 0

    shear = No_Shear()
    shear_sample_period = 0
    j = 0 

    t_bonds = No_T_Bonds()

    collisions = No_Collisions()

    bonded_interactions = No_Bonds()
    non_bonded_interactions = No_Non_Bonded()
    angle_interactions = No_Angles()

    nl = No_NL()

    sum = No_Sum()

    thermostat = No_Thermostat()
    barostat = No_Barostat()

    init_velocity = "Z"

    ∆t = 0.0 
    T = 0.0
    p = 0.0
    κ = 0.0 

    box = [0.0, 0.0, 0.0]

    n_steps = 0
    traj_period = 0
    log_period = 0 
    
    data_file = "input.data"
    traj_file = "trj.xyz"
    log_file = "logger_MD.txt"
    shear_file = "velocity_profile.txt"
    momentum_file = "momentum.txt"
    restart_file = "restart.data"
    
    trajectorydump = No_TrajectoryDump()
    logger = No_Logger()
    restart = false

    first_step = true

    return System(system_type, cutoff, velocity_rescale, shear, shear_sample_period, j, t_bonds, collisions, bonded_interactions, non_bonded_interactions, angle_interactions, nl, sum, thermostat, barostat, init_velocity, ∆t, T, p, κ, box, n_steps, traj_period, log_period, data_file, traj_file, log_file, shear_file, momentum_file, restart_file, trajectorydump, logger, restart, first_step)
end
