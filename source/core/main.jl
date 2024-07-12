export run_MD

include("system.jl")

# Main simulation loop.
# In this script parameters and functionalities are also initianized.
#
# Parameters are:
# energy -> potential energy of the system
# virial -> virial of the system
# curr_p -> current pressure of the system
# stress -> stress tensor of the system
# E_kin -> kinetic energy tensor of the system
#
# Functionalities are:
# - particle field representation
# - Neighbor list (with skim)
# - slip-springs (for CG polymers)
# - shear flow (viscosity calculation)

function run_MD(args::System)
    
    print_info(args)

    # Parameters:
    energy = 0.0e0
    virial = 0.0e0
    curr_p = 0.0e0
    # Stress tensor: xx, yy, zz, xy, xz, yz.
    stress = Vector(zeros(6))
    # Kinetic energy: xx, yy, zz, xy, xz, yz.
    E_kin = Vector(zeros(6)) 
    
    # Dummy mesh, used in hPF:
    mesh = init_mesh(args.box, [10, 10, 10])
    
    # Declare system type:
    if typeof(args.bonded_interactions) != No_Bonds
        args.system_type = Molecular()
    end
    # Use neighbor list, if LAT is active.
    if typeof(args.thermostat) == Lowe_Andersen_Thermostat || typeof(args.non_bonded_interactions) == LJ_12_6 || typeof(args.non_bonded_interactions) == DPD
        #args.nl = standard_NL()
        args.nl = CL_NL()
    end
    # Set cutoff, if LAT is active.
    if typeof(args.thermostat) == Lowe_Andersen_Thermostat
        global args.cutoff = args.thermostat.r_c
    end
    if typeof(args.thermostat) == Nose_Hoover_Thermostat
        args.thermostat.Q = args.T / args.thermostat.τ^2
    end

################################################################################

    # Read data

    particle_vec, c_l, velocity_vec, bond_vec, angle_vec, args.box = read_lammps_x(args.data_file, args.system_type)

################################################################################

    N = length(particle_vec) # Total number of particles N.
    
    # For neighbor list:
    skim = 0.5
    skim2 = skim^2
    nl_cutoff = args.cutoff + skim
    # For cell list (x,y and z):
    n_cells = floor.(Int64, args.box ./ nl_cutoff)

    #c_l = Int(N / particle_vec[N].i_mol) # Chain length.
    
    # Slip-spring initialization.
    if typeof(args.t_bonds) != No_T_Bonds
        # Initialize without using restart data:
        if args.t_bonds.restart == "n"
            t_bond_vec = init_t_bonds!(args, c_l, particle_vec, bond_vec, args.t_bonds)
        # Read temporary bonds (slip-springs)v from input file.
        elseif args.t_bonds.restart == "y"
            t_bond_vec = read_t_bonds()
        else
            error("Only y or n as slip-spring restart argument are allowed.")
        end
    end

    # Shear initialization.
    if typeof(args.shear) != No_Shear
        d_slab, z_bins = initialize_shear!(args, args.shear)
    end

    # hPF initialization.
    if typeof(args.non_bonded_interactions) == original_hPF_Interactions
        mesh = init_mesh(args.box, args.non_bonded_interactions.mesh_points)
        args.non_bonded_interactions = original_hPF_Interactions(args.non_bonded_interactions.κ, args.non_bonded_interactions.mesh_points)
    #elseif typeof(args.non_bonded_interactions) == spectral_hPF_Interactions
    #    mesh = init_mesh(args.box, args.non_bonded_interactions.mesh_points)
    #    args.non_bonded_interactions = spectral_hPF_Interactions(args.non_bonded_interactions.κ, args.non_bonded_interactions.σ, args.non_bonded_interactions.mesh_points)
    end

    vertex = init_vertex_matrix(mesh)
    cell_vertices = init_cell_vertices(mesh, vertex)

    first_step = args.first_step

    # Open output files.
    log_file = open(args.log_file, "w")
    traj_file = open(args.traj_file, "w")
    # Shear flow:
    if typeof(args.shear) != No_Shear
        shear_file = open(args.shear_file, "w")
        momentum_file = open(args.momentum_file, "w")
    end
    # Slip-springs:
    if typeof(args.t_bonds) != No_T_Bonds
        entanglement_file = open("entanglement.out", "w")
    end
    # force/property output
    f_sample_steps = false
    f_sample_freq = 1
    if f_sample_steps == true
        force_file = open("forces.txt", "w")
    end
    # Calculate momentum drift
    calc_drift = false
    P_tot_prev = 0
    momentum_drift = []

################################################################################

    # First simulation step.

    if first_step == true

        # Initialize neighbor list.
        neighbor_list = [Vector{Int}() for _ in 1:length(particle_vec)]
        # Create neighbor list.
        if typeof(args.nl) != No_NL
            update_NL = false
            prev_positions = [copy(p.pos) for p in particle_vec]
            create_neighbor_list!(args, particle_vec, neighbor_list, n_cells, nl_cutoff, args.nl)
        end

        # Initialize veclocities.
        initialize_velocities!(args, particle_vec, velocity_vec)
        # Initialize force vector.
        #force_vec = Vector{Force}(undef, length(particle_vec))
        #for i = 1:length(particle_vec)
        #    force_vec[i] = Force(zeros(3))
        #end
        force_vec = [zeros(3) for _ in 1:length(particle_vec)]
        # Pressure calculation.
        E_kin = calc_kinetic_energy(particle_vec, velocity_vec, E_kin)
        curr_p = calc_pressure(stress, E_kin, args.box)
        # Logger.
        apply_logger(log_file, first_step, 0, velocity_vec, energy, curr_p, particle_vec, args.logger)
        apply_dump(traj_file, 0, args, particle_vec, force_vec, velocity_vec, args.trajectorydump)        
        clear_kin_energy_tensor!(E_kin)
        first_step = false
    end

################################################################################

################################################################################

    # Main simulation loop.

    println("The main simulation loop starts now.")

    @showprogress for step_i = 1:args.n_steps
    #for step_i = 1:args.n_steps

        # Neigbor list update:
        if typeof(args.nl) != No_NL
            # Neighbor list is only created, if args.nl != No_NL.
            for i in 1:length(particle_vec)
                dist_moved2 = sum((particle_vec[i].pos .- prev_positions[i]).^2)
                if dist_moved2 > skim2
                    update_NL = true
                    break
                end
            end
            # Neighbor list is only updated, if at least one particle traveled more than skim.
            if update_NL == true
                create_neighbor_list!(args, particle_vec, neighbor_list, n_cells, nl_cutoff, args.nl)
                prev_positions = [copy(p.pos) for p in particle_vec]
                update_NL = false
            end
        end

        # Shear flow (Momentum exchange). 
        if typeof(args.shear) != No_Shear && step_i%args.shear.freq == 0.0
            apply_shear!(args, particle_vec, velocity_vec, args.shear)
        end

        # Slip-spring chain end relocation and migration.
        if typeof(args.t_bonds) != No_T_Bonds && step_i%args.t_bonds.freq_mc == 0.0
            migration!(first_step, args, c_l, particle_vec, t_bond_vec, args.t_bonds)
        end
        
        # Apply first integration.
        if typeof(args.thermostat) == Langevin_Thermostat
            update_velocities_langevin!(args, particle_vec, velocity_vec, force_vec, args.thermostat)
        else
            update_velocities!(args, particle_vec, velocity_vec, force_vec)
        end

        # Multi particle collision dynamics (MPCD).
        if typeof(args.collisions) != No_Collisions && step_i%args.collisions.δt == 0.0
            apply_collisions!(args, particle_vec, velocity_vec, args.collisions)
        end

        #update_velocities!(args, particle_vec, velocity_vec, force_vec)
        update_positions!(args, particle_vec, velocity_vec, force_vec)
        
        clear_force!(force_vec)

        # Apply interactions.
        energy, stress = apply_bonded_interactions!(args, particle_vec, force_vec, energy, stress, bond_vec, args.bonded_interactions)
        energy1 = energy
        energy, stress = apply_angle_interactions!(args, particle_vec, force_vec, energy, stress, angle_vec, args.angle_interactions)
        energy2 = energy - energy1 
        #energy, stress = apply_non_bonded_interactions!(args, args.system_type, particle_vec, velocity_vec, neighbor_list, c_l, force_vec, energy, stress, mesh, vertex, cell_vertices, args.non_bonded_interactions)
        if typeof(args.non_bonded_interactions) != No_Non_Bonded
            #energy, stress = apply_non_bonded_interactions!(args, args.system_type, particle_vec, velocity_vec, neighbor_list, c_l, force_vec, energy, stress, mesh, vertex, cell_vertices, args.non_bonded_interactions)
            apply_non_bonded_interactions!(args, args.system_type, particle_vec, velocity_vec, neighbor_list, c_l, force_vec, energy, stress, mesh, vertex, cell_vertices, args.non_bonded_interactions)
        end
        
        energy3 = energy - energy1 - energy2
        # Slip-spring interactions, if active.
        if typeof(args.t_bonds) != No_T_Bonds
            energy, stress = apply_bonded_interactions!(args, particle_vec, force_vec, energy, stress, t_bond_vec, args.t_bonds)
            energy4 = energy - energy1 - energy2 - energy3
        end
        # energy1: bonded
        # energy2: anngles
        # energy3: non-bonded
        # energy4: slip-springs


        #println("Forces: ", force_vec)
        # Apply second integration.
        if typeof(args.thermostat) == Langevin_Thermostat
            update_velocities_langevin!(args, particle_vec, velocity_vec, force_vec, args.thermostat)
        else
            update_velocities!(args, particle_vec, velocity_vec, force_vec)
        end
        #update_velocities!(args, particle_vec, velocity_vec, force_vec)
        if args.velocity_rescale != 0 && step_i%args.velocity_rescale == 0.0 
            velocity_rescale!(particle_vec, velocity_vec)
        end
        #print(particle_vec)

        # Calcylate virial, kinetic energy and pressure.
        E_kin = calc_kinetic_energy(particle_vec, velocity_vec, E_kin)
        curr_p = calc_pressure(stress, E_kin, args.box)
        # Apply thermostat and barostat.
        apply_thermostat!(args, particle_vec, velocity_vec, neighbor_list, args.thermostat)   
        apply_barostat!(args, particle_vec, curr_p, args.barostat)

        # Output.
        if step_i%args.log_period == 0.0
            apply_logger(log_file, first_step, step_i, velocity_vec, energy, curr_p, particle_vec, args.logger)
        end
        if step_i%args.traj_period == 0.0
            apply_dump(traj_file, step_i, args, particle_vec, force_vec, velocity_vec, args.trajectorydump)
        end
        if args.shear_sample_period != 0 && step_i%args.shear_sample_period == 0.0 && typeof(args.shear) != No_Shear
            sample_velocity_profile!(args, particle_vec, velocity_vec, d_slab, z_bins, args.shear)
        end
        if step_i%args.n_steps == 0.0 && typeof(args.shear) != No_Shear
            apply_shear_log(shear_file, momentum_file, step_i, args, d_slab, z_bins, args.shear)
        end
        if step_i%args.n_steps == 0.0 && args.restart == true
            restart_file = open(args.restart_file, "w")
            apply_restart(restart_file, c_l, args, particle_vec, velocity_vec, bond_vec, angle_vec, args.system_type)
            close(restart_file)
        end
        if step_i%args.n_steps == 0.0 && typeof(args.t_bonds) != No_T_Bonds
            apply_t_bond_dump(entanglement_file, t_bond_vec)
        end
        f_sample_final = false
        if step_i%args.n_steps == 0.0 && f_sample_final == true
            write_final_forces(force_vec)
        end
        if step_i%f_sample_freq == 0.0 && f_sample_steps == true
            write_forces(force_file, force_vec)
        end
        
        # Calculate & output momentum drift

        P_tot = calc_totalmomentum(particle_vec, velocity_vec)
        if step_i > 1 && calc_drift == true
            momentum_drift!(momentum_drift, P_tot, P_tot_prev)
        end
        P_tot_prev = P_tot
        if step_i%args.n_steps == 0.0 && calc_drift == true
            write_drift(momentum_drift)
        end

        # Clear energies after simulation step ends. 
	    energy = 0.0e0
        #clear_energy!(energy_vec)
        clear_stress_tensor!(stress)
        clear_kin_energy_tensor!(E_kin) 
    end

################################################################################

end