function apply_restart(restart_file, l_chain, args::System, particle_vec::Vector{Particle}, velocity_vec, bond_vec::Vector{Bond}, angle_vec::Vector{Angle}, ::No_Mol)

    write(restart_file, "LAMMPS data file generated via init_lammps.jl\n \n")
    write(restart_file, "$(length(particle_vec)) atoms\n")
    n_spec = 1
    write(restart_file, "$n_spec atom types\n")
        
    write(restart_file, "$( - args.box[1] / 2) $(args.box[1] / 2) xlo xhi\n")
    write(restart_file, "$( - args.box[2] / 2) $(args.box[2] / 2) ylo yhi\n")
    write(restart_file, "$( - args.box[3] / 2) $(args.box[3] / 2) zlo zhi\n")
    write(restart_file, "\n")
    write(restart_file, "Masses\n")
    write(restart_file, "\n")
    for i in 1:n_spec
        write(restart_file, "$i $(particle_vec[i].mass)\n")
    end
    write(restart_file, "\n")
    write(restart_file, "Atoms\n")
    write(restart_file, "\n")
        
    i_counter = 0

    for i in 1:length(particle_vec)
        i_counter += 1
        write(restart_file, "$i_counter $(particle_vec[i_counter].i_mol) $(particle_vec[i_counter].spec) $(particle_vec[i_counter].charge) $(particle_vec[i_counter].pos[1]) $(particle_vec[i_counter].pos[2]) $(particle_vec[i_counter].pos[3]) $(particle_vec[i_counter].image[1]) $(particle_vec[i_counter].image[2]) $(particle_vec[i_counter].image[3])\n")
    end

    write(restart_file, "\n")
    write(restart_file, "Velocities\n")
    write(restart_file, "\n")
        
    i_counter = 0

    for i in 1:length(particle_vec)
        i_counter += 1
        write(restart_file, "$i_counter $(velocity_vec[i_counter][1]) $(velocity_vec[i_counter][2]) $(velocity_vec[i_counter][3])\n")
    end

end

function apply_restart(restart_file, l_chain, args::System, particle_vec::Vector{Particle}, velocity_vec, bond_vec::Vector{Bond}, angle_vec::Vector{Angle}, ::Molecular)
    
    n_mol = length(particle_vec) / l_chain

    write(restart_file, "LAMMPS data file generated via init_lammps.jl\n \n")
    write(restart_file, "$(length(particle_vec)) atoms\n")
    n_spec = 1
    write(restart_file, "$n_spec atom types\n")
    write(restart_file, "$(length(bond_vec)) bonds\n")
    n_spec = 1
    write(restart_file, "$n_spec bond types\n")
    write(restart_file, "$(length(angle_vec)) angles\n")
    n_spec = 1
    write(restart_file, "$n_spec angle types\n")
    write(restart_file, "\n")
        
    write(restart_file, "$( - args.box[1] / 2) $(args.box[1] / 2) xlo xhi\n")
    write(restart_file, "$( - args.box[2] / 2) $(args.box[2] / 2) ylo yhi\n")
    write(restart_file, "$( - args.box[3] / 2) $(args.box[3] / 2) zlo zhi\n")
    write(restart_file, "\n")
    write(restart_file, "Masses\n")
    write(restart_file, "\n")
    for i in 1:n_spec
        write(restart_file, "$i $(particle_vec[i].mass)\n")
    end
    write(restart_file, "\n")
    write(restart_file, "Atoms\n")
    write(restart_file, "\n")
        
    i_counter = 0

    for i in 1:n_mol
        for j in 1:l_chain
            i_counter += 1
            write(restart_file, "$i_counter $(particle_vec[i_counter].i_mol) $(particle_vec[i_counter].spec) $(particle_vec[i_counter].charge) $(particle_vec[i_counter].pos[1]) $(particle_vec[i_counter].pos[2]) $(particle_vec[i_counter].pos[3]) $(particle_vec[i_counter].image[1]) $(particle_vec[i_counter].image[2]) $(particle_vec[i_counter].image[3])\n")
        end
    end

    write(restart_file, "\n")
    write(restart_file, "Velocities\n")
    write(restart_file, "\n")
        
    i_counter = 0

    for i in 1:n_mol
        for j in 1:l_chain
            i_counter += 1
            write(restart_file, "$i_counter $(velocity_vec[i_counter][1]) $(velocity_vec[i_counter][2]) $(velocity_vec[i_counter][3])\n")
        end
    end
        
    write(restart_file, "\n")
    write(restart_file, "Bonds\n")
    write(restart_file, "\n")
        
    cl_counter = 1 # Chain length counter.
    global i_counter = 0 # Index counter.
    global i_counter2 = 0 # Bond/Angle counter.
    bond_type = 1

    for i in 1:n_mol
        global cl_counter = 1
        for j in 1:l_chain
            global i_counter += 1
            if cl_counter != l_chain
                global cl_counter += 1
                global i_counter2 += 1
                write(restart_file, "$i_counter2 $bond_type $i_counter $(i_counter + 1)\n")
            end
        end
    end

    write(restart_file, "\n")
    write(restart_file, "Angles\n")
    write(restart_file, "\n")
        
    global i_counter = 0 # Index counter.
    global i_counter2 = 0 # Bond/Angle counter.
    angle_type = 1

    for i in 1:n_mol
        global cl_counter = 1
        for j in 1:l_chain
            global i_counter += 1
            if cl_counter < l_chain - 1
                global cl_counter += 1
                global i_counter2 += 1
                write(restart_file, "$i_counter2 $angle_type $i_counter $(i_counter + 1) $(i_counter + 2)\n")
            end
        end
    end

end    