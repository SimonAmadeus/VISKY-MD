#include("system.jl")

# Read input for molecular system.
function read_lammps_x(filename::AbstractString, ::No_Mol)    

    xhi = 0 
    xlo = 0
    yhi = 0
    ylo = 0
    zhi = 0
    zlo = 0
    box = [0.0, 0.0, 0.0]
    c_l = 1
    number_particles = 0

    counter = 0
    counter_check = "check"

    open(filename) do file
        for line in eachline(file)
            
            a = split(line) # Create array from line.

            # Read header
            if length(a) > 1
                if a[2] == "atoms"
                    number_particles = parse(Int64, a[1])
                    println("$number_particles particles in the system.")
                    # Define arrays for particles, velocities, bonds, angles.
                    global particle_vec = Vector{Particle}(undef, number_particles)
                    global velocity_vec = [zeros(3) for _ in 1:number_particles]
                    global bond_vec = Vector{Bond}(undef, number_particles)
                    global angle_vec = Vector{Angle}(undef, number_particles)
                end
            end
            if length(a) > 2
                if a[2] == "atom" && a[3] == "types"
                    global ntypes_particles = a[1] 
                end
            end

            if length(a) > 3
                if a[3] == "xlo"
                    xlo = parse(Float64, a[1])
                    xhi = parse(Float64, a[2])
                end
            end
            if length(a) > 3
                if a[3] == "ylo"
                    ylo = parse(Float64, a[1])
                    yhi = parse(Float64, a[2])
                end
            end
            if length(a) > 3
                if a[3] == "zlo"
                    zlo = parse(Float64, a[1])
                    zhi = parse(Float64, a[2])
                    box = [(xhi - xlo), (yhi - ylo), (zhi - zlo)]
                end
            end

            # Read Masses
            if counter == 2 && counter_check == "Masses"
                global mass = parse(Float64, a[2])
                #global particle.mass = a[2]
                counter = 0
                counter_check = "check"
            end
            if counter == 1 && counter_check == "Masses"
                counter = 2
            end
            if isempty(a) == false && a[1] == "Masses"
                counter_check = "Masses"
                counter = 1
            end 

            # Read atom information -> index, mol. nr., type, charge, pos(x,y,z), image
            if counter >= 2 && counter_check == "Atoms"
                counter += 1
                index = parse(Int64, a[1])
                i_mol = parse(Int64, a[2])
                spec_i = parse(Int64, a[3])
                charge_i = parse(Float64, a[4])
                pos_i = [parse(Float64, a[5]), parse(Float64, a[6]), parse(Float64, a[7])]
                image_i = [parse(Int64, a[8]), parse(Int64, a[9]), parse(Int64, a[10])]
                global particle_i = Particle(i_mol, spec_i, mass, charge_i, pos_i, image_i)
                particle_vec[index] = particle_i
                if counter == (number_particles + 2)
                    counter = 0
                    counter_check = "check"
                end
            end
            if counter == 1 && counter_check == "Atoms"
                counter = 2
            end
            if isempty(a) == false && a[1] == "Atoms"
                counter_check = "Atoms"
                counter = 1
            end   
            
            if counter >= 2 && counter_check == "Velocities"
                counter += 1
                index = parse(Int64, a[1])
                vel_i = [parse(Float64, a[2]), parse(Float64, a[3]), parse(Float64, a[4])]
                velocity_vec[index] = vel_i
                if counter == (number_particles + 2)
                    counter = 0
                    counter_check = "check"
                end
            end
            if counter == 1 && counter_check == "Velocities"
                counter = 2
            end
            if isempty(a) == false && a[1] == "Velocities"
                counter_check = "Velocities"
                counter = 1
            end 

        end
    end
    return particle_vec, c_l, velocity_vec, bond_vec, angle_vec, box
end

# Read input for molecular system.
function read_lammps_x(filename::AbstractString, ::Molecular)    
    
    c_l = 0 # Chain length (1 for atoms).
    number_particles = 0

    xhi = 0 
    xlo = 0
    yhi = 0
    ylo = 0
    zhi = 0
    zlo = 0
    box = [0.0, 0.0, 0.0]

    counter = 0
    counter_check = "check"

    open(filename) do file
        for line in eachline(file)

            a = split(line) # Create array from line.

            # Read header.
            if length(a) > 1
                if a[2] == "atoms"
                    number_particles = parse(Int64, a[1])
                    println("$number_particles particles in the system.")
                    # Define arrays for particles, velocities, bonds, angles.
                    global particle_vec = Vector{Particle}(undef, number_particles)
                    global velocity_vec = [zeros(3) for _ in 1:number_particles]
                end
                if a[2] == "bonds"
                    global number_bonds = parse(Int64, a[1]) 
                    global bond_vec = Vector{Bond}(undef, number_bonds)
                end
                if a[2] == "angles"
                    global number_angles = parse(Int64, a[1])
                    global angle_vec = Vector{Angle}(undef, number_angles) 
                end
            end
            if length(a) > 2
                if a[2] == "atom" && a[3] == "types"
                    global ntypes_particles = a[1] 
                end
                if a[2] == "bond" && a[3] == "types"
                    global ntypes_bonds = a[1] 
                end
                if a[2] == "angle" && a[3] == "types"
                    global ntypes_angles = a[1] 
                end
            end
            
            # Box information.
            if length(a) > 3
                if a[3] == "xlo"
                    xlo = parse(Float64, a[1])
                    xhi = parse(Float64, a[2])
                end
            end
            if length(a) > 3
                if a[3] == "ylo"
                    ylo = parse(Float64, a[1])
                    yhi = parse(Float64, a[2])
                end
            end
            if length(a) > 3
                if a[3] == "zlo"
                    zlo = parse(Float64, a[1])
                    zhi = parse(Float64, a[2])
                    box = [(xhi - xlo), (yhi - ylo), (zhi - zlo)]
                end
            end


            # Read Masses.
            if counter == 2 && counter_check == "Masses"
                global mass = parse(Float64, a[2])
                #global particle.mass = a[2]
                counter = 0
                counter_check = "check"
            end
            if counter == 1 && counter_check == "Masses"
                counter = 2
            end
            if isempty(a) == false && a[1] == "Masses"
                counter_check = "Masses"
                counter = 1
            end 

            # Read atom information -> index, mol. nr., type, charge, pos(x,y,z), image.
            if counter >= 2 && counter_check == "Atoms"
                counter += 1
                index = parse(Int64, a[1])
                i_mol = parse(Int64, a[2])
                spec_i = parse(Int64, a[3])
                charge_i = parse(Float64, a[4])
                pos_i = [parse(Float64, a[5]), parse(Float64, a[6]), parse(Float64, a[7])]
                image_i = [parse(Int64, a[8]), parse(Int64, a[9]), parse(Int64, a[10])]
                global particle_i = Particle(i_mol, spec_i, mass, charge_i, pos_i, image_i)
                particle_vec[index] = particle_i
                if counter == (number_particles + 2)
                    c_l = Int(number_particles / i_mol)
                    counter = 0
                    counter_check = "check"
                end
            end
            if counter == 1 && counter_check == "Atoms"
                counter = 2
            end
            if isempty(a) == false && a[1] == "Atoms"
                counter_check = "Atoms"
                counter = 1
            end   

            # Read velocity information.
            if counter >= 2 && counter_check == "Velocities"
                counter += 1
                index = parse(Int64, a[1])
                vel_i = [parse(Float64, a[2]), parse(Float64, a[3]), parse(Float64, a[4])]
                velocity_vec[index] = vel_i
                if counter == (number_particles + 2)
                    counter = 0
                    counter_check = "check"
                end
            end
            if counter == 1 && counter_check == "Velocities"
                counter = 2
            end
            if isempty(a) == false && a[1] == "Velocities"
                counter_check = "Velocities"
                counter = 1
            end   

            # Read bond information
            if counter >= 2 && counter_check == "Bonds"
                counter += 1
                index = parse(Int64, a[1])
                type_i = parse(Int64, a[2])
                p_1_i = parse(Int64, a[3])
                p_2_i = parse(Int64, a[4])
                global bond_i = Bond(type_i, p_1_i, p_2_i)
                bond_vec[index] = bond_i
                if counter == (number_bonds + 2)
                    counter = 0
                    counter_check = "check"
                end
            end
            if counter == 1 && counter_check == "Bonds"
                counter = 2
            end
            if isempty(a) == false && a[1] == "Bonds"
                counter_check = "Bonds"
                counter = 1
            end   

            # Read angle information
            if counter >= 2 && counter_check == "Angles"
                counter += 1
                index = parse(Int64, a[1])
                type_i = parse(Int64, a[2])
                p_1_i = parse(Int64, a[3])
                p_2_i = parse(Int64, a[4])
                p_3_i = parse(Int64, a[5])
                global angle_i = Angle(type_i, p_1_i, p_2_i, p_3_i)
                angle_vec[index] = angle_i
                if counter == (number_angles + 2)
                    counter = 0
                    counter_check = "check"
                end
            end
            if counter == 1 && counter_check == "Angles"
                counter = 2
            end
            if isempty(a) == false && a[1] == "Angles"
                counter_check = "Angles"
                counter = 1
            end
        end
    end
    return particle_vec, c_l, velocity_vec, bond_vec, angle_vec, box
end