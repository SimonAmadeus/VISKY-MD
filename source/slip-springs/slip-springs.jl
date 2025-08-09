export T_Bonds_CTRL

struct T_Bonds_CTRL <: T_Bonds
    n_t_bonds::Int64 # Number of slip-springs.
    freq_mc::Int64 # Frequeny od Monte-Carlo migration.
    n_mc::Int64 # Number of migration steps.
    k::Float64 # Slip-spring force constant.
    r_0::Float64 # Eq. slip-spring bond length.
    restart::AbstractString # Start from restart file (typ. entanglements.in).
end

function select_random_chain_end(N, c_l)
    ind = Int(ceil(rand()*N))
    while ind != 1 && mod(ind, c_l) != 0 && mod(ind - 1, c_l) != 0
        ind = Int(ceil(rand()*N))
    end
    return ind 
end

function init_t_bonds!(args::System, c_l::Int64, particle_vec::Vector{Particle}, bond_vec::Vector{Bond}, ::No_T_Bonds) end

function read_t_bonds()

end

function init_t_bonds!(args::System, c_l::Int64, particle_vec::Vector{Particle}, bond_vec::Vector{Bond}, t_bonds::T_Bonds_CTRL) 
    
    # This function initializes slip-springs and connects two randon beads 
    # within a distance criterion. Beads which are connected by bonds, can not
    # be connected by slip-speings at the same time.

    max_capture = 3

    t_bond_vec = Vector{T_Bond}(undef, t_bonds.n_t_bonds)
    i_counter = 1
    check = true

    while i_counter <= length(t_bond_vec)    
        p_i = Int(ceil(rand() * length(particle_vec)))
        p_j = Int(ceil(rand() * length(particle_vec)))
        # Neighbbors already connected by bonds, will not be connected with
        # slip-springs.

        # Check for particle i being at chainend(end).
        if p_i != p_j && mod(p_i, c_l) == 0 && p_i - p_j != 1
            # Check, if particles are already connected by slip-springs.
            for i = 1:i_counter - 1
                if (p_i == t_bond_vec[i].p_1 && p_j == t_bond_vec[i].p_2) || (p_i == t_bond_vec[i].p_2 && p_j == t_bond_vec[i].p_1)
                    check = false
                end
            end
            if check == true
                ∆d = particle_vec[p_i].pos .- particle_vec[p_j].pos
                ∆d .-= args.box .* floor.((∆d .+ args.box / 2) ./ args.box) # PBCs
                r = norm(∆d)
                if r <= max_capture
                    t_bond_i = T_Bond(1, p_i, p_j)
                    t_bond_vec[i_counter] = t_bond_i
                    i_counter += 1
                end
            end
            check = true
        # Check for particle j being at chain end (end).
        elseif p_i != p_j && mod(p_j, c_l) == 0 && p_j - p_i != 1
            for i = 1:i_counter - 1
                if (p_i == t_bond_vec[i].p_1 && p_j == t_bond_vec[i].p_2) || (p_i == t_bond_vec[i].p_2 && p_j == t_bond_vec[i].p_1)
                    check = false
                end
            end
            if check == true
                ∆d = particle_vec[p_i].pos .- particle_vec[p_j].pos
                ∆d .-= args.box .* floor.((∆d .+ args.box / 2) ./ args.box) # PBCs
                r = norm(∆d)
                if r <= max_capture
                    t_bond_i = T_Bond(1, p_i, p_j)
                    t_bond_vec[i_counter] = t_bond_i
                    i_counter += 1
                end
            end
            check = true
        # Check for particle i being at chain end (start).
        elseif p_i != p_j && (p_i == 1 || mod(p_i - 1, c_l) == 0) && p_j - p_i != 1
            for i = 1:i_counter - 1
                if (p_i == t_bond_vec[i].p_1 && p_j == t_bond_vec[i].p_2) || (p_i == t_bond_vec[i].p_2 && p_j == t_bond_vec[i].p_1)
                    check = false
                end
            end
            if check == true
                ∆d = particle_vec[p_i].pos .- particle_vec[p_j].pos
                ∆d .-= args.box .* floor.((∆d .+ args.box / 2) ./ args.box) # PBCs
                r = norm(∆d)
                if r <= max_capture
                    t_bond_i = T_Bond(1, p_i, p_j)
                    t_bond_vec[i_counter] = t_bond_i
                    i_counter += 1
                end
            end
            check = true
        # Check for particle j being at chain end (start).
        elseif p_i != p_j && (p_j == 1 || mod(p_j - 1, c_l) == 0) && p_i - p_j != 1
            for i = 1:i_counter - 1
                if (p_i == t_bond_vec[i].p_1 && p_j == t_bond_vec[i].p_2) || (p_i == t_bond_vec[i].p_2 && p_j == t_bond_vec[i].p_1)
                    check = false
                end
            end
            if check == true
                ∆d = particle_vec[p_i].pos .- particle_vec[p_j].pos
                ∆d .-= args.box .* floor.((∆d .+ args.box / 2) ./ args.box) # PBCs
                r = norm(∆d)
                if r <= max_capture
                    t_bond_i = T_Bond(1, p_i, p_j)
                    t_bond_vec[i_counter] = t_bond_i
                    i_counter += 1
                end
            end
            check = true
        elseif p_i != p_j && abs(p_i - p_j) != 1
            for i = 1:i_counter - 1
                if (p_i == t_bond_vec[i].p_1 && p_j == t_bond_vec[i].p_2) || (p_i == t_bond_vec[i].p_2 && p_j == t_bond_vec[i].p_1)
                    check = false
                end
            end
            if check == true
                ∆d = particle_vec[p_i].pos .- particle_vec[p_j].pos
                ∆d .-= args.box .* floor.((∆d .+ args.box / 2) ./ args.box) # PBCs
                r = norm(∆d)
                if r <= max_capture
                    t_bond_i = T_Bond(1, p_i, p_j)
                    t_bond_vec[i_counter] = t_bond_i
                    i_counter += 1
                end
            end
            check = true
        else
            check = true
        end
    end
    return t_bond_vec
end

function migration!(first_step, args::System, c_l::Int64, particle_vec::Vector{Particle}, t_bond_vec::Vector{T_Bond}, t_bonds::T_Bonds_CTRL) 
    
    # Slip-spring migration. If a slip-spring is anchored at a chain end, it
    # might move to a randomly chosen chain end under the condition that the
    # metropolis criterion is fullfilled. 
    
    ∆i = 0
    ∆j = 0
    
    if first_step == false
        for i in 1:length(t_bond_vec)
            p_i = t_bond_vec[i].p_1
            p_j = t_bond_vec[i].p_2

            # Check anchor point 1 for chain end.
            if p_i == 1 || mod(p_i, c_l) == 0 || mod(p_i - 1, c_l) == 0
                ind = select_random_chain_end(length(particle_vec), c_l)
                while ind == p_i || ind == p_j
                    ind = select_random_chain_end(length(particle_vec), c_l)
                end
                energy_before = bond_energy(args, p_i, p_j, particle_vec)
                energy_after = bond_energy(args, ind, p_j, particle_vec) 
                if energy_after < energy_before
                    t_bond_vec[i].p_1 = ind
                else
                    random_var = rand()
                    if exp(energy_before - energy_after) > random_var
                        t_bond_vec[i].p_1 = ind
                    end
                end    
            end

            # Anchor point 2.
            if p_j == 1 || mod(p_j, c_l) == 0 || mod(p_j - 1, c_l) == 0
                ind = select_random_chain_end(length(particle_vec), c_l)
                while ind == p_i || ind == p_j
                    ind = select_random_chain_end(length(particle_vec), c_l)
                end
                energy_before = bond_energy(args, p_i, p_j, particle_vec)
                energy_after = bond_energy(args, p_i, ind, particle_vec) 
                if energy_after < energy_before
                    t_bond_vec[i].p_2 = ind
                else
                    random_var = rand()
                    if exp(energy_before - energy_after) > random_var
                        t_bond_vec[i].p_2 = ind
                    end
                end    
            end
        
        end
    end

    for i_mc in 1:t_bonds.n_mc # Number of Monte Carlo steps, set in the control file.
        for i in 1:length(t_bond_vec)
            p_i = t_bond_vec[i].p_1
            p_j = t_bond_vec[i].p_2

            # Determine direction of entanglementt movement.
            # If slip-spring is at chain end ∆i is + 1/- 1, otherwise random.

            if p_i ==  1 || mod(p_i - 1, c_l) == 0 
                ∆i = 1
            elseif mod(p_i, c_l) == 0
                ∆i = - 1
            else
                random_var = rand()
                if random_var > 0.5
                    ∆i = 1
                else 
                    ∆i = -1
                end
            end
            
            if p_j ==  1 || mod(p_j - 1, c_l) == 0
                ∆j = 1
            elseif mod(p_j, c_l) == 0
                ∆j = - 1
            else 
                random_var = rand()
                if random_var > 0.5
                    ∆j = 1
                else 
                    ∆j = -1
                end
            end

            if p_i + ∆i == p_j + ∆j
                # Remove and reate new slip-spring.
                p_i = Int(ceil(rand() * length(particle_vec)))
                p_j = Int(ceil(rand() * length(particle_vec)))
                if p_i != p_j
                    t_bond_vec[i].p_1 = p_i
                    t_bond_vec[i].p_2 = p_j
                end 
            else
                # Compare bond energies before and after slip-spring migrgation.
                energy_before = bond_energy(args, p_i, p_j, particle_vec)
                energy_after = bond_energy(args, p_i + ∆i, p_j + ∆j, particle_vec) 

                if energy_after < energy_before
                    t_bond_vec[i].p_1 = p_i + ∆i
                    t_bond_vec[i].p_2 = p_j + ∆j
                else
                    random_var = rand()
                    if exp(energy_before - energy_after) > random_var
                        t_bond_vec[i].p_1 = p_i + ∆i
                        t_bond_vec[i].p_2 = p_j + ∆j    
                    end
                end
            end
        end
    end
end

function apply_bonded_interactions!(args::System, particle_vec::Vector{Particle}, force_vec, energy::Float64, stress::Vector{Float64}, t_bond_vec::Vector{T_Bond}, t_bonds::T_Bonds_CTRL) 
    
    # 0.5×k×(r-r0)^2
    for i = 1:length(t_bond_vec)
        pos_1 = particle_vec[t_bond_vec[i].p_1].pos .+ particle_vec[t_bond_vec[i].p_1].image .* args.box
        pos_2 = particle_vec[t_bond_vec[i].p_2].pos .+ particle_vec[t_bond_vec[i].p_2].image .* args.box
        ∆d = pos_2 .- pos_1
        r = norm(∆d)
        f = - t_bonds.k * (r - t_bonds.r_0) 
        force_divr = f / r
        force_vec[t_bond_vec[i].p_1] .-= ∆d .* force_divr 
        force_vec[t_bond_vec[i].p_2] .+= ∆d .* force_divr
        energy += 0.5 * t_bonds.k * (r - t_bonds.r_0)^2
    end
    return energy, stress
end

function apply_t_bond_dump(out_file, t_bond_vec::Vector{T_Bond})
    string_out = ""
    string_out *= string(length(t_bond_vec))*" Slip-springs"
    string_out *= "\n"
    write(out_file, string_out)

    for i in 1:length(t_bond_vec)
        string_out = ""
        string_out = string(t_bond_vec[i].type)*" "*string(t_bond_vec[i].p_1)*" "*string(t_bond_vec[i].p_2)
        string_out *= "\n"
        write(out_file, string_out)
    end
end