function create_cell_list(args::System, particle_vec::Vector{Particle}, cell_size::Float64)
    #cell_size = r_cutoff
    n_cells = ceil.(Int, args.box ./ cell_size)
    n_total_cells = prod(n_cells)

    cell_list = Vector{Vector{Int}}(undef, n_total_cells)
    for i in 1:n_total_cells
        cell_list[i] = []
    end

    particle_to_cell = Vector{Int}(undef, length(particle_vec))

    for (i, particle) in enumerate(particle_vec)
        cell_indices = floor.(Int, (particle.pos .+ args.box / 2) ./ cell_size) .+ 1
        cell_index = sub2ind(n_cells, cell_indices...)
        push!(cell_list[cell_index], i)
        particle_to_cell[i] = cell_index
    end

    return cell_list, particle_to_cell
end

function create_neighbor_list(args::System, particle_vec::Vector{Particle}, r_cutoff::Float64)
    neighbor_list = Vector{Vector{Int}}(undef, length(particle_vec))
    for i in 1:length(particle_vec)
        neighbor_list[i] = []
    end

    cell_list, particle_to_cell = create_cell_list(args, particle_vec, r_cutoff)

    for i = 1:length(particle_vec)
        current_cell = particle_to_cell[i]

        for j in (i + 1):length(particle_vec)
            if current_cell == particle_to_cell[j]
                # Calculate distance only if particles are in the same or neighboring cells
                ∆d = particle_vec[i].pos .- particle_vec[j].pos
                ∆d .-= args.box .* floor.((∆d .+ args.box / 2) ./ args.box) # PBCs
                r = norm(∆d)

                if r < r_cutoff
                    push!(neighbor_list[i], j)
                    push!(neighbor_list[j], i)
                end
            end
        end
    end

    return neighbor_list
end