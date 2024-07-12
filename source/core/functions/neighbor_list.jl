export standard_NL, CL_NL

struct standard_NL <: NL end
struct CL_NL <: NL end # Cell lsit based neighborlist.

function create_neighbor_list!(args::System, particle_vec::Vector{Particle}, neighbor_list::Vector{Vector{Int}}, n_cells::Vector{Int64}, r_cutoff::Float64, ::No_NL) end

function create_cell_list(args::System, particle_vec::Vector{Particle}, cell_size::Vector{Float64}, n_cells::Vector{Int64})
    #cell_size = r_cutoff
    n_total_cells = prod(n_cells)
    cell_list = Vector{Vector{Int}}(undef, n_total_cells)
    for i in 1:n_total_cells
        cell_list[i] = []
    end
    particle_to_cell = Vector{Int}(undef, length(particle_vec))
    for (i, particle) in enumerate(particle_vec)
        cell_indices = min.(floor.(Int, (particle.pos .+ args.box / 2) ./ cell_size) .+ 1, n_cells)
        cell_index = (cell_indices[3] - 1) * n_cells[1] * n_cells[2] + (cell_indices[2] - 1) * n_cells[1] + cell_indices[1]
        if cell_index < 0
            println(particle_vec[i])
        end
        push!(cell_list[cell_index], i)
        particle_to_cell[i] = cell_index
    end
    return cell_list, particle_to_cell
end

function create_neighbor_list!(args::System, particle_vec::Vector{Particle}, neighbor_list::Vector{Vector{Int}}, n_cells::Vector{Int64}, r_cutoff::Float64, NL::CL_NL)
    cell_size = args.box ./ n_cells
    cell_list, particle_to_cell = create_cell_list(args, particle_vec, cell_size, n_cells)
    
    for i in 1:length(neighbor_list)
        empty!(neighbor_list[i])
    end
    
    for (i, particle) in enumerate(particle_vec)
        current_cell_index = particle_to_cell[i]
        current_cell_indices = [(current_cell_index - 1) % n_cells[1] + 1, 
                                ((current_cell_index - 1) ÷ n_cells[1]) % n_cells[2] + 1, 
                                (current_cell_index - 1) ÷ (n_cells[1] * n_cells[2]) + 1]
        for dx in -1:1, dy in -1:1, dz in -1:1
            neighbor_cell_indices = mod.(current_cell_indices .+ [dx, dy, dz] .- 1, n_cells) .+ 1
            neighbor_cell_index = (neighbor_cell_indices[3] - 1) * n_cells[1] * n_cells[2] + 
                                  (neighbor_cell_indices[2] - 1) * n_cells[1] + neighbor_cell_indices[1]
            for j in cell_list[neighbor_cell_index]
                if j != i
                    ∆r = particle_vec[i].pos .- particle_vec[j].pos
                    ∆r .-= args.box .* floor.((∆r .+ args.box / 2) ./ args.box)
                    r = norm(∆r)
                    if r < r_cutoff
                        push!(neighbor_list[i], j)
                    end
                end
            end
        end
    end
end

function create_neighbor_list!(args::System, particle_vec::Vector{Particle}, neighbor_list::Vector{Vector{Int}}, n_cells::Vector{Int64}, r_cutoff::Float64, NL::standard_NL)
    # First clear neighbor list, then recalculate.
    #neighbor_list = empty!(neighbor_list)
    for i in 1:length(neighbor_list)
        empty!(neighbor_list[i])
    end

    for i = 2:length(particle_vec)
        for j = 1:i - 1
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
