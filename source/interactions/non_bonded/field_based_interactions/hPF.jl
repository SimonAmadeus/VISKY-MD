export original_hPF_Interactions, spectral_hPF_Interactions

struct original_hPF_Interactions <: Non_Bonded_Interactions 
    κ::Float64
    mesh_points::Vector{Int64}
end

function apply_non_bonded_interactions!(args::System, ::System_Type, particle_vec::Vector{Particle}, velocity_vec, neighbor_list::Any, c_l::Int64, force_vec, energy::Float64, stress::Vector{Float64}, mesh::Mesh, vertex::Array{Int64, 3}, cell_vertices::Array{Int64, 2}, HPF::Non_Bonded_Interactions)
    σ = 1 # Parameter for Gaussian cloud in cell (particle diameter).
    clear_mesh!(mesh)
    d_pos = zeros(Float64, mesh.NC_total) 
    for i = 1:length(particle_vec)
        cloudincell!(particle_vec[i].pos, mesh, cell_vertices, i)
        #gaussian_distribution!(particle_vec[i].pos, mesh, i, σ)
    end

    d_pos = positional_density(particle_vec, mesh, cell_vertices)
    
    if !mesh.use_staggered
        energy = Particle_Field_Interactions!(particle_vec, force_vec, energy, mesh, vertex, cell_vertices, d_pos, HPF)
    else
        energy = Particle_Field_Interactions_Staggered_Lattice!(particle_vec, force_vec, energy, mesh, vertex, cell_vertices, d_pos, HPF)
    end

    return energy, stress
end

function Particle_Field_Interactions!(particle_vec::Vector{Particle}, force_vec, energy::Float64, mesh::Mesh, vertex::Array{Int64, 3}, cell_vertices::Array{Int64, 2}, d_pos::Vector{Float64}, HPF::original_hPF_Interactions)
    # Calculate forces and energies through interpolation of the gradients
    # average density
    avg_den = length(particle_vec) / prod(mesh.boxsize)
    # constants for force calculation
    ikomp = 1 / HPF.κ
    # inverse density per cell
    idpc = 1 / (avg_den * prod(mesh.edge_size))

    pic = zeros(Float64, 8)

    # force on the grids: -1 * derivative of the potential on grids
    grid_grad_x, grid_grad_y, grid_grad_z = Grad_DensVertex(vertex, mesh)

    for i = 1:length(particle_vec)

        position_x = particle_vec[i].pos[1] + 0.5 * mesh.boxsize[1]
        position_y = particle_vec[i].pos[2] + 0.5 * mesh.boxsize[2]
        position_z = particle_vec[i].pos[3] + 0.5 * mesh.boxsize[3]
        
        ix::Int64 = floor(position_x / mesh.edge_size[1]) + 1
        iy::Int64 = floor(position_y / mesh.edge_size[2]) + 1
        iz::Int64 = floor(position_z / mesh.edge_size[3]) + 1
        
        cell_index = getCellIndex(ix, iy, iz, mesh.NC)

        v1 = cell_vertices[cell_index, 1]
        v2 = cell_vertices[cell_index, 2]
        v3 = cell_vertices[cell_index, 3]
        v4 = cell_vertices[cell_index, 4]
        v5 = cell_vertices[cell_index, 5]
        v6 = cell_vertices[cell_index, 6]
        v7 = cell_vertices[cell_index, 7]
        v8 = cell_vertices[cell_index, 8]
        
        δx = position_x - floor(position_x / mesh.edge_size[1]) * mesh.edge_size[1]
        δy = position_y - floor(position_y / mesh.edge_size[2]) * mesh.edge_size[2]
        δz = position_z - floor(position_z / mesh.edge_size[3]) * mesh.edge_size[3]

        pic[1] = (mesh.edge_size[1] - δx) * (mesh.edge_size[2] - δy) * (mesh.edge_size[3] - δz) / prod(mesh.edge_size)
        pic[2] = (δx) * (mesh.edge_size[2] - δy) * (mesh.edge_size[3] - δz) / prod(mesh.edge_size)
        pic[3] = (mesh.edge_size[1] - δx) * (δy) * (mesh.edge_size[3] - δz) / prod(mesh.edge_size)
        pic[4] = (δx) * (δy) * (mesh.edge_size[3] - δz) / prod(mesh.edge_size)
        pic[5] = (mesh.edge_size[1] - δx) * (mesh.edge_size[2] - δy) * (δz) / prod(mesh.edge_size)
        pic[6] = (δx) * (mesh.edge_size[2] - δy) * (δz) / prod(mesh.edge_size)
        pic[7] = (mesh.edge_size[1] - δx) * (δy) * (δz) / prod(mesh.edge_size)
        pic[8] = (δx) * (δy) * (δz) / prod(mesh.edge_size)

        # grid grad calculation of the forces
        force_vec[i][1] += - grid_grad_x[v1] * pic[1] * ikomp * idpc
        force_vec[i][1] += - grid_grad_x[v2] * pic[2] * ikomp * idpc
        force_vec[i][1] += - grid_grad_x[v3] * pic[3] * ikomp * idpc
        force_vec[i][1] += - grid_grad_x[v4] * pic[4] * ikomp * idpc
        force_vec[i][1] += - grid_grad_x[v5] * pic[5] * ikomp * idpc
        force_vec[i][1] += - grid_grad_x[v6] * pic[6] * ikomp * idpc  
        force_vec[i][1] += - grid_grad_x[v7] * pic[7] * ikomp * idpc
        force_vec[i][1] += - grid_grad_x[v8] * pic[8] * ikomp * idpc

        force_vec[i][2] += - grid_grad_y[v1] * pic[1] * ikomp * idpc
        force_vec[i][2] += - grid_grad_y[v2] * pic[2] * ikomp * idpc
        force_vec[i][2] += - grid_grad_y[v3] * pic[3] * ikomp * idpc
        force_vec[i][2] += - grid_grad_y[v4] * pic[4] * ikomp * idpc
        force_vec[i][2] += - grid_grad_y[v5] * pic[5] * ikomp * idpc
        force_vec[i][2] += - grid_grad_y[v6] * pic[6] * ikomp * idpc
        force_vec[i][2] += - grid_grad_y[v7] * pic[7] * ikomp * idpc
        force_vec[i][2] += - grid_grad_y[v8] * pic[8] * ikomp * idpc

        force_vec[i][3] += - grid_grad_z[v1] * pic[1] * ikomp * idpc
        force_vec[i][3] += - grid_grad_z[v2] * pic[2] * ikomp * idpc
        force_vec[i][3] += - grid_grad_z[v3] * pic[3] * ikomp * idpc
        force_vec[i][3] += - grid_grad_z[v4] * pic[4] * ikomp * idpc
        force_vec[i][3] += - grid_grad_z[v5] * pic[5] * ikomp * idpc
        force_vec[i][3] += - grid_grad_z[v6] * pic[6] * ikomp * idpc
        force_vec[i][3] += - grid_grad_z[v7] * pic[7] * ikomp * idpc
        force_vec[i][3] += - grid_grad_z[v8] * pic[8] * ikomp * idpc

        # gaussian gradients
#=        
        force_vec[i].force[1] += - mesh.gradient_grid[v1, 1] * pic[1] * ikomp * idpc
        force_vec[i].force[1] += - mesh.gradient_grid[v2, 1] * pic[2] * ikomp * idpc
        force_vec[i].force[1] += - mesh.gradient_grid[v3, 1] * pic[3] * ikomp * idpc
        force_vec[i].force[1] += - mesh.gradient_grid[v4, 1] * pic[4] * ikomp * idpc
        force_vec[i].force[1] += - mesh.gradient_grid[v5, 1] * pic[5] * ikomp * idpc
        force_vec[i].force[1] += - mesh.gradient_grid[v6, 1] * pic[6] * ikomp * idpc
        force_vec[i].force[1] += - mesh.gradient_grid[v7, 1] * pic[7] * ikomp * idpc
        force_vec[i].force[1] += - mesh.gradient_grid[v8, 1] * pic[8] * ikomp * idpc

        force_vec[i].force[2] += - mesh.gradient_grid[v1, 2] * pic[1] * ikomp * idpc
        force_vec[i].force[2] += - mesh.gradient_grid[v2, 2] * pic[2] * ikomp * idpc
        force_vec[i].force[2] += - mesh.gradient_grid[v3, 2] * pic[3] * ikomp * idpc
        force_vec[i].force[2] += - mesh.gradient_grid[v4, 2] * pic[4] * ikomp * idpc
        force_vec[i].force[2] += - mesh.gradient_grid[v5, 2] * pic[5] * ikomp * idpc
        force_vec[i].force[2] += - mesh.gradient_grid[v6, 2] * pic[6] * ikomp * idpc
        force_vec[i].force[2] += - mesh.gradient_grid[v7, 2] * pic[7] * ikomp * idpc
        force_vec[i].force[2] += - mesh.gradient_grid[v8, 2] * pic[8] * ikomp * idpc

        force_vec[i].force[3] += - mesh.gradient_grid[v1, 3] * pic[1] * ikomp * idpc
        force_vec[i].force[3] += - mesh.gradient_grid[v2, 3] * pic[2] * ikomp * idpc
        force_vec[i].force[3] += - mesh.gradient_grid[v3, 3] * pic[3] * ikomp * idpc
        force_vec[i].force[3] += - mesh.gradient_grid[v4, 3] * pic[4] * ikomp * idpc
        force_vec[i].force[3] += - mesh.gradient_grid[v5, 3] * pic[5] * ikomp * idpc
        force_vec[i].force[3] += - mesh.gradient_grid[v6, 3] * pic[6] * ikomp * idpc
        force_vec[i].force[3] += - mesh.gradient_grid[v7, 3] * pic[7] * ikomp * idpc
        force_vec[i].force[3] += - mesh.gradient_grid[v8, 3] * pic[8] * ikomp * idpc
=#
        energy += d_pos[i] * ikomp * idpc

    end
    return energy
end

function Particle_Field_Interactions_Staggered_Lattice!(particle_vec::Vector{Particle}, force_vec, energy::Float64, mesh::Mesh, vertex::Array{Int64, 3}, cell_vertices::Array{Int64, 2}, d_pos::Vector{Float64}, HPF::original_hPF_Interactions)
    # Calculate forces and energies through interpolation of the gradients, 
    # which are defined on a staggered lattice
    avg_den = length(particle_vec) / prod(mesh.boxsize)
    # constants for force calculation
    ikomp = 1 / HPF.κ
    # dpc = densitry per cell (i = inverse)
    # idpc = constant1 in occam
    idpc = 1 / (avg_den * prod(mesh.edge_size))

    pic = zeros(Float64, 8)
    gdx = zeros(Float64, length(particle_vec))        
    gdy = zeros(Float64, length(particle_vec))
    gdz = zeros(Float64, length(particle_vec))

    # Half the edge length of a cell and the inverse
    es_x = mesh.edge_size[1]
    es_y = mesh.edge_size[2]
    es_z = mesh.edge_size[3]
    es_ix = 1 / es_x
    es_iy = 1 / es_y
    es_iz = 1 / es_z
    es_x2 = es_x * 0.5
    es_y2 = es_y * 0.5
    es_z2 = es_z * 0.5
    es_ix2 = 1 / es_x2
    es_iy2 = 1 / es_y2
    es_iz2 = 1 / es_z2

    # force on the grids: -1 * derivative of the potential on grids
    Grad_staggered_lattice!(vertex, mesh)

    for i = 1:length(particle_vec)

        position_x = particle_vec[i].pos[1] + 0.5 * mesh.boxsize[1]
        position_y = particle_vec[i].pos[2] + 0.5 * mesh.boxsize[2]
        position_z = particle_vec[i].pos[3] + 0.5 * mesh.boxsize[3]
        
        ix::Int64 = floor(position_x / mesh.edge_size[1]) + 1
        iy::Int64 = floor(position_y / mesh.edge_size[2]) + 1
        iz::Int64 = floor(position_z / mesh.edge_size[3]) + 1
        
        cell_index = getCellIndex(ix, iy, iz, mesh.NC)

        v1 = cell_vertices[cell_index, 1]
        v2 = cell_vertices[cell_index, 2]
        v3 = cell_vertices[cell_index, 3]
        v4 = cell_vertices[cell_index, 4]
        v5 = cell_vertices[cell_index, 5]
        v6 = cell_vertices[cell_index, 6]
        v7 = cell_vertices[cell_index, 7]
        v8 = cell_vertices[cell_index, 8]
        
        δx = position_x - floor(position_x / mesh.edge_size[1]) * mesh.edge_size[1]
        δy = position_y - floor(position_y / mesh.edge_size[2]) * mesh.edge_size[2]
        δz = position_z - floor(position_z / mesh.edge_size[3]) * mesh.edge_size[3]
        
        # determination of closest vertex in cell for gradient interpolation
        igx::Int64 = floor(δx * es_ix2)
        igy::Int64 = floor(δy * es_iy2)
        igz::Int64 = floor(δz * es_iz2)

        ig_scell = 1 + igx + 2 * igy + 4 * igz
        ig_cell = cell_vertices[cell_index, ig_scell]

        shift_x = es_x2 - igx * es_x
        shift_y = es_y2 - igy * es_y
        shift_z = es_z2 - igz * es_z

        xt1 = (δx + shift_x) * es_ix
        yt1 = (δy + shift_y) * es_iy
        zt1 = (δz + shift_z) * es_iz

        xt2 = 1 - xt1
        yt2 = 1 - yt1
        zt2 = 1 - zt1

        # Tri-linear weighting factors for density interpolation in the cell
        volinv = 1.0 / prod(mesh.edge_size)
        pic[1] = (mesh.edge_size[1] - δx)*(mesh.edge_size[2] - δy)*(mesh.edge_size[3] - δz)*volinv
        pic[2] = δx*(mesh.edge_size[2] - δy)*(mesh.edge_size[3] - δz)*volinv
        pic[3] = (mesh.edge_size[1] - δx)*δy*(mesh.edge_size[3] - δz)*volinv
        pic[4] = δx*δy*(mesh.edge_size[3] - δz)*volinv
        pic[5] = (mesh.edge_size[1] - δx)*(mesh.edge_size[2] - δy)*δz*volinv
        pic[6] = δx*(mesh.edge_size[2] - δy)*δz*volinv
        pic[7] = (mesh.edge_size[1] - δx)*δy*δz*volinv
        pic[8] = δx*δy*δz*volinv

        gdx = xt1 * mesh.gradient_grid[ig_cell, 1] + xt2 * mesh.gradient_grid[ig_cell, 2]
        gdy = yt1 * mesh.gradient_grid[ig_cell, 3] + yt2 * mesh.gradient_grid[ig_cell, 4]
        gdz = zt1 * mesh.gradient_grid[ig_cell, 5] + zt2 * mesh.gradient_grid[ig_cell, 6]

        # Energy and force calculation
        
        # staggered lattice force calculation
        force_vec[i][1] += - gdx * ikomp * idpc
        force_vec[i][2] += - gdy * ikomp * idpc
        force_vec[i][3] += - gdz * ikomp * idpc
  
        energy += d_pos[i] * ikomp * idpc

    end

    return energy

end