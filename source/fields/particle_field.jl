# particle field subroutines:
# - initialize mesh, clear mesh, various pbc functions
# - CIC, gaussian density
# - gradient interpolation 

# Example of a system with 3x3 cells:
# we have 64 vertices (4x4), but effectively only 27 (3x3)

#   1-----2-----3-----1         10----11----12----10     19----20----21----19     1-----2-----3-----1
#   |     |     |     |         |     |     |     |      |     |     |     |      |     |     |     |
#   |  7  |  8  |  9  |         | 16  | 17  | 18  |      | 25  | 26  | 27  |      | 25  | 26  | 27  |
#   |     |     |     |         |     |     |     |      |     |     |     |      |     |     |     |
#   7-----8-----9-----7        16----17----18-----16     25----26----27----25     7-----8-----9-----7
#   |     |     |     |  upper  |     |     |     |      |     |     |     |      |     |     |     | 
#   |  4  |  5  |  6  |     --> | 13  | 14  | 15  | -->  | 22  | 23  | 24  | -->  | 22  | 23  | 24  |
#   |     |     |     |  floor  |     |     |     |      |     |     |     |      |     |     |     |
#   4-----5-----6-----4        13----14----15-----13    22----23----24-----22     4-----5-----6-----4
#   |     |     |     |         |     |     |     |      |     |     |     |      |     |     |     | 
#   |  1  |  2  |  3  |         | 10  | 11  | 12  |      | 19  | 20  | 21  |      | 19  | 20  | 21  |
#   |     |     |     |         |     |     |     |      |     |     |     |      |     |     |     |
#   1-----2-----3-----1        10----11----12-----10    19----20----21-----19     1-----2-----3-----1

function init_mesh(boxsize::Vector{Float64}, NC::Vector{Int64})
    # edge is the length of a cell
    edge_x = boxsize[1] / NC[1]
    edge_y = boxsize[2] / NC[2]
    edge_z = boxsize[3] / NC[3]
    edge_size = [edge_x, edge_y, edge_z]
    # NV is the number of vertices in each dimension
    # NC is the number of cells
    NC_total::Int64 = NC[1] * NC[2] * NC[3]
    # number of vertices (practically, due to PBCs: NV_total = NC_total)
    NV_total::Int64 = (NC[1] + 1) * (NC[2] + 1) * (NC[3] + 1)
    # the particle grid maps the densities, the gradient grid the gradients on 
    # a staggered lattice
    particle_grid = zeros(Float64, NC_total)
    gradient_grid = zeros(Float64, NC_total, 6)
    use_staggered = false
    return Mesh(edge_size, NC, NC_total, NV_total, particle_grid, gradient_grid, boxsize, use_staggered)
end

function init_vertex_matrix(mesh::Mesh)
    # initialize vertex array -> gives you a vertex number for given indices 

    vertex = zeros(Int, mesh.NC[1] + 1, mesh.NC[2] + 1, mesh.NC[3] + 1)
    index_v = 1
    for k in 1:mesh.NC[3]
        for j in 1:mesh.NC[2]
            for i in 1:mesh.NC[1]
                vertex[i, j, k] = index_v
                index_v += 1
            end
        end
    end

    x_max = mesh.NC[1] + 1
    y_max = mesh.NC[2] + 1
    z_max = mesh.NC[3] + 1
    
    # corners
    # all corner points refer to the same point (#1)
    for k in [1, z_max]
        for j in [1, y_max]
            for i in [1, x_max]
                vertex[i, j, k] = 1
            end
        end
    end
    # edges
    # each point on a edge exists four times
    # we need three loops, to go along the four edges of each dimension 

    for k in [1, z_max]
        for j in [1, y_max]
            for i in range(2, mesh.NC[1])
                vertex[i, j, k] = vertex[i, 1, 1]
            end
        end
    end
    for k in [1, z_max]
        for j in range(2, mesh.NC[2])
            for i in [1, x_max]
                vertex[i, j, k] = vertex[1, j, 1]
            end
        end
    end
    for k in range(2, mesh.NC[3])
        for j in [1, y_max]
            for i in [1, x_max]
                vertex[i, j, k] = vertex[1, 1, k]
            end
        end
    end
    
    # surfaces
    # each point on a surface exists only two times
    for k in range(2, mesh.NC[3])
        for j in range(2, mesh.NC[2])
            vertex[x_max, j, k] = vertex[1, j, k]
        end
    end
    for k in range(2, mesh.NC[3])
        for i in range(2, mesh.NC[1])
            vertex[i, y_max, k] = vertex[i, 1, k]
        end
    end
    for j in range(2, mesh.NC[2])
        for i in range(2, mesh.NC[1])
            vertex[i, j, z_max] = vertex[i, j, 1]
        end
    end
    
    return vertex
end

function init_cell_vertices(mesh::Mesh, vertex::Array{Int64, 3})
# initialize cell vertices -> gives indices for all eight vertices in a given cell

    cell_vertices = zeros(Int, mesh.NC_total, 8)
    for k in 1:mesh.NC[3]
        for j in 1:mesh.NC[2]
            for i in 1:mesh.NC[1]
                current_index = vertex[i, j, k]
                cell_vertices[current_index, 1] = vertex[i, j, k]
                cell_vertices[current_index, 2] = vertex[i + 1, j, k]
                cell_vertices[current_index, 3] = vertex[i, j + 1, k]
                cell_vertices[current_index, 4] = vertex[i + 1, j + 1, k ]
                cell_vertices[current_index, 5] = vertex[i, j, k + 1]
                cell_vertices[current_index, 6] = vertex[i + 1, j, k + 1]
                cell_vertices[current_index, 7] = vertex[i, j + 1, k + 1]
                cell_vertices[current_index, 8] = vertex[i + 1, j + 1, k + 1]
            end
        end
    end
    return cell_vertices
end

function clear_mesh!(mesh::Mesh)
    for i = 1:mesh.NC_total
        mesh.particle_grid[i] = 0.0
    end
end

function getCellIndex(ix::Int64, iy::Int64, iz::Int64, NC::Vector{Int64})
    cellindex::Int64 = ix + (iy - 1) * NC[1] + (iz - 1) * NC[1] * NC[2]
    return cellindex
end

function getiXYZfromCellIndex(cellindex::Int64, N::Vector{Int64})
    ixyz = zeros(Int, 3)
    ixyz[3] = floor((cellindex - 1) / (N[1] * N[2])) + 1
    ixyz[2] = floor((cellindex - (ixyz[3] - 1) * N[1] * N[2] - 1) / N[1]) + 1
    ixyz[1] = (cellindex - (ixyz[3] - 1) * N[1] * N[2] - (ixyz[2] - 1) * N[1] - 1) % N[1] + 1
    return ixyz
end

function pbc_mesh(x_index::Int64, N_mesh::Int64)
    return mod(x_index - 1, N_mesh) + 1
end

function pbc_particle(position::Vector{Float64}, mesh::Mesh)
    position .-= mesh.boxsize .* floor.( position ./ mesh.boxsize)
    return position
end

function pbc_particle(position::Float64, boxsize::Float64)
    position -= boxsize * floor( position / boxsize)
    return position
end

function cloudincell(wrapped_pos::Array{Float64,1}, gridedge::Array{Float64,1}, meshnumber::Array{Int64,1})
    
    δx_lower = wrapped_pos[1] - floor(wrapped_pos[1] / edge_size[1]) * edge_size[1]
    δy_lower = wrapped_pos[2] - floor(wrapped_pos[2] / edge_size[2]) * edge_size[2]
    δz_lower = wrapped_pos[3] - floor(wrapped_pos[3] / edge_size[3]) * edge_size[3]

    δx_upper = edge_size[1] - δx_lower
    δy_upper = edge_size[1] - δy_lower
    δz_upper = edge_size[1] - δz_lower

end

function cloudincell!(position::Vector{Float64}, mesh::Mesh, cell_vertices::Array{Int64, 2}, i)
    # distribute particle mass over eight vertices in its cell

    # cell index = 1 + floor(x/cell_sizex) + floor(y/cell_sizey)*Nx + floor(z/cell_sizez)*Nx*Ny 
    # Bottom (x ->, y ^):
    # 3           4 
    # -------------
    # |  |  |  |  | 
    # |  |  |  |  | 
    # |  |  |  |  |
    # -------------
    # 1           2 
    # Top (+ z):
    # 7           8 
    # -------------
    # |  |  |  |  | 
    # |  |  |  |  | 
    # |  |  |  |  |
    # -------------
    # 5           6

    position_x = position[1] + 0.5 * mesh.boxsize[1]
    position_y = position[2] + 0.5 * mesh.boxsize[2]
    position_z = position[3] + 0.5 * mesh.boxsize[3]
    ix::Int64 = floor(position_x / mesh.edge_size[1]) + 1
    iy::Int64 = floor(position_y / mesh.edge_size[2]) + 1
    iz::Int64 = floor(position_z / mesh.edge_size[3]) + 1
    δx = position_x - floor(position_x / mesh.edge_size[1]) * mesh.edge_size[1]
    δy = position_y - floor(position_y / mesh.edge_size[2]) * mesh.edge_size[2]
    δz = position_z - floor(position_z / mesh.edge_size[3]) * mesh.edge_size[3]
    
    cell_index = getCellIndex(ix, iy, iz, mesh.NC)

    v1 = cell_vertices[cell_index, 1]
    v2 = cell_vertices[cell_index, 2]
    v3 = cell_vertices[cell_index, 3]
    v4 = cell_vertices[cell_index, 4]
    v5 = cell_vertices[cell_index, 5]
    v6 = cell_vertices[cell_index, 6]
    v7 = cell_vertices[cell_index, 7]
    v8 = cell_vertices[cell_index, 8]

    mesh.particle_grid[v1] += (mesh.edge_size[1] - δx) * (mesh.edge_size[2] - δy) * (mesh.edge_size[3] - δz) / prod(mesh.edge_size)
    mesh.particle_grid[v2] += (δx) * (mesh.edge_size[2] - δy) * (mesh.edge_size[3] - δz) / prod(mesh.edge_size)
    mesh.particle_grid[v3] += (mesh.edge_size[1] - δx) * (δy) * (mesh.edge_size[3] - δz) / prod(mesh.edge_size)                 
    mesh.particle_grid[v4] += (δx) * (δy) * (mesh.edge_size[3] - δz) / prod(mesh.edge_size)
    mesh.particle_grid[v5] += (mesh.edge_size[1] - δx) * (mesh.edge_size[2] - δy) * (δz) / prod(mesh.edge_size)
    mesh.particle_grid[v6] += (δx) * (mesh.edge_size[2] - δy) * (δz) / prod(mesh.edge_size)
    mesh.particle_grid[v7] += (mesh.edge_size[1] - δx) * (δy) * (δz) / prod(mesh.edge_size)
    mesh.particle_grid[v8] += (δx) * (δy) * (δz) / prod(mesh.edge_size)                                                      
    
end

function gaussian_distribution!(position::Vector{Float64}, mesh::Mesh, i, σ)
    # Adjusted position
    position_x = position[1] + 0.5 * mesh.boxsize[1]
    position_y = position[2] + 0.5 * mesh.boxsize[2]
    position_z = position[3] + 0.5 * mesh.boxsize[3]

    # Cell indices: round to the nearest cell center instead of always rounding down
    cell_index_x = round(Int64, position_x / mesh.edge_size[1])
    cell_index_y = round(Int64, position_y / mesh.edge_size[2])
    cell_index_z = round(Int64, position_z / mesh.edge_size[3])

    # Gaussian coefficient
    gaussian_coeff = 1 / (σ * sqrt(2π))^3

    current_mass = 0.0
    # Temporary grid to store updates
    particle_grid_update = zeros(Float64, mesh.NC_total) 
    gradient_field_update = zeros(Float64, mesh.NC_total, 3)

    # Distribute the mass
    for dx in -1:1, dy in -1:1, dz in -1:1
        # Adjusted cell indices with periodic boundary conditions
        adj_x = pbc_mesh(cell_index_x + dx, mesh.NC[1])
        adj_y = pbc_mesh(cell_index_y + dy, mesh.NC[2])
        adj_z = pbc_mesh(cell_index_z + dz, mesh.NC[3])

        # Distance from the particle to the current cell center
        dist_x = ((adj_x) * mesh.edge_size[1] - position_x) % mesh.boxsize[1]
        dist_y = ((adj_y) * mesh.edge_size[2] - position_y) % mesh.boxsize[2]
        dist_z = ((adj_z) * mesh.edge_size[3] - position_z) % mesh.boxsize[3]

        # Apply minimum image convention
        dist_x -= mesh.boxsize[1] * round(dist_x / mesh.boxsize[1])
        dist_y -= mesh.boxsize[2] * round(dist_y / mesh.boxsize[2])
        dist_z -= mesh.boxsize[3] * round(dist_z / mesh.boxsize[3])

        # Gaussian distribution
        gaussian_dist = gaussian_coeff * exp(-(dist_x^2 + dist_y^2 + dist_z^2) / (2*σ^2))
        
        current_mass += gaussian_dist
       
        # compute gradients
        grad_x = -(dist_x / σ^2) * gaussian_dist
        grad_y = -(dist_y / σ^2) * gaussian_dist
        grad_z = -(dist_z / σ^2) * gaussian_dist

        # Update the mass in the mesh
        current_index = getCellIndex(adj_x, adj_y, adj_z, mesh.NC)
        particle_grid_update[current_index] += gaussian_dist
        gradient_field_update[current_index, 1] += grad_x
        gradient_field_update[current_index, 2] += grad_y
        gradient_field_update[current_index, 3] += grad_z
    end

    # Normalize so that the total mass of a particle is 1
    for dx in -1:1, dy in -1:1, dz in -1:1
        # Adjusted cell indices with periodic boundary conditions
        adj_x = pbc_mesh(cell_index_x + dx, mesh.NC[1])
        adj_y = pbc_mesh(cell_index_y + dy, mesh.NC[2])
        adj_z = pbc_mesh(cell_index_z + dz, mesh.NC[3])
        current_index = getCellIndex(adj_x, adj_y, adj_z, mesh.NC)
        mesh.particle_grid[current_index] += particle_grid_update[current_index] / current_mass
        mesh.gradient_grid[current_index, 1] += gradient_field_update[current_index, 1] / current_mass
        mesh.gradient_grid[current_index, 2] += gradient_field_update[current_index, 2] / current_mass
        mesh.gradient_grid[current_index, 3] += gradient_field_update[current_index, 3] / current_mass
    end
end

function positional_density(particle_vec::Vector{Particle}, mesh::Mesh, cell_vertices::Array{Int64, 2})
# this function calculates the local interpolated density at the particle positions
    d_pos = zeros(Float64, length(particle_vec)) 
    pic = zeros(Float64, 8)

    for i in 1:length(particle_vec)
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

        dv1 = mesh.particle_grid[v1] * pic[1]
        dv2 = mesh.particle_grid[v2] * pic[2]
        dv3 = mesh.particle_grid[v3] * pic[3]
        dv4 = mesh.particle_grid[v4] * pic[4]
        dv5 = mesh.particle_grid[v5] * pic[5]
        dv6 = mesh.particle_grid[v6] * pic[6]
        dv7 = mesh.particle_grid[v7] * pic[7]
        dv8 = mesh.particle_grid[v8] * pic[8]
        
        d_pos[i] = dv1 + dv2 + dv3 + dv4 + dv5 + dv6 + dv7 + dv8
    end
    return d_pos
end

function Grad_DensVertex(vertex::Array{Int64, 3}, mesh::Mesh)
# ca<lculate density gradients on the standard density grid
    grid_grad_x = zeros(Float64, mesh.NC_total)
    grid_grad_y = zeros(Float64, mesh.NC_total)
    grid_grad_z = zeros(Float64, mesh.NC_total)

    for k = 1:mesh.NC[3]
        for j = 1:mesh.NC[2]
            for i = 1:mesh.NC[1]

                index = vertex[i, j, k] 

                i_plus = pbc_mesh(i + 1, mesh.NC[1])
                j_plus = pbc_mesh(j + 1, mesh.NC[2])
                k_plus = pbc_mesh(k + 1, mesh.NC[3])

                i_minus = pbc_mesh(i - 1, mesh.NC[1])
                j_minus = pbc_mesh(j - 1, mesh.NC[2])
                k_minus = pbc_mesh(k - 1, mesh.NC[3])
                
                grad_x = 0.5 / mesh.edge_size[1] * (mesh.particle_grid[vertex[i_plus, j, k]] - mesh.particle_grid[vertex[i_minus, j, k]])
                grad_y = 0.5 / mesh.edge_size[2] * (mesh.particle_grid[vertex[i, j_plus, k]] - mesh.particle_grid[vertex[i, j_minus, k]])
                grad_z = 0.5 / mesh.edge_size[3] * (mesh.particle_grid[vertex[i, j, k_plus]] - mesh.particle_grid[vertex[i, j, k_minus]])
                
                grid_grad_x[index] += grad_x
                grid_grad_y[index] += grad_y
                grid_grad_z[index] += grad_z
            end
        end
    end
    return grid_grad_x, grid_grad_y, grid_grad_z
end

function Grad_staggered_lattice!(vertex::Array{Int64, 3}, mesh::Mesh)
# calculate the density gradients on a staggered lattice

    x_max = mesh.NC[1] + 1
    y_max = mesh.NC[2] + 1
    z_max = mesh.NC[3] + 1

    for k = 1:mesh.NC[3]
        for j = 1:mesh.NC[2]
            for i = 1:mesh.NC[1]

                current_index = vertex[i, j, k]
                
                if i == 1
                    iv_minus = vertex[mesh.NC[1], j, k]
                else
                    iv_minus = vertex[i - 1, j, k]
                end
                
                if i == mesh.NC[1]
                    iv_plus = vertex[1, j, k]
                else 
                    iv_plus = vertex[i + 1, j, k]
                end

                if j == 1
                    jv_minus = vertex[i, mesh.NC[2], k]
                else
                    jv_minus = vertex[i, j - 1, k]
                end

                if j == mesh.NC[2]
                    jv_plus = vertex[i, 1, k]
                else
                    jv_plus = vertex[i, j + 1, k]
                end

                if k == 1
                    kv_minus = vertex[i, j, mesh.NC[3]]
                else
                    kv_minus = vertex[i, j, k - 1]
                end

                if k == mesh.NC[3]
                    kv_plus = vertex[i, j, 1]
                else
                    kv_plus = vertex[i, j, k + 1]
                end

                mesh.gradient_grid[current_index, 1] = ( mesh.particle_grid[iv_plus] - mesh.particle_grid[current_index] ) / mesh.edge_size[1]
                mesh.gradient_grid[current_index, 2] = ( mesh.particle_grid[current_index] - mesh.particle_grid[iv_minus] ) / mesh.edge_size[1]
                mesh.gradient_grid[current_index, 3] = ( mesh.particle_grid[jv_plus] - mesh.particle_grid[current_index] ) / mesh.edge_size[2]
                mesh.gradient_grid[current_index, 4] = ( mesh.particle_grid[current_index] - mesh.particle_grid[jv_minus] ) / mesh.edge_size[2]
                mesh.gradient_grid[current_index, 5] = ( mesh.particle_grid[kv_plus] - mesh.particle_grid[current_index] ) / mesh.edge_size[3]
                mesh.gradient_grid[current_index, 6] = ( mesh.particle_grid[current_index] - mesh.particle_grid[kv_minus] ) / mesh.edge_size[3]
            end
        end
    end
end

# short pbcs for vertex matrix, but doesn't work in this form,
# if you want only your array to go from 1 to N[1]*N[2]*N[3] 
#for k in 1 mesh.NV[3]
    #    for j in 1:mesh.NV[2]
    #        vertex[mesh.NV[1] + 1, j, k] = vertex[1, j, k]
    #    end
    #end
    #for k in 1 mesh.NV[3]
    #    for i in 1:mesh.NV[1]
    #        vertex[i, mesh.NV[2] + 1, k] = vertex[i, 1, k]
    #    end
    #end
    #for j in 1 mesh.NV[2]
    #    for i in 1:mesh.NV[1]
    #        vertex[i, j, mesh.NV[3] + 1] = vertex[i, j, 1]
    #    end
    #end