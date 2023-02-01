struct Grid
    num_vertices::Int
    triangles::Vector{Tuple{Int, Int, Int}}
    edges::Vector{Tuple{Int, Int}}
    dirichlet_vertices::Vector{Int}
    neumann_vertices::Vector{Int}
end