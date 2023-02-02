struct Grid
    num_vertices::Int
    num_dirichlet_vertices::Int
    neumann_vertices::Set{Int}
    triangles::Set{Set{Int}}
    x0::Vector{Float64}

    function Grid(num_vertices, num_dirichlet_vertices, neumann_vertices, triangles, x0)
        num_dirichlet_vertices <= num_vertices || throw(ArgumentError("the number of dirichlet vertices cannot be larger than the total number of vertices"))
        vertices = Set(1:num_vertices)
        neumann_vertices ⊆ vertices || throw(ArgumentError("neumann_vertices must be a subset of all vertices"))
        dirichlet_vertices = Set(1:num_dirichlet_vertices)
        isempty(neumann_vertices ∩ dirichlet_vertices) || throw(ArgumentError("vertex cannot be both in dirichlet_vertices and neumann_vertices")) 

        for t in triangles
            length(t) == 3 || throw(ArgumentError("triangles must contain exactly three vertices"))
            t ⊆ vertices || throw(ArgumentError("invalid triangle indices"))
        end

        not_covered_vertices = Set(1:num_vertices)
        for t in triangles
            setdiff!(not_covered_vertices, t)
        end
        not_covered_string = ""
        for i in not_covered_vertices
            not_covered_string *= ", $i"
        end
        not_covered_string = not_covered_string[3:end]

        isempty(not_covered_vertices) || throw(ArgumentError("vertices {$not_covered_string} not covered by triangles"))

        length(x0) == num_vertices || throw(ArgumentError("x0 must have length $num_vertices (vector of length $(length(x0)) was provided)"))

        new(num_vertices, num_dirichlet_vertices, neumann_vertices, triangles, x0)
    end
end

Grid(4, 2, Set([3, 4]), Set([Set([1, 2, 3]), Set([2, 3, 4])]), ones(5))
