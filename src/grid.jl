struct Grid
    num_vertices::Int
    triangles::Set{Set{Int}}
    dirichlet_vertices::Set{Int}
    neumann_vertices::Set{Int}
    x0::Vector{Float64}

    function Grid(n, triangles, dirichlet_vertices, neumann_vertices, x0)
        for t in triangles
            length(t) == 3 || throw(ArgumentError("triangles must contain exactly three vertices (provided triangle with $(length(t)) provided)"))
            for i in t
                1 <= i <= n || throw(ArgumentError("vertices must be between 1 and $n (vertex $i provided)"))
            end
        end

        isempty(dirichlet_vertices ∩ neumann_vertices) || throw(ArgumentError("vertex cannot be in both dirichlet_vertices and neumann_vertices"))

        not_covered_vertices = Set(1:n)
        for t in triangles
            setdiff!(not_covered_vertices, t)
        end
        isempty(not_covered_vertices) || throw(ArgumentError("vertices $not_covered_vertices not covered by triangles"))

        length(x0) == n || throw(ArgumentError("x0 must have length $n (vector of length $(length(x0)) was provided)"))

        new(n, triangles, dirichlet_vertices, neumann_vertices, x0)
    end
end
