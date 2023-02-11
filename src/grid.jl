@enum VertexType external boundary internal undetermined

mutable struct Vertex
    x::Float64
    y::Float64
    type::VertexType

    Vertex(x, y, type=undetermined) = new(x, y, type)
end

Triangle = Tuple{Vertex,Vertex,Vertex}

function isosceles_right_triangulation(width::Float64, height::Float64)::Vector{Triangle}
    num_squares_hor = ceil(Int, width)
    num_squares_ver = ceil(Int, height)
    num_vertices_hor = num_squares_hor + 1
    num_vertices_ver = num_squares_ver + 1
    vertices = Matrix{Vertex}(undef, num_vertices_hor, num_vertices_ver)
    for i in 1:num_vertices_hor
        for j in 1:num_vertices_ver
            vertices[i, j] = Vertex(i - 1, j - 1)
        end
    end

    num_triangles = 2 * num_squares_hor * num_squares_ver
    triangles = Vector{Triangle}(undef, num_triangles)
    k = 1
    for i in 1:num_squares_hor
        for j in 1:num_squares_ver
            triangles[k] = (vertices[i, j], vertices[i+1, j], vertices[i+1, j+1])
            triangles[k+1] = (vertices[i, j], vertices[i, j+1], vertices[i+1, j+1])
            k += 2
        end
    end

    triangles
end

struct SimulationGrid
    num_vertices::Int
    num_dirichlet_vertices::Int
    num_neumann_vertices::Int
    x::Vector{Float64}
    y::Vector{Float64}
    triangles::Vector{Tuple{Int,Int,Int}}
end

function SimulationGrid(width, height, triangulation, isinternal, isdirichlet)
    covering_triangulation = triangulation(width, height)

    vertex_isinternal(v::Vertex) = isinternal(v.x, v.y)
    triangle_isinternal(triangle::Triangle) = all(vertex_isinternal(v) for v in triangle)
    vertex_isdirichlet(v::Vertex) = isdirichlet(v.x, v.y)

    # determine if a vertex is internal, external, or on the boundary
    for T in covering_triangulation
        containing_triangle_isinternal = triangle_isinternal(T)
        for v in T
            update_vertex_type(v, containing_triangle_isinternal)
        end
    end

    # group vertices by their type
    internal_vertices = Set{Vertex}()
    dirichlet_vertices = Set{Vertex}()
    neumann_vertices = Set{Vertex}()
    for T in covering_triangulation
        for v in T
            if v.type == internal
                push!(internal_vertices, v)
            elseif v.type == boundary
                if vertex_isdirichlet(v)
                    push!(dirichlet_vertices, v)
                else
                    push!(neumann_vertices, v)
                end
            end
        end
    end

    internal_triangles = Vector{Triangle}()
    for T in covering_triangulation
        if triangle_isinternal(T)
            push!(internal_triangles, T)
        end
    end
end

function update_vertex_type(v::Vertex, containing_triangle_isinternal::Bool)
    if v.type == boundary
        return
    end

    if v.type == undetermined
        v.type = containing_triangle_isinternal ? internal : external
        return
    end

    if (v.type == external && containing_triangle_isinternal) || (v.type == internal && !containing_triangle_isinternal)
        v.type = boundary
    end
end

radius = 3.0
dirichlet_arc_angle = pi / 4
width = 2 * radius
height = 2 * radius
triangulation = isosceles_right_triangulation
isinternal(x, y) = (x - radius)^2 + (y - radius)^2 <= radius^2
isdirichlet(x, y) = x - radius <= -cos(dirichlet_arc_angle)
SimulationGrid(width, height, triangulation, isinternal, isdirichlet)