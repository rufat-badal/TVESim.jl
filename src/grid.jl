import CairoMakie

@enum VertexType external boundary internal undetermined

const triangle_side_length = 1.0

mutable struct Vertex
    x::Float64
    y::Float64
    type::VertexType

    Vertex(x, y, type=undetermined) = new(x, y, type)
end

Triangle = Tuple{Vertex,Vertex,Vertex}

function isosceles_right_triangulation(width::Number, height::Number)
    num_squares_hor = ceil(Int, width / triangle_side_length)
    num_squares_ver = ceil(Int, height / triangle_side_length)
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

function equilateral_triangulation(width::Number, height::Number)
    triangle_height = sqrt(3) / 2 * triangle_side_length
    num_layers = ceil(Int, height / triangle_height)
    num_upward_triangles_per_layer = ceil(Int, (width - triangle_side_length / 2) / triangle_side_length) + 1
    vertices = Matrix{Vertex}(undef, num_upward_triangles_per_layer + 1, num_layers + 1)
    for i in 1:(num_upward_triangles_per_layer+1)
        for j in 1:(num_layers+1)
            vertices[i, j] = Vertex(-(j % 2) * triangle_side_length / 2 + (i - 1) * triangle_side_length, (j - 1) * triangle_height)
        end
    end

    num_triangles = 2 * num_upward_triangles_per_layer * num_layers
    triangles = Vector{Triangle}(undef, num_triangles)
    k = 1
    for j in 1:num_layers
        for i in 1:num_upward_triangles_per_layer
            if j % 2 == 0
                triangles[k] = (vertices[i, j], vertices[i+1, j+1], vertices[i, j+1])
                triangles[k+1] = (vertices[i, j], vertices[i+1, j], vertices[i+1, j+1])
            else
                triangles[k] = (vertices[i, j], vertices[i+1, j], vertices[i, j+1])
                triangles[k+1] = (vertices[i+1, j], vertices[i+1, j+1], vertices[i, j+1])
            end

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
    θ::Vector{Float64}
    triangles::Vector{Tuple{Int,Int,Int}}
end

function SimulationGrid(width, height, initial_temperature, triangulation, isinternal, isdirichlet)
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

    # vertices are grouped as follows: dirichlet, neumann, internal
    vertices = Vector{Vertex}()
    append!(vertices, dirichlet_vertices)
    append!(vertices, neumann_vertices)
    append!(vertices, internal_vertices)
    x = [v.x for v in vertices]
    y = [v.y for v in vertices]
    vertex_to_id = Dict(v => i for (i, v) in enumerate(vertices))
    triangles = [(vertex_to_id[T[1]], vertex_to_id[T[2]], vertex_to_id[T[3]]) for T in internal_triangles]

    θ = [initial_temperature for _ in vertices]

    SimulationGrid(length(vertices), length(dirichlet_vertices), length(neumann_vertices), x, y, θ, triangles)
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

function plot(grid::SimulationGrid, temp_range; show_edges=false)
    plot_width = 1000 # in pixels!
    fontsize = 20
    min_x, max_x = minimum(grid.x), maximum(grid.x)
    min_y, max_y = minimum(grid.y), maximum(grid.y)
    width = max_x - min_x
    height = max_y - min_y
    aspect = width / height
    plot_height = plot_width / aspect
    fig = CairoMakie.Figure(resolution=(plot_width, plot_height), fontsize=fontsize)
    ax = CairoMakie.Axis(fig[1, 1], limits=(min_x, max_x, min_y, max_y), aspect=aspect)
    CairoMakie.hidedecorations!(ax)
    CairoMakie.hidespines!(ax)

    vertices = [grid.x grid.y]
    faces = Matrix{Int}(undef, length(grid.triangles), 3)
    for (i, T) in enumerate(grid.triangles)
        for (j, vertex) in enumerate(T)
            faces[i, j] = vertex
        end
    end

    if show_edges
        CairoMakie.poly!(
            vertices, faces,
            color=grid.θ, colormap=:plasma, colorrange=temp_range, strokewidth=1, shading=true)
    else
        CairoMakie.mesh!(
            vertices, faces, color=grid.θ,
            colormap=:plasma, colorrange=temp_range)
    end

    fig
end


radius = 10
dirichlet_arc_angle = pi / 4
width = 2 * radius
height = 2 * radius
triangulation = equilateral_triangulation
isinternal(x, y) = (x - radius)^2 + (y - radius)^2 <= radius^2
isdirichlet(x, y) = x - radius <= -cos(dirichlet_arc_angle)
grid = SimulationGrid(width, height, 0.0, triangulation, isinternal, isdirichlet)
plot(grid, (0.0, 1.0), show_edges=true)