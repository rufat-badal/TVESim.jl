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
    # order: dirichlet vertices, neumann vertices, internal vertices
    dirichlet_vertices
    neumann_vertices
    boundary_vertices
    x::Vector{Float64}
    y::Vector{Float64}
    θ::Vector{Float64}
    triangles::Vector{Tuple{Int,Int,Int}}
    area_factors::Vector{Float64}

    function SimulationGrid(
        num_vertices, num_dirichlet_vertices, num_neumann_vertices,
        x, y, θ, triangles
    )
        dirichlet_vertices = 1:num_dirichlet_vertices
        neumann_vertices = num_dirichlet_vertices+1:num_dirichlet_vertices+num_neumann_vertices
        boundary_vertices = 1:num_dirichlet_vertices+num_neumann_vertices
        area_factors = [
            area_factor([x[i1], y[i1]], [x[i2], y[i2]], [x[i3], y[i3]])
            for (i1, i2, i3) in triangles
        ]

        new(
            num_vertices, num_dirichlet_vertices, num_neumann_vertices,
            dirichlet_vertices, neumann_vertices, boundary_vertices,
            x, y, θ, triangles, area_factors
        )
    end
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

    SimulationGrid(
        length(vertices), length(dirichlet_vertices), length(neumann_vertices),
        x, y, θ, triangles
    )
end

function SimulationGrid(boundary_points, initial_temperature, isdirichlet)
    min_angle = 10

    min_segment_length_sqrd = Inf
    for (p, q) in zip(eachcol(boundary_points[:, 1:end-1]), eachcol(boundary_points[:, 2:end]))
        segment_length_sqrd = (p[1] - q[1])^2 + (p[2] - q[2])^2
        if segment_length_sqrd < min_segment_length_sqrd
            min_segment_length_sqrd = segment_length_sqrd
        end
    end
    max_area = sqrt(3) / 4 * min_segment_length_sqrd # area of an equilateral triangle

    num_initial_boundary_points = size(boundary_points, 2)
    triin = Triangulate.TriangulateIO()
    triin.pointlist = boundary_points
    triin.segmentlist = [1:num_initial_boundary_points [2:num_initial_boundary_points...; 1]]'
    triout, _ = Triangulate.triangulate("pa$(max_area)q$(min_angle)Q", triin)

    points = triout.pointlist
    markers = triout.pointmarkerlist
    num_boundary_vertices = 0
    num_dirichlet_vertices = 0
    for (p, m) in zip(eachcol(points), markers)
        if m == 0
            continue
        end

        num_boundary_vertices += 1
        x, y = p
        if isdirichlet(x, y)
            num_dirichlet_vertices += 1
        end
    end
    vertex_new_id = Dict{Int,Int}()
    dirichlet_points = Vector{Vector{Float64}}()
    neumann_points = Vector{Vector{Float64}}()
    internal_points = Vector{Vector{Float64}}()
    next_dirichlet_id = 1
    next_neumann_id = num_dirichlet_vertices + 1
    next_internal_id = num_boundary_vertices + 1
    for (i, (p, m)) in enumerate(zip(eachcol(points), markers))
        if m == 1
            x, y = p
            if isdirichlet(x, y)
                push!(dirichlet_points, p)
                vertex_new_id[i] = next_dirichlet_id
                next_dirichlet_id += 1
            else
                push!(neumann_points, p)
                vertex_new_id[i] = next_neumann_id
                next_neumann_id += 1
            end
        else
            push!(internal_points, p)
            vertex_new_id[i] = next_internal_id
            next_internal_id += 1
        end
    end

    x = [p[1] for p in dirichlet_points]
    append!(x, [p[1] for p in neumann_points])
    append!(x, [p[1] for p in internal_points])
    y = [p[2] for p in dirichlet_points]
    append!(y, [p[2] for p in neumann_points])
    append!(y, [p[2] for p in internal_points])
    θ = initial_temperature * ones(1:size(points, 2))
    triangles = [
        (vertex_new_id[i], vertex_new_id[j], vertex_new_id[k])
        for (i, j, k) in eachcol(triout.trianglelist)
    ]
    area_factors = [
        area_factor([x[i1], y[i1]], [x[i2], y[i2]], [x[i3], y[i3]])
        for (i1, i2, i3) in triangles
    ]

    SimulationGrid(
        size(points, 2), num_dirichlet_vertices, num_boundary_vertices - num_dirichlet_vertices,
        x, y, θ, triangles
    )
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
    num_horizontal_pixels = 2500 # in pixels!
    strokewidth = 2.5 / 1000 * num_horizontal_pixels
    markersize = 2 / 100 * num_horizontal_pixels

    min_x, max_x = minimum(grid.x), maximum(grid.x)
    min_y, max_y = minimum(grid.y), maximum(grid.y)
    width = max_x - min_x
    height = max_y - min_y
    aspect = width / height
    plot_height = num_horizontal_pixels / aspect
    fig = CairoMakie.Figure(resolution=(num_horizontal_pixels, plot_height))
    length_pixel = width / num_horizontal_pixels
    padding = (strokewidth + markersize / 2) * length_pixel
    ax = CairoMakie.Axis(
        fig[1, 1],
        limits=(
            min_x - padding, max_x + padding,
            min_y - padding, max_y + padding),
        aspect=aspect)
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
            color=grid.θ, colormap=:plasma, colorrange=temp_range, strokewidth=strokewidth, shading=true)
    else
        CairoMakie.mesh!(
            vertices, faces, color=grid.θ,
            colormap=:plasma, colorrange=temp_range)
    end

    dirichlet_ids = 1:grid.num_dirichlet_vertices
    CairoMakie.scatter!(
        grid.x[dirichlet_ids], grid.y[dirichlet_ids],
        color=:red, markersize=markersize, strokewidth=strokewidth
    )
    neumann_ids = grid.num_dirichlet_vertices+1:grid.num_dirichlet_vertices+grid.num_neumann_vertices
    CairoMakie.scatter!(
        grid.x[neumann_ids], grid.y[neumann_ids],
        color=:green, markersize=markersize, strokewidth=strokewidth
    )

    fig
end

function plot(triangulation::Vector{Triangle})
    num_horizontal_pixels = 1000
    strokewidth = 1
    padding_perc = 0.5

    vertices = Set{Vertex}()
    for T in triangulation
        for v in T
            push!(vertices, v)
        end
    end

    x = [v.x for v in vertices]
    y = [v.y for v in vertices]
    vertex_coords = [x y]

    vertex_to_id = Dict(v => i for (i, v) in enumerate(vertices))
    faces = Matrix{Int}(undef, length(triangulation), 3)
    for (i, T) in enumerate(triangulation)
        for (j, vertex) in enumerate(T)
            faces[i, j] = vertex_to_id[vertex]
        end
    end

    min_x, max_x = minimum(x), maximum(x)
    min_y, max_y = minimum(y), maximum(y)
    width = max_x - min_x
    height = max_y - min_y
    aspect = width / height
    plot_height = num_horizontal_pixels / aspect
    fig = CairoMakie.Figure(resolution=(num_horizontal_pixels, plot_height))
    horizontal_padding = padding_perc / 100 * width
    vertical_padding = padding_perc / 100 * height
    ax = CairoMakie.Axis(
        fig[1, 1],
        limits=(
            min_x - horizontal_padding, max_x + vertical_padding,
            min_y - horizontal_padding, max_y + horizontal_padding),
        aspect=aspect)
    CairoMakie.hidedecorations!(ax)
    CairoMakie.hidespines!(ax)
    CairoMakie.poly!(vertex_coords, faces, color=:transparent, strokewidth=strokewidth, shading=true)

    fig
end

function plot(triangulation::Vector{Triangle}, width, height)
    fig = plot(triangulation)
    CairoMakie.poly!(
        CairoMakie.Point2f[(0, 0), (width, 0), (width, height), (0, height)],
        color=:transparent, strokecolor=:blue, strokewidth=2)
    fig
end

function plot(triangulation::Vector{Triangle}, side_length)
    plot(triangulation, side_length, side_length)
end

function circle_boundary_points(radius, num_points)
    θ = range(0, 2π, length=num_points + 1)
    x = radius * (1 .+ cos.(θ))
    pop!(x)
    y = radius * (1 .+ sin.(θ))
    pop!(y)

    [x y]'
end

area_factor(a, b, c) = abs(LinearAlgebra.det([b - a c - a]))

function test_crystalline_simulationgrid(triangulation)
    radius = 20
    initial_temperature = 0
    dirichlet_arc_angle = pi / 4
    width = 2 * radius
    height = 2 * radius
    isinternal(x, y) = (x - radius)^2 + (y - radius)^2 <= radius^2
    isdirichlet(x, y) = x - radius <= -radius * cos(dirichlet_arc_angle)
    grid = SimulationGrid(width, height, initial_temperature, triangulation, isinternal, isdirichlet)
    plot(grid, (0.0, 1.0), show_edges=true)
end

function test_cdt_simulationgrid()
    radius = 20
    initial_temperature = 0.0
    num_boundary_points = 100
    dirichlet_arc_angle = 45
    isdirichlet(x, y) = x - radius <= -radius * cos(dirichlet_arc_angle / 360 * 2pi)
    grid = SimulationGrid(circle_boundary_points(radius, num_boundary_points), initial_temperature, isdirichlet)
    plot(grid, (0.0, 1.0), show_edges=true)
end
