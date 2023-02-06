@enum VertexType external boundary internal undetermined

mutable struct Vertex
    x::Float64
    y::Float64
    type::VertexType

    Vertex(x, y, type=undetermined) = new(x, y, type)
end

struct Triangle
    a::Vertex
    b::Vertex
    c::Vertex
end

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
            triangles[k] = Triangle(vertices[i, j], vertices[i+1, j], vertices[i+1, j+1])
            triangles[k+1] = Triangle(vertices[i, j], vertices[i, j+1], vertices[i+1, j+1])
            k += 2
        end
    end

    triangles
end

for T in isosceles_right_triangulation(2.3, 3.1)
    println(T)
end
