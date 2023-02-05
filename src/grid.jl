@enum VertexType external boundary internal

struct Point
    x::Float64
    y::Float64
end

struct Triangle
    a::Point
    b::Point
    c::Point
end

struct IsoscelesRightTriangulation
    width::Float64
    height::Float64
    vertices::Matrix{Point}
    triangles::Vector{Triangle}

    function IsoscelesRightTriangulation(width::Float64, height::Float64)
        num_squares_hor = ceil(Int, width)
        num_squares_ver = ceil(Int, height)
        num_vertices_hor = num_squares_hor + 1
        num_vertices_ver = num_squares_ver + 1
        vertices = Matrix{Point}(undef, num_vertices_hor, num_vertices_ver)
        for i in 1:num_vertices_hor
            for j in 1:num_vertices_ver
                vertices[i, j] = Point(i - 1, j - 1)
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

        new(width, height, vertices)
    end
end
