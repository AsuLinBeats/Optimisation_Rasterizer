#pragma once

#include <vector>
#include <iostream>
#include <immintrin.h>
#include "Vec4_Optimised.h"
//#include "vec4.h"
// #include "matrix.h"
#include "matrix_optimised.h"
//#include "colour.h"
#include "colour_optimised.h"
#include<cmath>
#include<thread>
#include<mutex>
#include<array>

// Represents a vertex in a 3D mesh, including its position, normal, and color
struct Vertex {
    vec4 p;         // Position of the vertex in 3D space
    vec4 normal;    // Normal vector for the vertex
    colour rgb;     // Color of the vertex

    Vertex() : p(), normal(), rgb() {}  // Default constructor
    Vertex(const vec4& pos, const vec4& norm, const colour& col)
        : p(pos), normal(norm), rgb(col) {
    }

};

// Stores indices of vertices that form a triangle in a mesh
struct triIndices {
    unsigned int v[3]; // Indices into the vertex array
    triIndices() : v{ 0, 0, 0 } {} // add default constructor
    // Constructor to initialize the indices of a triangle
    triIndices(unsigned int v1, unsigned int v2, unsigned int v3) {
        v[0] = v1;
        v[1] = v2;
        v[2] = v3;
    }
};

// Class representing a 3D mesh made up of vertices and triangles
class Mesh {
public:
    colour col;       // Uniform color for the mesh
    float kd;         // Diffuse reflection coefficient
    float ka;         // Ambient reflection coefficient
    matrix world;     // Transformation matrix for the mesh
    std::vector<Vertex> vertices;       // List of vertices in the mesh
    std::vector<triIndices> triangles;  // List of triangles in the mesh

    void setColour(colour _c, float _ka, float _kd) {
        col = _c;
        ka = _ka;
        kd = _kd;
    }

    // Default constructor initializes default color and reflection coefficients
    Mesh() {
        col.set(1.0f, 1.0f, 1.0f);
        ka = kd = 0.75f;
    }

    void addVertex(const vec4& vertex, const vec4& normal) {
        Vertex v = { vertex, normal, col };
        vertices.push_back(v);
    }

    void addTriangle(int v1, int v2, int v3) {
        triangles.emplace_back(v1, v2, v3);
    }

    // Display the vertices and triangles of the mesh
    void display() const {
        std::cout << "Vertices and Normals:\n";
        for (size_t i = 0; i < vertices.size(); ++i) {
            std::cout << i << ": Vertex (" << vertices[i].p[0] << ", " << vertices[i].p[1] << ", " << vertices[i].p[2] << ", " << vertices[i].p[3] << ")"
                << " Normal (" << vertices[i].normal[0] << ", " << vertices[i].normal[1] << ", " << vertices[i].normal[2] << ", " << vertices[i].normal[3] << ")\n";
        }

        std::cout << "\nTriangles:\n";
        for (const auto& t : triangles) {
            std::cout << "(" << t.v[0] << ", " << t.v[1] << ", " << t.v[2] << ")\n";
        }
    }


    static Mesh makeRectangle(float x1, float y1, float x2, float y2) {
        //! multithread here
        Mesh mesh;
        mesh.vertices.resize(4);
        mesh.triangles.resize(2);
        // init
        __m128 xValues = _mm_setr_ps(x1, x2, x2, x1);
        __m128 yValues = _mm_setr_ps(y1, y1, y2, y2);
        __m128 zValues = _mm_setzero_ps();
        alignas(16) float xa[4], ya[4];
        _mm_store_ps(xa, xValues);
        _mm_store_ps(ya, yValues);

        // normal calculation
        const vec4 v1 = vec4(xa[0], ya[0], 0.0f, 1.0f);
        const vec4 v2 = vec4(xa[1], ya[1], 0.0f, 1.0f);
        const vec4 v4 = vec4(xa[3], ya[3], 0.0f, 1.0f);
        // subtraction, normal, cross has been optimised in vec4.h, so just keep them.
        const vec4 edge1 = v2 - v1;
        const vec4 edge2 = v4 - v1;

        vec4 normal = vec4::cross(edge1, edge2);
        normal.normalise();

        //set vertices
       //  const __m128 z_combined = _mm_unpacklo_ps(zValues); // [0,1,0,1]
        for (int i = 0; i < 4; ++i) {
            const __m128 pos = _mm_setr_ps(xa[i], ya[i], 0.0f, 1.0f);
            mesh.vertices[i].p = vec4(pos);
            mesh.vertices[i].normal = normal;
        }

        // add triangle
        mesh.addTriangle(0, 2, 1);
        mesh.addTriangle(0, 3, 2);

        return mesh;
    }


    static Mesh makeCube(float size) {
        Mesh mesh;
        constexpr int numFaces = 6;
        constexpr int vertsFace = 4;
        constexpr int trisFace = 2;

        mesh.vertices.resize(numFaces * vertsFace);
        mesh.triangles.resize(numFaces * trisFace);

        float halfSize = size / 2.0f;

        vec4 positions[8] = {
            vec4(-halfSize, -halfSize, -halfSize, 1.0f),
            vec4(halfSize, -halfSize, -halfSize, 1.0f),
            vec4(halfSize, halfSize, -halfSize, 1.0f),
            vec4(-halfSize, halfSize, -halfSize, 1.0f),
            vec4(-halfSize, -halfSize, halfSize, 1.0f),
            vec4(halfSize, -halfSize, halfSize, 1.0f),
            vec4(halfSize, halfSize, halfSize, 1.0f),
            vec4(-halfSize, halfSize, halfSize, 1.0f)
        };

        vec4 normals[6] = {
            vec4(0, 0, -1, 0),  // front
            vec4(0, 0, 1, 0),   // back
            vec4(-1, 0, 0, 0),  // left
            vec4(1, 0, 0, 0),   // right
            vec4(0, -1, 0, 0),  // bottom
            vec4(0, 1, 0, 0)    // top
        };

        constexpr int faceIndices[6][4] = {
            {1, 0, 3, 2},  // front
            {4, 5, 6, 7},  // back
            {3, 0, 4, 7},  // left
            {5, 1, 2, 6},  // right
            {0, 1, 5, 4},  // bottom
            {2, 3, 7, 6}   // top
        };

        std::vector<std::thread> threads;
        for (int i = 0; i < numFaces; ++i) {
            threads.emplace_back([i, &mesh, &positions, &normals, &faceIndices]() {
                std::array<Vertex, 4> localVertices;
               
                std::array<triIndices, 2> localTriangles = {
                    triIndices(0, 0, 0),  
                    triIndices(0, 0, 0)
                };
                for (int j = 0; j < 4; ++j) {
                    localVertices[j] = Vertex{ positions[faceIndices[i][j]], normals[i], mesh.col };
                }

                const int vertBase = i * vertsFace;
                const int triBase = i * trisFace;

                localTriangles = {
                    triIndices(vertBase, vertBase + 2, vertBase + 1),
                    triIndices(vertBase, vertBase + 3, vertBase + 2)
                };

                for (int j = 0; j < 4; ++j) {
                    mesh.vertices[vertBase + j] = localVertices[j];
                }
                mesh.triangles[triBase] = localTriangles[0];
                mesh.triangles[triBase + 1] = localTriangles[1];
                });
        }

        for (auto& t : threads) t.join();
        return mesh;
    }

    static Mesh makeSphere(float radius, int latitudeDivisions, int longitudeDivisions) {
        Mesh mesh;
        if (latitudeDivisions < 2 || longitudeDivisions < 3) {
            throw std::invalid_argument("Latitude divisions must be >= 2 and longitude divisions >= 3");
        }

        mesh.vertices.clear();
        mesh.triangles.clear();

        // Create vertices
        for (int lat = 0; lat <= latitudeDivisions; ++lat) {
            float theta = M_PI * lat / latitudeDivisions;
            float sinTheta = std::sin(theta);
            float cosTheta = std::cos(theta);

            for (int lon = 0; lon <= longitudeDivisions; ++lon) {
                float phi = 2 * M_PI * lon / longitudeDivisions;
                float sinPhi = std::sin(phi);
                float cosPhi = std::cos(phi);

                vec4 position(
                    radius * sinTheta * cosPhi,
                    radius * sinTheta * sinPhi,
                    radius * cosTheta,
                    1.0f
                );

                vec4 normal = position;
                normal.normalise();
                normal[3] = 0.f;

                mesh.addVertex(position, normal);
            }
        }

        // Create indices for triangles
        for (int lat = 0; lat < latitudeDivisions; ++lat) {
            for (int lon = 0; lon < longitudeDivisions; ++lon) {
                int v0 = lat * (longitudeDivisions + 1) + lon;
                int v1 = v0 + 1;
                int v2 = (lat + 1) * (longitudeDivisions + 1) + lon;
                int v3 = v2 + 1;

                mesh.addTriangle(v0, v1, v2);
                mesh.addTriangle(v1, v3, v2);
            }
        }
        return mesh;
    }


 //! work, but slower

//static Mesh makeSphere(float radius, int latitudeDivisions, int longitudeDivisions) {
//    Mesh mesh;
//    if (latitudeDivisions < 2 || longitudeDivisions < 3) {
//        throw std::invalid_argument("Latitude divisions must be >= 2 and longitude divisions >= 3");
//    }
//
//    // Pre-calculate total number of vertices and triangles.
//    int totalLat = latitudeDivisions + 1;                     // rows of vertices
//    int verticesPerLat = longitudeDivisions + 1;                // vertices per row
//    int totalVertices = totalLat * verticesPerLat;
//    int totalTriangles = 2 * latitudeDivisions * longitudeDivisions;
//
//    // Pre-allocate the arrays.
//    mesh.vertices.resize(totalVertices);
//    mesh.triangles.resize(totalTriangles);
//
//    // Determine how many threads to use.
//    int numThreads = std::thread::hardware_concurrency();
//    if (numThreads <= 0)
//        numThreads = 2;
//    // We partition by vertex rows.
//    // We never want more threads than there are rows.
//    int effectiveThreads = min(numThreads, totalLat);
//
//    std::vector<std::thread> threads;
//    threads.reserve(effectiveThreads);
//
//    // Partition the latitude rows among threads.
//    // Use an even partition with remainder.
//    int base = totalLat / effectiveThreads;
//    int rem = totalLat % effectiveThreads;
//    int currentLat = 0;
//
//    for (int t = 0; t < effectiveThreads; t++) {
//        int latCount = base + (t < rem ? 1 : 0);
//        int latStart = currentLat;
//        int latEnd = latStart + latCount;  // latEnd is exclusive
//        currentLat = latEnd;
//
//        // Launch a thread that processes vertex rows [latStart, latEnd)
//        threads.emplace_back([=, &mesh]() {
//            // For each latitude row in this thread, compute the vertices.
//            // Note: For the spherical parameterization, theta ranges from 0 to PI,
//            // and we use latitude index 'lat' in [0, latitudeDivisions].
//            for (int lat = latStart; lat < latEnd; lat++) {
//                float theta = M_PI * lat / latitudeDivisions;
//                float sinTheta = std::sin(theta);
//                float cosTheta = std::cos(theta);
//                for (int lon = 0; lon < longitudeDivisions + 1; lon++) {
//                    float phi = 2 * M_PI * lon / longitudeDivisions;
//                    float sinPhi = std::sin(phi);
//                    float cosPhi = std::cos(phi);
//
//                    vec4 position(
//                        radius * sinTheta * cosPhi,
//                        radius * sinTheta * sinPhi,
//                        radius * cosTheta,
//                        1.0f
//                    );
//
//                    vec4 normal = position;
//                    normal.normalise();
//                    normal[3] = 0.f;
//
//                    // Compute the global vertex index.
//                    int index = lat * verticesPerLat + lon;
//                    mesh.vertices[index] = Vertex(position, normal, mesh.col);
//                }
//            }
//
//            // Next, generate triangles for all rows that this thread covers
//            // except the very last row because triangles use row 'lat' and 'lat+1'.
//            // Triangle indices are computed globally:
//            // For a given (lat, lon), let:
//            //   v0 = lat * verticesPerLat + lon,
//            //   v1 = v0 + 1,
//            //   v2 = (lat+1) * verticesPerLat + lon,
//            //   v3 = v2 + 1.
//            // Then create two triangles: (v0, v1, v2) and (v1, v3, v2).
//            // Note that valid lat for triangles is 0 <= lat < latitudeDivisions.
//            int localTriangleStart = latStart;  // first row this thread might generate triangles for
//            // But clamp the upper bound to latitudeDivisions.
//            int latForTrianglesEnd = min(latEnd, latitudeDivisions);
//            for (int lat = localTriangleStart; lat < latForTrianglesEnd; lat++) {
//                for (int lon = 0; lon < longitudeDivisions; lon++) {
//                    int v0 = lat * verticesPerLat + lon;
//                    int v1 = v0 + 1;
//                    int v2 = (lat + 1) * verticesPerLat + lon;
//                    int v3 = v2 + 1;
//                    // Compute the global triangle array offset:
//                    int triIndex = lat * (2 * longitudeDivisions) + (lon * 2);
//                    mesh.triangles[triIndex] = triIndices(v0, v1, v2);
//                    mesh.triangles[triIndex + 1] = triIndices(v1, v3, v2);
//                }
//            }
//            });
//    }
//
//    // Join all threads.
//    for (auto& th : threads) {
//        th.join();
//    }
//    return mesh;
//}

//! miss sphere.....

//static Mesh makeSphere(float radius, int latitudeDivisions, int longitudeDivisions) {
//    Mesh mesh;
//    if (latitudeDivisions < 2 || longitudeDivisions < 3) {
//        throw std::invalid_argument("Latitude divisions must be >= 2 and longitude divisions >= 3");
//    }
//
//    // We want to generate (latitudeDivisions+1) rows of vertices.
//    int totalRows = latitudeDivisions + 1;
//    // We'll generate exactly 'longitudeDivisions' vertices per row (seam is wrapped via modulus)
//    int totalVertices = totalRows * longitudeDivisions;
//    // The number of triangles is: for each of the (totalRows-1) bands, 
//    // for each vertex in a row, we create two triangles.
//    int totalTriangles = (totalRows - 1) * longitudeDivisions * 2;
//
//    // Reserve (optional) but note: we will later merge local data.
//    mesh.vertices.reserve(totalVertices);
//    mesh.triangles.reserve(totalTriangles);
//
//    // Determine the number of threads and partition rows among threads.
//    int numThreads = std::thread::hardware_concurrency();
//    if (numThreads == 0)
//        numThreads = 4;
//    int rowsPerThread = totalRows / numThreads;
//    if (rowsPerThread == 0) {
//        rowsPerThread = 1;
//        numThreads = totalRows;
//    }
//
//    std::vector<std::thread> threads;
//    std::mutex meshMutex;  // To protect merging into the mesh
//
//    for (int t = 0; t < numThreads; ++t) {
//        int startRow = t * rowsPerThread;
//        int endRow = (t == numThreads - 1) ? totalRows : startRow + rowsPerThread;
//        threads.emplace_back([=, &mesh, &meshMutex]() {
//            std::vector<Vertex> localVertices;
//            std::vector<triIndices> localTriangles;
//
//            // Generate vertices for rows [startRow, endRow)
//            for (int lat = startRow; lat < endRow; ++lat) {
//                float theta = M_PI * lat / latitudeDivisions; // lat=0 => theta=0, lat=latitudeDivisions => theta=M_PI
//                float sinTheta = sin(theta);
//                float cosTheta = cos(theta);
//                for (int lon = 0; lon < longitudeDivisions; ++lon) {
//                    float phi = 2 * M_PI * lon / longitudeDivisions;
//                    vec4 position(
//                        radius * sinTheta * cos(phi),
//                        radius * sinTheta * sin(phi),
//                        radius * cosTheta,
//                        1.0f
//                    );
//                    vec4 normal = position;
//                    normal.normalise();
//                    localVertices.push_back({ position, normal, mesh.col });
//                }
//            }
//
//            // Generate triangles.
//            // We can form triangles only for rows that have a row below (i.e. lat < totalRows - 1)
//            // Note: localVertices holds only the vertices for rows [startRow, endRow).
//            // Thus, only process rows from startRow up to endRow-1.
//            for (int lat = startRow; lat < endRow - 1; ++lat) {
//                // The local row index for 'lat' is (lat - startRow).
//                for (int lon = 0; lon < longitudeDivisions; ++lon) {
//                    int currentRow = lat - startRow;      // local index for current row
//                    int nextRow = currentRow + 1;           // local index for next row
//                    // Wrap-around: the next vertex in the row is at (lon+1) mod longitudeDivisions.
//                    int v0 = currentRow * longitudeDivisions + lon;
//                    int v1 = currentRow * longitudeDivisions + ((lon + 1) % longitudeDivisions);
//                    int v2 = nextRow * longitudeDivisions + lon;
//                    int v3 = nextRow * longitudeDivisions + ((lon + 1) % longitudeDivisions);
//
//                    localTriangles.push_back(triIndices(v0, v1, v2));
//                    localTriangles.push_back(triIndices(v1, v3, v2));
//                }
//            }
//
//            // Now merge the local vertices and triangles into the global mesh.
//            // IMPORTANT: Adjust the triangle indices by the current global vertex offset.
//            std::lock_guard<std::mutex> lock(meshMutex);
//            int vertexOffset = mesh.vertices.size();
//            mesh.vertices.insert(mesh.vertices.end(), localVertices.begin(), localVertices.end());
//            for (auto& tri : localTriangles) {
//                tri.v[0] += vertexOffset;
//                tri.v[1] += vertexOffset;
//                tri.v[2] += vertexOffset;
//            }
//            mesh.triangles.insert(mesh.triangles.end(), localTriangles.begin(), localTriangles.end());
//            });
//    }
//
//    for (auto& t : threads) {
//        t.join();
//    }
//    return mesh;
//}



};
