#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>

#include "GamesEngineeringBase.h" // Include the GamesEngineeringBase header
#include <algorithm>
#include <chrono>

#include <cmath>
//#include "matrix.h"
#include "matrix_optimised.h"
//#include "colour.h"
#include "colour_optimised.h"
// #include "mesh.h"
#include "mesh_optimised.h"
// #include "zbuffer.h"
#include "zbuffer_optimised.h"
// #include "renderer.h"
#include "renderer_optimised.h"
#include "RNG.h"
#include "light.h"
// #include "triangle.h"
#include "triangle_optimised.h"

#include "TestTools.h" // introduce support functions to evaluate performance

void render(Renderer& renderer, Mesh* mesh, const matrix& camera, const Light& L) {
    // combined transformation matrix
    matrix p = renderer.perspective * camera * mesh->world;

    // Precompute transformed vertices
    std::vector<Vertex> transformedVertices(mesh->vertices.size());
    for (size_t i = 0; i < mesh->vertices.size(); ++i) {
        // Transform the vertex position and perform perspective division
        transformedVertices[i].p = p * mesh->vertices[i].p;
        transformedVertices[i].p.divideW();
        // Transform the normal into world space and normalize it (no perspective division needed)
        transformedVertices[i].normal = mesh->world * mesh->vertices[i].normal;
        transformedVertices[i].normal.normalise();
        transformedVertices[i].rgb = mesh->vertices[i].rgb;
    }

    // Get screen dimensions
    int screenWidth = renderer.canvas.getWidth();
    int screenHeight = renderer.canvas.getHeight();

    // Process each triangle in the mesh
    for (triIndices& ind : mesh->triangles) {
        // Retrieve the three vertices for the current triangle
        Vertex t[3] = {
            transformedVertices[ind.v[0]],
            transformedVertices[ind.v[1]],
            transformedVertices[ind.v[2]]
        };

        // Backface culling:
        //  - Compute a camera direction based on the first vertex position.
        //  - Compute two edge vectors and then the face normal.
        //  - If the triangle normal faces away from the camera, skip rendering.
        vec4 camDir = t[0].p;
        camDir.normalise();

        vec4 edge1 = t[1].p - t[0].p;
        vec4 edge2 = t[2].p - t[0].p;
        vec4 triangleNormal = vec4::cross(edge1, edge2);
        triangleNormal.normalise();

        if (triangleNormal.dot(camDir) > 0.0f)
            continue;

        // Map vertices to screen space
        for (int i = 0; i < 3; ++i) {
            t[i].p[0] = (t[i].p[0] + 1.f) * 0.5f * static_cast<float>(screenWidth);
            t[i].p[1] = (t[i].p[1] + 1.f) * 0.5f * static_cast<float>(screenHeight);
            t[i].p[1] = screenHeight - t[i].p[1]; // Invert y-axis for screen space
        }

        // Skip triangles off-screen
        bool allLeft = (t[0].p[0] < 0 && t[1].p[0] < 0 && t[2].p[0] < 0);
        bool allRight = (t[0].p[0] > screenWidth && t[1].p[0] > screenWidth && t[2].p[0] > screenWidth);
        bool allTop = (t[0].p[1] < 0 && t[1].p[1] < 0 && t[2].p[1] < 0);
        bool allBottom = (t[0].p[1] > screenHeight && t[1].p[1] > screenHeight && t[2].p[1] > screenHeight);
        if (allLeft || allRight || allTop || allBottom)
            continue;

        // Clip triangles outside the valid range [-1, 1]
        if (fabs(t[0].p[2]) > 1.0f || fabs(t[1].p[2]) > 1.0f || fabs(t[2].p[2]) > 1.0f)
            continue;

        // Draw triangles
        triangle tri(t[0], t[1], t[2]);
        tri.draw(renderer, L, mesh->ka, mesh->kd);
    }
}

void renderThread(Renderer& renderer, const std::vector<Mesh*>& meshes, const matrix& camera, const Light& L) {
    for (Mesh* mesh : meshes) {
        render(renderer, mesh, camera, L);
    }
}

void multithreadedRender(Renderer& renderer, const std::vector<Mesh*>& scene, const matrix& camera, const Light& L, int numThreads) {
    std::vector<std::thread> threads;
    size_t totalMeshes = scene.size();
    size_t meshesPerThread = (totalMeshes + numThreads - 1) / numThreads;  // Allocate meshes for each threads

    for (int i = 0; i < numThreads; i++) {
        size_t startIdx = i * meshesPerThread; // This is beginning of mesh subset of current thread
        size_t endIdx = min(startIdx + meshesPerThread, totalMeshes); // This make sure we will not beyond end of mesh
        // Check for validation
        if (startIdx >= endIdx)
            break; // No mesh left, exit
        std::vector<Mesh*> threadMeshes(scene.begin() + startIdx, scene.begin() + endIdx); // Meshes that current threads will process
        threads.emplace_back(renderThread, std::ref(renderer), threadMeshes, std::cref(camera), std::cref(L)); // spawn threads
    }
    for (auto& t : threads) {
        t.join();
    }
}

// Test scene function to demonstrate rendering with user-controlled transformations
// No input variables
void sceneTest() {
    Renderer renderer;
    // create light source {direction, diffuse intensity, ambient intensity}
    Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.1f, 0.1f, 0.1f) };
    // camera is just a matrix
    matrix camera = matrix::makeIdentity(); // Initialize the camera with identity matrix

    bool running = true; // Main loop control variable

    std::vector<Mesh*> scene; // Vector to store scene objects

    // Create a sphere and a rectangle mesh
    Mesh mesh = Mesh::makeSphere(1.0f, 10, 20);
    //Mesh mesh2 = Mesh::makeRectangle(-2, -1, 2, 1);

    // add meshes to scene
    scene.push_back(&mesh);
   // scene.push_back(&mesh2); 

    float x = 0.0f, y = 0.0f, z = -4.0f; // Initial translation parameters
    mesh.world = matrix::makeTranslation(x, y, z);
    //mesh2.world = matrix::makeTranslation(x, y, z) * matrix::makeRotateX(0.01f);
    int numThreads = 4;
    // Main rendering loop
    while (running) {
        renderer.canvas.checkInput(); // Handle user input
        renderer.clear(); // Clear the canvas for the next frame

        // Apply transformations to the meshes
     //   mesh2.world = matrix::makeTranslation(x, y, z) * matrix::makeRotateX(0.01f);
        mesh.world = matrix::makeTranslation(x, y, z);

        // Handle user inputs for transformations
        if (renderer.canvas.keyPressed(VK_ESCAPE)) break;
        if (renderer.canvas.keyPressed('A')) x += -0.1f;
        if (renderer.canvas.keyPressed('D')) x += 0.1f;
        if (renderer.canvas.keyPressed('W')) y += 0.1f;
        if (renderer.canvas.keyPressed('S')) y += -0.1f;
        if (renderer.canvas.keyPressed('Q')) z += 0.1f;
        if (renderer.canvas.keyPressed('E')) z += -0.1f;

        multithreadedRender(renderer, scene, camera, L, numThreads);

        renderer.present(); // Display the rendered frame
    }
}

// Utility function to generate a random rotation matrix
// No input variables
matrix makeRandomRotation() {
    RandomNumberGenerator& rng = RandomNumberGenerator::getInstance();
    unsigned int r = rng.getRandomInt(0, 3);

    switch (r) {
    case 0: return matrix::makeRotateX(rng.getRandomFloat(0.f, 2.0f * M_PI));
    case 1: return matrix::makeRotateY(rng.getRandomFloat(0.f, 2.0f * M_PI));
    case 2: return matrix::makeRotateZ(rng.getRandomFloat(0.f, 2.0f * M_PI));
    default: return matrix::makeIdentity();
    }
}

// Function to render a scene with multiple objects and dynamic transformations
// No input variables
void scene1() {
    Renderer renderer;
    matrix camera;
    Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.1f, 0.1f, 0.1f) };

    bool running = true;

    std::vector<Mesh*> scene;

    // Create a scene of 40 cubes with random rotations
    for (unsigned int i = 0; i < 20; i++) {
        Mesh* m = new Mesh();
        *m = Mesh::makeCube(1.f);
        m->world = matrix::makeTranslation(-2.0f, 0.0f, (-3 * static_cast<float>(i))) * makeRandomRotation();
        scene.push_back(m);
        m = new Mesh();
        *m = Mesh::makeCube(1.f);
        m->world = matrix::makeTranslation(2.0f, 0.0f, (-3 * static_cast<float>(i))) * makeRandomRotation();
        scene.push_back(m);
    }

    float zoffset = 8.0f; // Initial camera Z-offset
    float step = -0.1f;  // Step size for camera movement

    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    int cycle = 0;

    int numThreads = 20;
    // Main rendering loop
    while (running) {
        renderer.canvas.checkInput();
        renderer.clear();

        camera = matrix::makeTranslation(0, 0, -zoffset); // Update camera position

        // Rotate the first two cubes in the scene
        scene[0]->world = scene[0]->world * matrix::makeRotateXYZ(0.1f, 0.1f, 0.0f);
        scene[1]->world = scene[1]->world * matrix::makeRotateXYZ(0.0f, 0.1f, 0.2f);

        if (renderer.canvas.keyPressed(VK_ESCAPE)) break;

        zoffset += step;
        if (zoffset < -60.f || zoffset > 8.f) {
            step *= -1.f;
            if (++cycle % 2 == 0) {
                end = std::chrono::high_resolution_clock::now();
                std::cout << cycle / 2 << " :" << std::chrono::duration<double, std::milli>(end - start).count() << "ms\n";
                start = std::chrono::high_resolution_clock::now();
            }
        }

       multithreadedRender(renderer, scene, camera, L, numThreads);
        //for (Mesh* m : scene) {
        //    render(renderer, m, camera, L);
        //}
        renderer.present();
    }

    for (auto& m : scene)
        delete m;
}

// Scene with a grid of cubes and a moving sphere
// No input variables
void scene2() {
    Renderer renderer;
    matrix camera = matrix::makeIdentity();
    Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.1f, 0.1f, 0.1f) };

    std::vector<Mesh*> scene;

    struct rRot { float x; float y; float z; }; // Structure to store random rotation parameters
    std::vector<rRot> rotations;

    RandomNumberGenerator& rng = RandomNumberGenerator::getInstance();

    // Create a grid of cubes with random rotations
    for (unsigned int y = 0; y < 6; y++) {
        for (unsigned int x = 0; x < 8; x++) {
            Mesh* m = new Mesh();
            *m = Mesh::makeCube(1.f);
            scene.push_back(m);
            m->world = matrix::makeTranslation(-7.0f + (static_cast<float>(x) * 2.f), 5.0f - (static_cast<float>(y) * 2.f), -8.f);
            rRot r{ rng.getRandomFloat(-.1f, .1f), rng.getRandomFloat(-.1f, .1f), rng.getRandomFloat(-.1f, .1f) };
            rotations.push_back(r);
        }
    }

    // Create a sphere and add it to the scene
    Mesh* sphere = new Mesh();
    *sphere = Mesh::makeSphere(1.0f, 10, 20);
    scene.push_back(sphere);
    float sphereOffset = -6.f;
    float sphereStep = 0.1f;
    sphere->world = matrix::makeTranslation(sphereOffset, 0.f, -6.f);

    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    int cycle = 0;
    int numThreads = 20;
    bool running = true;
    while (running) {
        renderer.canvas.checkInput();
        renderer.clear();

        // Rotate each cube in the grid
        for (unsigned int i = 0; i < rotations.size(); i++)
            scene[i]->world = scene[i]->world * matrix::makeRotateXYZ(rotations[i].x, rotations[i].y, rotations[i].z);

        // Move the sphere back and forth
        sphereOffset += sphereStep;
        sphere->world = matrix::makeTranslation(sphereOffset, 0.f, -6.f);
        if (sphereOffset > 6.0f || sphereOffset < -6.0f) {
            sphereStep *= -1.f;
            if (++cycle % 2 == 0) {
                end = std::chrono::high_resolution_clock::now();
                std::cout << cycle / 2 << " :" << std::chrono::duration<double, std::milli>(end - start).count() << "ms\n";
                start = std::chrono::high_resolution_clock::now();
            }
        }

        if (renderer.canvas.keyPressed(VK_ESCAPE)) break;

        //multithreadedRender(renderer, scene, camera, L, numThreads);
		for (Mesh* m : scene) {
			render(renderer, m, camera, L);
		}
        renderer.present();
    }

    for (auto& m : scene)
        delete m;
}

void scene3() {
    Renderer renderer;
    matrix camera;
    Light L{ vec4(0.f, 5.f, 5.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.1f, 0.1f, 0.1f) };

    std::vector<Mesh*> scene;

    // 10000 cube
    const int GRID_SIZE = 100;
    const float SPACING = 2.0f;

    
    for (int y = 0; y < GRID_SIZE; ++y) {
        for (int x = 0; x < GRID_SIZE; ++x) {
            Mesh* m = new Mesh();
            *m = Mesh::makeCube(0.5f);  
            m->world = matrix::makeTranslation(
                (x - GRID_SIZE / 2) * SPACING,
                (y - GRID_SIZE / 2) * SPACING,
                -50.f  
            );
            scene.push_back(m);
        }
    }

    // camera moving is similar with scene1
    float zoffset = 50.f;
    float step = -0.2f;
    bool running = true;
    auto start = std::chrono::high_resolution_clock::now();
    int cycle = 0;
    int numThreads = 20;
    while (running) {
        renderer.canvas.checkInput();
        renderer.clear();

   
        camera = matrix::makeTranslation(0, 0, -zoffset);
        zoffset += step;

      
        if (zoffset < 10.f || zoffset > 50.f) {
            step *= -1.f;
            if (++cycle % 2 == 0) {
                auto end = std::chrono::high_resolution_clock::now();
                std::cout << "Cycle " << cycle / 2 << ": "
                    << std::chrono::duration<double, std::milli>(end - start).count()
                    << "ms" << std::endl;
                start = end;
            }
        }

        //multithreadedRender(renderer, scene, camera, L, numThreads);
        for (Mesh* m : scene) {
            render(renderer, m, camera, L);
        }
        renderer.present();
    }

    for (auto& m : scene) delete m;
}

void printFPS() {
    static std::chrono::time_point<std::chrono::steady_clock> oldTime = std::chrono::high_resolution_clock::now();
    static int fps; fps++;

    if (std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - oldTime) >= std::chrono::seconds{ 1 }) {
        oldTime = std::chrono::high_resolution_clock::now();
        std::cout << "FPS: " << fps << std::endl;
        fps = 0;
    }
}

// Entry point of the application
// No input variables
int main() {
    // Time test
    auto start = std::chrono::high_resolution_clock::now();
    scene1();
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Time doing some work : "
        << std::chrono::duration<double, std::milli>(end - start).count()
        << " ms\n";


    //scene3();
    //printFPS();
   // std::cout << sum << std::endl;
    // Uncomment the desired scene function to run
    // scene1();
   // scene2();
    //sceneTest(); 
    

    return 0;
}