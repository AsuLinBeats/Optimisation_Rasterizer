#pragma once

// #include "mesh.h"
#include "mesh_optimised.h"
//#include "colour.h"
#include "colour_optimised.h"
// #include "renderer.h"
#include "renderer_optimised.h"
#include "light.h"
#include <iostream>
#include<thread>
#include<mutex>


// Simple support class for a 2D vector
class vec2D {
public:
    float x, y;

    // Default constructor initializes both components to 0
    vec2D() { x = y = 0.f; };

    // Constructor initializes components with given values
    vec2D(float _x, float _y) : x(_x), y(_y) {}

    // Constructor initializes components from a vec4
    vec2D(vec4 v) {
        x = v[0];
        y = v[1];
    }

    // Display the vector components
    void display() { std::cout << x << '\t' << y << std::endl; }

    // Overloaded subtraction operator for vector subtraction
    vec2D operator- (vec2D& v) {
        vec2D q;
        q.x = x - v.x;
        q.y = y - v.y;
        return q;
    }
};

// Class representing a triangle for rendering purposes
class triangle {
    Vertex v[3];       // Vertices of the triangle
    float area;        // Area of the triangle
    colour col[3];     // Colors for each vertex of the triangle

public:

    triangle(const Vertex& v1, const Vertex& v2, const Vertex& v3) : v{v1, v2, v3} {
        // Calculate the 2D area of the triangle
        const vec2D e1 = vec2D(v[1].p - v[0].p);
        const vec2D e2 = vec2D(v[2].p - v[0].p);
        area = std::abs(e1.x * e2.y - e1.y * e2.x);
    }

    float getC(vec2D v1, vec2D v2, vec2D p) {
        vec2D e = v2 - v1;
        vec2D q = p - v1;
        return q.y * e.x - q.x * e.y;
    }


    bool getCoordinates(vec2D p, float& alpha, float& beta, float& gamma) {
        alpha = getC(vec2D(v[0].p), vec2D(v[1].p), p) / area;
        beta = getC(vec2D(v[1].p), vec2D(v[2].p), p) / area;
        gamma = getC(vec2D(v[2].p), vec2D(v[0].p), p) / area;

        if (alpha < 0.f || beta < 0.f || gamma < 0.f) return false;
        return true;
    }


    template <typename T>
    T interpolate(float alpha, float beta, float gamma, T a1, T a2, T a3) {
        return (a1 * alpha) + (a2 * beta) + (a3 * gamma);
    }

    void draw(Renderer& renderer, const Light& L, float ka, float kd)
    {

        // exit if the triangle area is negligible
        if (area < 1.f) return;

        // Determine the bounding box for this triangle
        vec2D minV, maxV;
        getBoundsWindow(renderer.canvas, minV, maxV);

        int minX = static_cast<int>(minV.x);
        int maxX = static_cast<int>(ceil(maxV.x));
        int minY = static_cast<int>(minV.y);
        int maxY = static_cast<int>(ceil(maxV.y));


        // Prepare SSE/AVX variables for lighting
        // Normalize the light direction
        vec4 normalizedLight = L.omega_i;
        normalizedLight.normalise();

        // Convert light direction and colors to SIMD vectors
        __m128 lightDir = _mm_set_ps(0.0f, normalizedLight.z, normalizedLight.y, normalizedLight.x);
        __m128 lightColor = _mm_set_ps(0.0f, L.L.b, L.L.g, L.L.r);
        __m128 ambientColor = _mm_set_ps(0.0f, L.ambient.b, L.ambient.g, L.ambient.r);
        __m128 kdSimd = _mm_set1_ps(kd);


        // Cache vertex information
        const vec2D& v0 = v[0].p;
        const vec2D& v1 = v[1].p;
        const vec2D& v2 = v[2].p;

        // Precompute reciprocal of the triangle area
        float invArea = 1.0f / area;


        // Fill the pixels within the bounding box
        for (int y = minY; y < maxY; y++)
        {
            for (int x = minX; x < maxX; x++)
            {

                // Compute barycentric coordinates
  
                vec2D p{ static_cast<float>(x), static_cast<float>(y) };
                float alpha = getC(v0, v1, p) * invArea;
                float beta = getC(v1, v2, p) * invArea;
                float gamma = getC(v2, v0, p) * invArea;

                // Check if the point is inside the triangle
                if (alpha >= 0.f && beta >= 0.f && gamma >= 0.f)
                {

                    // Interpolate color using SSE 
                    // Weights
                    __m128 w = _mm_set_ps(0.0f, alpha, beta, gamma);

                    // Triangle vertex colors
                    __m128 c1 = _mm_set_ps(0.0f, v[0].rgb.r, v[0].rgb.g, v[0].rgb.b);
                    __m128 c2 = _mm_set_ps(0.0f, v[1].rgb.r, v[1].rgb.g, v[1].rgb.b);
                    __m128 c3 = _mm_set_ps(0.0f, v[2].rgb.r, v[2].rgb.g, v[2].rgb.b);

                    // Combine colors: result = c1*alpha + c2*beta + c3*gamma
                    __m128 part1 = _mm_mul_ps(c1, _mm_shuffle_ps(w, w, _MM_SHUFFLE(0, 0, 0, 0)));
                    __m128 part2 = _mm_mul_ps(c2, _mm_shuffle_ps(w, w, _MM_SHUFFLE(1, 1, 1, 1)));
                    __m128 part3 = _mm_mul_ps(c3, _mm_shuffle_ps(w, w, _MM_SHUFFLE(2, 2, 2, 2)));
                    __m128 result = _mm_add_ps(part1, _mm_add_ps(part2, part3));

                    // Convert SSE result to colour struct
                    colour c;
                    c.r = _mm_cvtss_f32(_mm_shuffle_ps(result, result, _MM_SHUFFLE(0, 0, 0, 0)));
                    c.g = _mm_cvtss_f32(_mm_shuffle_ps(result, result, _MM_SHUFFLE(1, 1, 1, 1)));
                    c.b = _mm_cvtss_f32(_mm_shuffle_ps(result, result, _MM_SHUFFLE(2, 2, 2, 2)));

                    c.clampColour();


                    //Interpolate depth and normal

                    float depth = interpolate(beta, gamma, alpha, v[0].p[2], v[1].p[2], v[2].p[2]);
                    vec4 normal = interpolate(beta, gamma, alpha, v[0].normal, v[1].normal, v[2].normal);
                    normal.normalise();

                    // Depth testing and final shading
                    if (depth > 0.01f && renderer.zbuffer(x, y) > depth)
                    {
                        // Convert normal to SIMD
                        __m128 normalSimd = _mm_set_ps(0.0f, normal.z, normal.y, normal.x);

                        // Dot product between lightDir and the interpolated normal
                        __m128 dotSimd = _mm_dp_ps(lightDir, normalSimd, 0x71);
                        float dot = _mm_cvtss_f32(dotSimd);

                        dot = max(dot, 0.0f);

                        // Calculate shaded color: (lightColor*dot + ambientColor) * kd * base_color
                        __m128 cSimd = _mm_set_ps(0.0f, c.b, c.g, c.r);
                        __m128 litColor = _mm_add_ps(_mm_mul_ps(lightColor, _mm_set1_ps(dot)), ambientColor);
                        __m128 shadedColor = _mm_mul_ps(_mm_mul_ps(kdSimd, litColor), cSimd);

                        // Convert SIMD color to unsigned char
                        unsigned char r = static_cast<unsigned char>(_mm_cvtss_f32(
                            _mm_shuffle_ps(shadedColor, shadedColor, _MM_SHUFFLE(0, 0, 0, 0))) * 255);
                        unsigned char g = static_cast<unsigned char>(_mm_cvtss_f32(
                            _mm_shuffle_ps(shadedColor, shadedColor, _MM_SHUFFLE(1, 1, 1, 1))) * 255);
                        unsigned char b = static_cast<unsigned char>(_mm_cvtss_f32(
                            _mm_shuffle_ps(shadedColor, shadedColor, _MM_SHUFFLE(2, 2, 2, 2))) * 255);

                        renderer.canvas.draw(x, y, r, g, b);
                        renderer.zbuffer(x, y) = depth;
                    }
                }
            }
        }
    }

    void getBounds(vec2D& minV, vec2D& maxV) {
        // Load the first two components (x and y) of each vertex into SSE registers.
        // We use _mm_setr_ps to set the lower two floats to (x,y) and the upper two to 0.
        __m128 a = _mm_setr_ps(v[0].p[0], v[0].p[1], 0.f, 0.f);
        __m128 b = _mm_setr_ps(v[1].p[0], v[1].p[1], 0.f, 0.f);
        __m128 c = _mm_setr_ps(v[2].p[0], v[2].p[1], 0.f, 0.f);

        // Compute the component‑wise minimum:
        __m128 min_ab = _mm_min_ps(a, b);
        __m128 min_abc = _mm_min_ps(min_ab, c);

        // Compute the component‑wise maximum:
        __m128 max_ab = _mm_max_ps(a, b);
        __m128 max_abc = _mm_max_ps(max_ab, c);

        // Extract the lower two floats (x and y) from each result.
        // (Since we only care about the first two elements, we can use a temporary array.)
        float bounds[4];
        _mm_storeu_ps(bounds, min_abc);
        minV.x = bounds[0];
        minV.y = bounds[1];

        _mm_storeu_ps(bounds, max_abc);
        maxV.x = bounds[0];
        maxV.y = bounds[1];
    }



    void getBoundsWindow(GamesEngineeringBase::Window& canvas, vec2D& minV, vec2D& maxV) {
        // First compute the raw bounds.
        getBounds(minV, maxV);

        // Load the computed min and max bounds into SSE registers.
        __m128 minVec = _mm_setr_ps(minV.x, minV.y, 0.f, 0.f);
        __m128 maxVec = _mm_setr_ps(maxV.x, maxV.y, 0.f, 0.f);

        // Create a zero vector for clamping the minimum.
        __m128 zero = _mm_setzero_ps();
        // Clamp minVec to be no less than zero.
        minVec = _mm_max_ps(minVec, zero);

        // Prepare a vector with the canvas dimensions (as floats).
        float canvasW = static_cast<float>(canvas.getWidth());
        float canvasH = static_cast<float>(canvas.getHeight());
        __m128 canvasMax = _mm_setr_ps(canvasW, canvasH, 0.f, 0.f);
        // Clamp maxVec to be no greater than the canvas dimensions.
        maxVec = _mm_min_ps(maxVec, canvasMax);

        // Store the clamped results back to minV and maxV.
        float bounds[4];
        _mm_storeu_ps(bounds, minVec);
        minV.x = bounds[0];
        minV.y = bounds[1];

        _mm_storeu_ps(bounds, maxVec);
        maxV.x = bounds[0];
        maxV.y = bounds[1];
    }

    void drawBounds(GamesEngineeringBase::Window& canvas) {
        vec2D minV, maxV;
        getBounds(minV, maxV);

        int xs = static_cast<int>(minV.x); // start from minimum x
        int xe = static_cast<int>(maxV.x); // end at max x
        int ys = static_cast<int>(minV.y); // start from minimum y
        int ye = static_cast<int>(maxV.y); // end at max y

        int rows = ye - ys;

        // threads set
        unsigned int numThreads = std::thread::hardware_concurrency();
        int rowsPerThread = (rows + numThreads - 1) / numThreads; // Ensure correct distribution of rows
        std::vector<std::thread> threads;
        threads.reserve(numThreads);

        for (unsigned int t = 0; t < numThreads; t++) {
            int sliceYStart = ys + t * rowsPerThread;
            int sliceYEnd = min(ye, sliceYStart + rowsPerThread); // Ensure we don't go out of bounds

            threads.emplace_back([sliceYStart, sliceYEnd, xs, xe, &canvas]() {
                for (int y = sliceYStart; y < sliceYEnd; y++) {
                    for (int x = xs; x < xe; x++) {
                        canvas.draw(x, y, 255, 0, 0);
                    }
                }
                });
        }

        // Wait for all threads to finish.
        for (auto& t : threads) {
            t.join();
        }
    }


    // Debugging utility to display the coordinates of the triangle vertices
    void display() {
        for (unsigned int i = 0; i < 3; i++) {
            v[i].p.display();
        }
        std::cout << std::endl;
    }
};


