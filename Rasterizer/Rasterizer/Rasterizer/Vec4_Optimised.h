#pragma once

#include <iostream>
#include <immintrin.h>
#include"TestTools.h"
// The `vec4O` class represents a 4D vector and provides operations such as scaling, addition, subtraction, 
// normalization, and vector products (dot and cross).
class vec4O {
    union {
        struct {
            float x, y, z, w; // Components of the vector
        };
        alignas(32) float v[4];           // Aligned Array representation of the vector components
    };

public:
    // Constructor to initialize the vector with specified values.
    // Default values: x = 0, y = 0, z = 0, w = 1.
    // Input Variables:
    // - _x: X component of the vector
    // - _y: Y component of the vector
    // - _z: Z component of the vector
    // - _w: W component of the vector (default is 1.0)
    vec4O(float _x = 0.f, float _y = 0.f, float _z = 0.f, float _w = 1.f)
        : x(_x), y(_y), z(_z), w(_w) {
    }

    // Displays the components of the vector in a readable format.
    void display() {
        std::cout << x << '\t' << y << '\t' << z << '\t' << w << std::endl;
    }

    //void init() {
    //    alignas(32) float v1[8];
    //    alignas(32) float v2[8];
    //    float result[8];
    //}

    vec4O add_SSE(vec4O &other) {
        vec4O result;
        // check if aligned to 32 bits
        if (reinterpret_cast<uintptr_t>(this->v) % 32 == 0 && reinterpret_cast<uintptr_t>(other.v) % 32 == 0) {
            std::cout << "aligned" << "\n";
            // Load 8 floats into 128-bit registers
            __m128 va = _mm_load_ps(this->v); // use array representation
            __m128 vb = _mm_load_ps(other.v);
            // Perform element-wise addition
            __m128 vc = _mm_add_ps(va, vb);
            // Store the result back to memory
            _mm_store_ps(result.v, vc);
            return result;
        }
    }

    vec4O sub_SSE(vec4O& other) {
        vec4O result;
        // check if aligned to 32 bits
        if (reinterpret_cast<uintptr_t>(this->v) % 16 == 0 && reinterpret_cast<uintptr_t>(other.v) % 32 == 0) {
            std::cout << "aligned" << "\n";
            // Load 8 floats into 128-bit registers
            __m128 va = _mm_load_ps(this->v); // use array representation
            __m128 vb = _mm_load_ps(other.v);
            // Perform element-wise addition
            __m128 vc = _mm_sub_ps(va, vb);
            // Store the result back to memory
            _mm_store_ps(result.v, vc);
            return result;
        }
    }

    vec4O mul_SSE(vec4O& other) {
        // multiplication between vectors
        vec4O result;
        // check if aligned to 16 bits
        if (reinterpret_cast<uintptr_t>(this->v) % 16 == 0 && reinterpret_cast<uintptr_t>(other.v) % 32 == 0) {
            std::cout << "aligned" << "\n";
            // Load 8 floats into 128-bit registers
            __m128 va = _mm_load_ps(this->v); // use array representation
            __m128 vb = _mm_load_ps(other.v);
            // Perform element-wise addition
            __m128 vc = _mm_mul_ps(va, vb);
            // Store the result back to memory
            _mm_store_ps(result.v, vc);
            return result;
        }
    }

    vec4O div_SSE(vec4O& other) {
        // division between vectors
        vec4O result;
        // check if aligned to 32 bits
        if (reinterpret_cast<uintptr_t>(this->v) % 16 == 0 && reinterpret_cast<uintptr_t>(other.v) % 32 == 0) {
            std::cout << "aligned" << "\n";
            // Load 8 floats into 128-bit registers
            __m128 va = _mm_load_ps(this->v); // use array representation
            __m128 vb = _mm_load_ps(other.v);
            // Perform element-wise addition
            __m128 vc = _mm_div_ps(va, vb);
            // Store the result back to memory
            _mm_store_ps(result.v, vc);
            return result;
        }
    }

    float dot_SSE(vec4O &other) {
        vec4O result;
        __m128 va = _mm_load_ps(this->v);
        __m128 vb = _mm_load_ps(other.v);
        __m128 vc = _mm_mul_ps(va, vb);
        _mm_store_ps(result.v, vc);
        // combine results
        float sum;
        for (int i = 0; i < 4; ++i){
            // because we are dealing with 4 floats, 
            sum += result[i];
        }
        return sum;
    }

    static vec4O cross_SSE(const vec4O& v1, const vec4O& v2) {
        vec4O result;
        __m128 va = _mm_set_ps(v1.z,v1.y,v1.x,0.f);
        __m128 vb = _mm_set_ps(v2.x,v2.z,v2.y,0.f);

        __m128 vc = _mm_setr_ps(v1.z, v1.x, v1.y, 0.f);
        __m128 vd = _mm_setr_ps(v2.y, v2.z, v2.x, 0.f);
        __m128 mul1 = _mm_mul_ps(va, vb); 
        __m128 mul2 = _mm_mul_ps(vc, vd);
        __m128 temp = _mm_sub_ps(mul1, mul2);

        _mm_storeu_ps(result.v, temp);
        vec4O final(result[0], result[1], result[2], 0.f);
        return final;
    }
    // Scales the vector by a scalar value.
    // Input Variables:
    // - scalar: Value to scale the vector by
    // Returns a new scaled `vec4O`.
    vec4O operator*(float scalar) const {
        return { x * scalar, y * scalar, z * scalar, w * scalar };
    }

    // Divides the vector by its W component and sets W to 1.
    // Useful for normalizing the W component after transformations.
    void divideW() {
        x /= w;
        y /= w;
        z /= w;
        w = 1.f;
    }

    // Accesses a vector component by index.
    // Input Variables:
    // - index: Index of the component (0 for x, 1 for y, 2 for z, 3 for w)
    // Returns a reference to the specified component.
    float& operator[](const unsigned int index) {
        return v[index];
    }

    // Accesses a vector component by index (const version).
    // Input Variables:
    // - index: Index of the component (0 for x, 1 for y, 2 for z, 3 for w)
    // Returns the specified component value.
    float operator[](const unsigned int index) const {
        return v[index];
    }

    // Subtracts another vector from this vector.
    // Input Variables:
    // - other: The vector to subtract
    // Returns a new `vec4O` resulting from the subtraction.
    vec4O operator-(vec4O& other) {
        return sub_SSE(other);
    }

    // Adds another vector to this vector.
    // Input Variables:
    // - other: The vector to add
    // Returns a new `vec4O` resulting from the addition.
    vec4O operator+(vec4O& other) {
        return add_SSE(other);
    }

    // Computes the cross product of two vectors.
    // Input Variables:
    // - v1: The first vector
    // - v2: The second vector
    // Returns a new `vec4O` representing the cross product.
    static vec4O cross(const vec4O& v1, const vec4O& v2) {
        return vec4O(
            v1.y * v2.z - v1.z * v2.y,
            v1.z * v2.x - v1.x * v2.z,
            v1.x * v2.y - v1.y * v2.x,
            0.0f // The W component is set to 0 for cross products
        );
    }

    // Computes the dot product of two vectors.
    // Input Variables:
    // - v1: The first vector
    // - v2: The second vector
    // Returns the dot product as a float.
    static float dot(const vec4O& v1, const vec4O& v2) {
        return (v1.x * v2.x) + (v1.y * v2.y) + (v1.z * v2.z);
    }

    // Normalizes the vector to make its length equal to 1.
    // This operation does not affect the W component.
    void normalise() {
        float length = std::sqrt(x * x + y * y + z * z);
        x /= length;
        y /= length;
        z /= length;
    }
};

void testAVX2() {
    vec4O a(1, 2, 3, 4);
    vec4O b(5,6,7,8);
    a.add_SSE(b);

}

void normalTest() {
    vec4O a(1, 2, 3, 4);
    vec4O b(5, 6, 7, 8);
    a = a + b;
}