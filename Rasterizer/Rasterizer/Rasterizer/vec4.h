#pragma once

#include <iostream>

// The `vec4` class represents a 4D vector and provides operations such as scaling, addition, subtraction, 
// normalization, and vector products (dot and cross).
class vec4 {
public:
    union {
        struct {
            float x, y, z, w; // Components of the vector
        };
        float v[4];           // Array representation of the vector components
    };

public:

    vec4(float _x = 0.f, float _y = 0.f, float _z = 0.f, float _w = 1.f)
        : x(_x), y(_y), z(_z), w(_w) {}

    // Displays the components of the vector in a readable format.
    void display() {
        std::cout << x << '\t' << y << '\t' << z << '\t' << w << std::endl;
    }

    vec4 operator*(float scalar) const {
        return { x * scalar, y * scalar, z * scalar, w * scalar };
    }

    void divideW() {
        x /= w;
        y /= w;
        z /= w;
        w = 1.f;
    }


    float& operator[](const unsigned int index) {
        return v[index];
    }

    float operator[](const unsigned int index) const {
        return v[index];
    }

   
    vec4 operator-(const vec4& other) const {
        return vec4(x - other.x, y - other.y, z - other.z, 0.0f);
    }

    vec4 operator+(const vec4& other) const {
        return vec4(x + other.x, y + other.y, z + other.z, 0.0f);
    }

    static vec4 cross(const vec4& v1, const vec4& v2) {
        return vec4(
            v1.y * v2.z - v1.z * v2.y,
            v1.z * v2.x - v1.x * v2.z,
            v1.x * v2.y - v1.y * v2.x,
            0.0f // The W component is set to 0 for cross products
        );
    }

    static float dot(const vec4& v1, const vec4& v2) {
        return (v1.x * v2.x) + (v1.y * v2.y) + (v1.z * v2.z);
    }


    void normalise() {
        float length = std::sqrt(x * x + y * y + z * z);
        x /= length;
        y /= length;
        z /= length;
    }
};
