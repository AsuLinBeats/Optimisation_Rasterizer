#pragma once
#pragma once

#include <cmath>
#include <cstdint>
#include <immintrin.h>
// The `colour` class represents an RGB colour with floating-point precision.
// It provides various utilities for manipulating and converting colours.
class alignas(16) colour {
    union {
        struct {
            // add a so it can align with sse registor
            float r, g, b, a; // Red, Green, and Blue components of the colour
            __m128 simd;
        };
        float rgb[4];     // Array representation of the RGB components
    };

public:
    // Enum for indexing the RGB components
    enum Colour { RED = 0, GREEN = 1, BLUE = 2 };

    // Constructor to initialize the colour with specified RGB values.
    // Default values are 0 (black).
    // Input Variables:
    // - _r: Red component (default 0.0f)
    // - _g: Green component (default 0.0f)
    // - _b: Blue component (default 0.0f)
    colour(float _r = 0, float _g = 0, float _b = 0) : r(_r), g(_g), b(_b) {}

    // Sets the RGB components of the colour.
    // Input Variables:
    // - _r: Red component
    // - _g: Green component
    // - _b: Blue component
    void set(float _r, float _g, float _b) { simd = _mm_set_ps(0.0f, _b, _g, _r); }

    // Accesses the specified component of the colour by index.
    // Input Variables:
    // - c: Index of the component (RED, GREEN, or BLUE)
    // Returns a reference to the specified component.
    float& operator[] (Colour c) { return rgb[c]; }

    // Assigns the values of another colour to this one.
    // Input Variables:
    // - c: The source color
    colour operator = (colour c) {
        simd = c.simd;
        return *this;
    }

    // Clamps the RGB components of the colour to the range [0, 1].
    void clampColour() {
        __m128 zero = _mm_setzero_ps();
        __m128 one = _mm_set1_ps(1.0f);
        simd = _mm_min_ps(_mm_max_ps(simd, zero), one);

    }

    // Converts the floating-point RGB values to integer values (0-255).
    // Output Variables:
    // - cr: Red component as an unsigned char
    // - cg: Green component as an unsigned char
    // - cb: Blue component as an unsigned char
    void toRGB(unsigned char& cr, unsigned char& cg, unsigned char& cb) {
        __m128 scaled = _mm_mul_ps(simd, _mm_set1_ps(255.0f));
        //! not familiar with this
        cr = static_cast<unsigned char>(std::floor(r * 255));
        cg = static_cast<unsigned char>(std::floor(g * 255));
        cb = static_cast<unsigned char>(std::floor(b * 255));
    }

    // Scales the RGB components of the colour by a scalar value.
    // Input Variables:
    // - scalar: The scaling factor
    // Returns a new `colour` object with scaled components.
    colour operator * (const float scalar) {
        colour c;
        c.simd = _mm_mul_ps(simd, _mm_set1_ps(scalar));
        return c;
    }

    // Multiplies the RGB components of this colour with another colour.
    // Input Variables:
    // - col: The other color to multiply with
    // Returns a new `colour` object with multiplied components.
    colour operator * (const colour& col) {
        colour c;
        c.simd = _mm_mul_ps(simd, col.simd);
        return c;
    }

    // Adds the RGB components of another colour to this one.
    // Input Variables:
    // - _c: The other colour to add
    // Returns a new `colour` object with added components.
    colour operator + (const colour& _c) {
        colour c;
        c.simd = _mm_add_ps(simd, _c.simd);
        return c;
    }
};
