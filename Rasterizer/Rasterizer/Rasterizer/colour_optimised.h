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
            float r, g, b, a; // a is unused
            __m128 simd;
        };
        float rgb[4];     // Array representation of the RGB components
    };

public:
    // Enum for indexing the RGB components
    enum Colour { RED = 0, GREEN = 1, BLUE = 2 };


    colour(float _r = 0, float _g = 0, float _b = 0) : r(_r), g(_g), b(_b),a(1.f) {}

    void set(float _r, float _g, float _b) { simd = _mm_set_ps(1.0f, _b, _g, _r); }


    float& operator[] (Colour c) { return rgb[c]; }


    colour operator = (colour c) {
        simd = c.simd;
        return *this;
    }

    // Clamps the RGB components of the colour to the range [0, 1].
    void clampColour() {
        __m128 one = _mm_set1_ps(1.0f);
        simd = _mm_min_ps(simd, one);

    }


    void toRGB(unsigned char& cr, unsigned char& cg, unsigned char& cb) {
        __m128 scaled = _mm_mul_ps(simd, _mm_set1_ps(255.0f));
        //! not familiar with this
        cr = static_cast<unsigned char>(std::floor(r * 255));
        cg = static_cast<unsigned char>(std::floor(g * 255));
        cb = static_cast<unsigned char>(std::floor(b * 255));
    }


    colour operator * (const float scalar) {
        colour c;
        c.simd = _mm_mul_ps(simd, _mm_set1_ps(scalar));
        return c;
    }


    colour operator * (const colour& col) {
        colour c;
        c.simd = _mm_mul_ps(simd, col.simd);
        return c;
    }


    colour operator + (const colour& _c) {
        colour c;
        c.simd = _mm_add_ps(simd, _c.simd);
        return c;
    }
};
