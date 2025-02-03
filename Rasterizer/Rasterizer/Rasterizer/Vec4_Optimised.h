#pragma once

#include <iostream>
#include <immintrin.h>
#include"TestTools.h"

class alignas(16) vec4 {
public:
    union {
        struct {
            float x, y, z, w; // Components of the vector
        };
       float v[4];           // Aligned directly when define-----can remove alignment check later
       __m128 simd;       
    };

public:

    vec4(float _x = 0.f, float _y = 0.f, float _z = 0.f, float _w = 1.f)
        {
        simd = _mm_setr_ps(_x, _y, _z, _w);
    }
    
    vec4(__m128 m) {
        simd = m;
    }

    vec4(__m128 m, float w) {
        simd = m;
        //// case when w need to be restore
        float arr[4];
        arr[3] = w;
        simd = _mm_load_ps(arr);
    }


    // Displays the components of the vector in a readable format.
    void display() {
        std::cout << x << '\t' << y << '\t' << z << '\t' << w << std::endl;
    }

    vec4 mul_SSE(const vec4& other) const {
        // multiplication between vectors
        vec4 result;
        __m128 va = _mm_load_ps(this->v); // use array representation
        __m128 vb = _mm_load_ps(other.v);
        // Perform element-wise addition
        __m128 vc = _mm_mul_ps(va, vb);
        // Store the result back to memory
        _mm_store_ps(result.v, vc);
        return result;
    }

    vec4 div_SSE(const vec4& other) const {
        // division between vectors
        vec4 result;

        __m128 va = _mm_load_ps(this->v); // use array representation
        __m128 vb = _mm_load_ps(other.v);
        // Perform element-wise addition
        __m128 vc = _mm_div_ps(va, vb);
        // Store the result back to memory
        _mm_store_ps(result.v, vc);
        return result;
        //}
    }


    static float dot(const vec4& a, const vec4& b) {
        // better version using horizontal add(so do not need to use temp array result
        __m128 va = _mm_load_ps(a.v);
        __m128 vb = _mm_load_ps(b.v);
        __m128 vc = _mm_mul_ps(va, vb);

        vc = _mm_hadd_ps(vc, vc);
        vc = _mm_hadd_ps(vc, vc);
        return _mm_cvtss_f32(vc);
    }

    float dot(const vec4& other) const {
        __m128 va = _mm_load_ps(v);           
        __m128 vb = _mm_load_ps(other.v);        
        __m128 prod = _mm_mul_ps(va, vb);

        prod = _mm_hadd_ps(prod, prod);
        prod = _mm_hadd_ps(prod, prod);
        return _mm_cvtss_f32(prod);
    }
    

    static vec4 cross(const vec4& v1, const vec4& v2) {
        // va = [v1.z, v1.x, v1.y, 0]
        __m128 va = _mm_shuffle_ps(_mm_load_ps(v1.v), _mm_load_ps(v1.v), _MM_SHUFFLE(3, 0, 2, 1));
        // vb = [v2.y, v2.z, v2.x, 0]
        __m128 vb = _mm_shuffle_ps(_mm_load_ps(v2.v), _mm_load_ps(v2.v), _MM_SHUFFLE(3, 1, 0, 2));
        __m128 vc = _mm_mul_ps(va, vb);
        // va_new = [v1.y, v1.z, v1.x, 0]
        va = _mm_shuffle_ps(va, va, _MM_SHUFFLE(3, 0, 2, 1));
        // vb_new = [v2.z, v2.x, v2.y, 0]
        vb = _mm_shuffle_ps(vb, vb, _MM_SHUFFLE(3, 1, 0, 2));
        __m128 vd = _mm_mul_ps(va, vb);

        return vec4(_mm_sub_ps(vc, vd));
    }

    vec4 operator*(const vec4& other) const {
        return mul_SSE(other);
    }

    vec4 operator*(float scalar) const {
        __m128 s = _mm_set1_ps(scalar);
        return vec4(_mm_mul_ps(simd, s));
    }

    void divideW() {
        float temp = 1.f / w;
        __m128 tempa = _mm_setr_ps(temp, temp, temp, 1.f);
        // Multiply x, y, z by recip. (The fourth component will be multiplied by 1.)
        simd = _mm_mul_ps(simd, tempa);
        // set w separately
        w = 1.f;
    }

    float& operator[](const unsigned int index) {
        return v[index];
    }


    float operator[](const unsigned int index) const {
        return v[index];
    }


    vec4 operator-(const vec4& other) const {
        __m128 dif = _mm_sub_ps(simd, other.simd);
        return vec4(dif);
    }


    vec4 operator+(const vec4& other) const {
        __m128 sum = _mm_add_ps(simd, other.simd);
        return vec4(sum);
    }

    vec4 operator/(const vec4& other) const {
        return div_SSE(other);
    }


    void normalise() {
        // Compute squared length of the (x,y,z) part.
        float lenSq = dot(*this, *this);

        float invLen = 1.f / std::sqrt(lenSq);
        __m128 inv = _mm_set1_ps(invLen);
  
        __m128 norm = _mm_mul_ps(simd, inv);
        //  restore w
        float arr[4];
        _mm_store_ps(arr, norm);
        arr[3] = w;
        simd = _mm_load_ps(arr);
    }

};
