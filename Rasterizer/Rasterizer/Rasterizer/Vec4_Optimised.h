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
        //float v[4];           
    };

public:

    vec4(float _x = 0.f, float _y = 0.f, float _z = 0.f, float _w = 1.f)
        : x(_x), y(_y), z(_z), w(_w) {
    }

    // Displays the components of the vector in a readable format.
    void display() {
        std::cout << x << '\t' << y << '\t' << z << '\t' << w << std::endl;
    }



    vec4 add_SSE(const vec4& other)const {
        vec4 result;

        __m128 va = _mm_load_ps(this->v); // use array representation
        __m128 vb = _mm_load_ps(other.v);
        // Perform element-wise addition
        __m128 vc = _mm_add_ps(va, vb);
        // Store the result back to memory
        _mm_store_ps(result.v, vc);
        return result;
        // }
    }

    vec4 sub_SSE(const vec4& other) const {
        vec4 result;

        __m128 va = _mm_load_ps(this->v); // use array representation
        __m128 vb = _mm_load_ps(other.v);
        // Perform element-wise addition
        __m128 vc = _mm_sub_ps(va, vb);
        // Store the result back to memory
        _mm_store_ps(result.v, vc);
        return result;
        //}
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
        // }
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

    //static float dot(vec4& first, vec4& other) {
    //    vec4 result;
    //    __m128 va = _mm_load_ps(first.v);
    //    __m128 vb = _mm_load_ps(other.v);
    //    __m128 vc = _mm_mul_ps(va, vb);
    //    _mm_store_ps(result.v, vc);
    //    // combine results
    //    float sum = 0;
    //    for (int i = 0; i < 4; ++i) {
    //        // because we are dealing with 4 floats, 
    //        sum += result[i];
    //    }
    //    return sum;
    //}

    static float dot(const vec4& a, const vec4& b) {
        // better version using horizontal add(so do not need to use temp array result
        __m128 va = _mm_load_ps(a.v);
        __m128 vb = _mm_load_ps(b.v);
        __m128 vc = _mm_mul_ps(va, vb);

        vc = _mm_hadd_ps(vc, vc);
        vc = _mm_hadd_ps(vc, vc);
        return _mm_cvtss_f32(vc);
    }

    //static vec4 cross(const vec4& v1, const vec4& v2) {
    //    vec4 result;
    //    __m128 va = _mm_set_ps(v1.z, v1.y, v1.x, 0.f);
    //    __m128 vb = _mm_set_ps(v2.x, v2.z, v2.y, 0.f);

    //    __m128 vc = _mm_setr_ps(v1.z, v1.x, v1.y, 0.f);
    //    __m128 vd = _mm_setr_ps(v2.y, v2.z, v2.x, 0.f);
    //    __m128 mul1 = _mm_mul_ps(va, vb);
    //    __m128 mul2 = _mm_mul_ps(vc, vd);
    //    __m128 temp = _mm_sub_ps(mul1, mul2);

    //    _mm_storeu_ps(result.v, temp);
    //    vec4 final(result[0], result[1], result[2], 0.f);
    //    return final;
    //}

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

        vec4 result;
        _mm_store_ps(result.v, _mm_sub_ps(vc, vd));
        result.w = 0.f;
        return result;
    }

    // Scales the vector by a scalar value.
    // Input Variables:
    // - scalar: Value to scale the vector by
    // Returns a new scaled `vec4`.
    vec4 operator*(const vec4& other) const {
        return mul_SSE(other);
    }

    vec4 operator*(float scalar) const {
        __m128 vs = _mm_set1_ps(scalar);
        __m128 va = _mm_load_ps(v);
        vec4 result;
        _mm_store_ps(result.v, _mm_mul_ps(va, vs));
        return result;
    }

    void divideW() {
       // __m128 vv = _mm_set_ps(1, z, y, x);
        __m128 vv = _mm_load_ps(this->v);
        __m128 vw = _mm_set_ps(w, w, w, w);
        _mm_store_ps(this->v, _mm_div_ps(vv, vw));
    }

    //void divideW() {
    //    x /= w;
    //    y /= w;
    //    z /= w;
    //    w = 1.f;
    //}

    float& operator[](const unsigned int index) {
        return v[index];
    }


    float operator[](const unsigned int index) const {
        return v[index];
    }


    vec4 operator-(const vec4& other) const {
        return sub_SSE(other);
    }


    vec4 operator+(const vec4& other) const {
        return add_SSE(other);
    }

    vec4 operator/(const vec4& other) const {
        return div_SSE(other);
    }

    //void normalise() {
    //    float length = std::sqrt(x * x + y * y + z * z);

    //    x /= length;
    //    y /= length;
    //    z /= length;
    //}

    void normalise() {
        // set w to 0
        __m128 va = _mm_set_ps(0.f, z, y, x);
        // z^2, y^2, x^2
        __m128 vsq = _mm_mul_ps(va, va);
        // z^2, 
        vsq = _mm_hadd_ps(vsq, vsq);
        vsq = _mm_hadd_ps(vsq, vsq);
        __m128 vlen = _mm_sqrt_ps(vsq);
        
        __m128 vnorm = _mm_div_ps(va, vlen);
        alignas(16) float tmp[4];
        _mm_store_ps(tmp, vnorm);
        x = tmp[0];
        y = tmp[1];
        z = tmp[2];
    }

};

//! initial version
//class alignas(16) vec4 {
//public:
//    union {
//        struct {
//            float x, y, z, w; // Components of the vector
//        };
//        float v[4];           // Aligned directly when define-----can remove alignment check later
//        //float v[4];           
//    };
//
//public:
//    vec4(float _x = 0.f, float _y = 0.f, float _z = 0.f, float _w = 1.f)
//        : x(_x), y(_y), z(_z), w(_w) {
//    }
//    // Displays the components of the vector in a readable format.
//    void display() {
//        std::cout << x << '\t' << y << '\t' << z << '\t' << w << std::endl;
//    }
//
//    vec4 add_SSE(const vec4& other)const {
//        vec4 result;
//
//        __m128 va = _mm_load_ps(this->v); // use array representation
//        __m128 vb = _mm_load_ps(other.v);
//        // Perform element-wise addition
//        __m128 vc = _mm_add_ps(va, vb);
//        // Store the result back to memory
//        _mm_store_ps(result.v, vc);
//        return result;
//        // }
//    }
//
//    vec4 sub_SSE(const vec4& other) const {
//        vec4 result;
//
//        __m128 va = _mm_load_ps(this->v); // use array representation
//        __m128 vb = _mm_load_ps(other.v);
//        // Perform element-wise addition
//        __m128 vc = _mm_sub_ps(va, vb);
//        // Store the result back to memory
//        _mm_store_ps(result.v, vc);
//        return result;
//        //}
//    }
//
//    vec4 mul_SSE(const vec4& other) const {
//        // multiplication between vectors
//        vec4 result;
//
//        __m128 va = _mm_load_ps(this->v); // use array representation
//        __m128 vb = _mm_load_ps(other.v);
//        // Perform element-wise addition
//        __m128 vc = _mm_mul_ps(va, vb);
//        // Store the result back to memory
//        _mm_store_ps(result.v, vc);
//        return result;
//        // }
//    }
//
//    vec4 div_SSE(const vec4& other) const {
//        // division between vectors
//        vec4 result;
//
//        __m128 va = _mm_load_ps(this->v); // use array representation
//        __m128 vb = _mm_load_ps(other.v);
//        // Perform element-wise addition
//        __m128 vc = _mm_div_ps(va, vb);
//        // Store the result back to memory
//        _mm_store_ps(result.v, vc);
//        return result;
//        //}
//    }
//
//    static float dot(vec4& first, vec4& other) {
//        vec4 result;
//        __m128 va = _mm_load_ps(first.v);
//        __m128 vb = _mm_load_ps(other.v);
//        __m128 vc = _mm_mul_ps(va, vb);
//        _mm_store_ps(result.v, vc);
//        // combine results
//        float sum = 0;
//        for (int i = 0; i < 4; ++i) {
//            // because we are dealing with 4 floats, 
//            sum += result[i];
//        }
//        return sum;
//    }
//
//    static vec4 cross(const vec4& v1, const vec4& v2) {
//        vec4 result;
//        __m128 va = _mm_set_ps(v1.z, v1.y, v1.x, 0.f);
//        __m128 vb = _mm_set_ps(v2.x, v2.z, v2.y, 0.f);
//
//        __m128 vc = _mm_setr_ps(v1.z, v1.x, v1.y, 0.f);
//        __m128 vd = _mm_setr_ps(v2.y, v2.z, v2.x, 0.f);
//        __m128 mul1 = _mm_mul_ps(va, vb);
//        __m128 mul2 = _mm_mul_ps(vc, vd);
//        __m128 temp = _mm_sub_ps(mul1, mul2);
//
//        _mm_storeu_ps(result.v, temp);
//        vec4 final(result[0], result[1], result[2], 0.f);
//        return final;
//    }
//
//    vec4 operator*(const vec4& other) const {
//        return mul_SSE(other);
//    }
//
//    vec4 operator*(float scalar) const {
//        __m128 vs = _mm_set1_ps(scalar);
//        __m128 va = _mm_load_ps(v);
//        vec4 result;
//        _mm_store_ps(result.v, _mm_mul_ps(va, vs));
//        return result;
//    }
//
//    void divideW() {
//        x /= w;
//        y /= w;
//        z /= w;
//        w = 1.f;
//    }
//
//    float& operator[](const unsigned int index) {
//        return v[index];
//    }
//
//    float operator[](const unsigned int index) const {
//        return v[index];
//    }
//
//    vec4 operator-(const vec4& other) const {
//        return sub_SSE(other);
//    }
//
//    vec4 operator+(const vec4& other) const {
//        return add_SSE(other);
//    }
//
//    vec4 operator/(const vec4& other) const {
//        return div_SSE(other);
//    }
//
//    void normalise() {
//        float length = std::sqrt(x * x + y * y + z * z);
//
//        x /= length;
//        y /= length;
//        z /= length;
//    }
//};
//
