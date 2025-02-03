#pragma once

#include <iostream>
#include <vector>
//#include "vec4.h"
#include "Vec4_Optimised.h"
#include <immintrin.h>
#include <cmath>
// Matrix class for 4x4 transformation matrices
class alignas(16) matrix {
    union {
        float m[4][4]; // 2D array representation of the matrix
        float a[16];   // 1D array representation of the matrix for linear access
    };

public:
    // Default constructor initializes the matrix as an identity matrix
    matrix() {
        identity();
    }

    // Access matrix elements by row and column
    float& operator()(unsigned int row, unsigned int col) { return m[row][col]; }

    // Display the matrix elements in a readable format
    void display() {
        for (unsigned int i = 0; i < 4; i++) {
            for (unsigned int j = 0; j < 4; j++)
                std::cout << m[i][j] << '\t';
            std::cout << std::endl;
        }
    }


    vec4 operator * (const vec4& v) const {
        vec4 result;
        __m128 vec = _mm_load_ps(v.v);

        for (unsigned int i = 0; i < 4; i++) {
            __m128 row = _mm_set_ps(m[i][3], m[i][2], m[i][1], m[i][0]); // represent each row
            __m128 product = _mm_dp_ps(row, vec, 0xFF); // well, avx has built in dot product 
            _mm_store_ss(&result[i], product);
        }
        return result;
    }


    matrix operator * (const matrix& mx) const {
        matrix ret;
 
        // __m256 row0 = _mm256_load_ps();
        for (unsigned int i = 0; i < 4; i++) {
            __m128 row = _mm_load_ps(&this->m[i][0]); 
            for (unsigned int j = 0; j < 4; j++) {
                __m128 col = _mm_set_ps(mx.m[3][j], mx.m[2][j], mx.m[1][j], mx.m[0][j]); 
                __m128 product = _mm_dp_ps(row, col, 0xFF); // Dot product
                _mm_store_ss(&ret.m[i][j], product);
            }
        }
        return ret;
    }


    static matrix makePerspective(float fov, float aspect, float n, float f) {
        
        matrix m;
        m.zero();

        float tanHalfFov = tan(fov / 2.0f);
        float range = f - n;

        m.m[0][0] = 1.0f / (aspect * tanHalfFov);
        m.m[1][1] = 1.0f / tanHalfFov;
        m.m[2][2] = -(f + n) / range;
        m.m[2][3] = -2.0f * f * n / range;
        m.m[3][2] = -1.0f; 

        return m;
    }


    static matrix makeTranslation(float tx, float ty, float tz) {
        matrix result;

        _mm_store_ps(result.a + 0, _mm_setr_ps(1.0f, 0.0f, 0.0f, tx));
        _mm_store_ps(result.a + 4, _mm_setr_ps(0.0f, 1.0f, 0.0f, ty));
        _mm_store_ps(result.a + 8, _mm_setr_ps(0.0f, 0.0f, 1.0f, tz));
        _mm_store_ps(result.a + 12, _mm_setr_ps(0.0f, 0.0f, 0.0f, 1.0f));
        return result;
    }


    static matrix makeRotateZ(float aRad) {
        matrix m;
        // Reduce multiple calculating
        float sinA = sin(aRad);
        float cosA = cos(aRad);
        
        __m128 row0 = _mm_set_ps(0.0f, 0.0f, -sinA, cosA);
        _mm_store_ps(&m.m[0][0], row0); // Third row

        // last row
        __m128 row1 = _mm_set_ps(0.0f, 0.0f, cosA, sinA);
        _mm_store_ps(&m.m[1][0], row1);

        return m;
    }


    static matrix makeRotateX(float aRad) {

        matrix m;
        // Reduce multiple calculating
        float sinA = sin(aRad);
        float cosA = cos(aRad);

        __m128 row1 = _mm_set_ps(0.0f, -sinA, cosA, 0.0f);
        _mm_store_ps(&m.m[1][0], row1);
        
        __m128 row2 = _mm_set_ps(0.0f, cosA, sinA, 0.0f);
        _mm_store_ps(&m.m[2][0], row2);

        return m;
    }


    static matrix makeRotateY(float aRad) {
        matrix m;
        // Reduce multiple calculating
        float sinA = sin(aRad);
        float cosA = cos(aRad);

        __m128 row0 = _mm_set_ps(0.0f, sinA, 0.0f, cosA);
        _mm_store_ps(&m.m[0][0], row0);

        __m128 row2 = _mm_set_ps(0.0f, cosA, 0.0f, -sinA);
        _mm_store_ps(&m.m[2][0], row2);

        return m;

    }


    static matrix makeRotateXYZ(float x, float y, float z) {
        // the original one will invoke 3 functions and there are many overlaps.
       //! Change to this for now
            matrix m;
            float sx = sin(x), cx = cos(x);
            float sy = sin(y), cy = cos(y);
            float sz = sin(z), cz = cos(z);

            
            m.a[0] = cy * cz;
            m.a[1] = -cy * sz;
            m.a[2] = sy;
            m.a[3] = 0.0f;

            m.a[4] = sx * sy * cz + cx * sz;
            m.a[5] = -sx * sy * sz + cx * cz;
            m.a[6] = -sx * cy;
            m.a[7] = 0.0f;

            m.a[8] = -cx * sy * cz + sx * sz;
            m.a[9] = cx * sy * sz + sx * cz;
            m.a[10] = cx * cy;
            m.a[11] = 0.0f;

            m.a[12] = 0.0f;
            m.a[13] = 0.0f;
            m.a[14] = 0.0f;
            m.a[15] = 1.0f;
            
            return m;
        
        /*return matrix::makeRotateX(x) * matrix::makeRotateY(y) * matrix::makeRotateZ(z);*/
    }

    static matrix makeScale(float s) {
        matrix m;
        s = max(s, 0.01f); // Ensure scaling factor is not too small
        // try not invoke identity()
        _mm_store_ps(&m.a[0], _mm_setzero_ps());
        _mm_store_ps(&m.a[4], _mm_setzero_ps());
        _mm_store_ps(&m.a[8], _mm_setzero_ps());
        _mm_store_ps(&m.a[12], _mm_setzero_ps());
        _mm_store_ss(&m.a[0], _mm_set_ss(s));  // a[0] = s
        _mm_store_ss(&m.a[5], _mm_set_ss(s));  // a[5] = s
        _mm_store_ss(&m.a[10], _mm_set_ss(s)); // a[10] = s
        m.a[15] = 1.0f;
        return m;
    }

    // Create an identity matrix
    // Returns an identity matrix
    static matrix makeIdentity() {
        matrix m;
        _mm256_store_ps(&m.a[0], _mm256_setzero_ps());
        _mm256_store_ps(&m.a[8], _mm256_setzero_ps());
        m.a[0] = 1.0f;
        m.a[5] = 1.0f;
        m.a[10] = 1.0f;
        m.a[15] = 1.0f;
        return m;
    }

private:
    // Set all elements of the matrix to 0
    void zero() {
         _mm256_store_ps(&a[0], _mm256_setzero_ps());
        _mm256_store_ps(&a[8], _mm256_setzero_ps());
    }

    // Set the matrix as an identity matrix
    void identity() {
        // To avoid memory jumping by invoking function, don't invoke zero()
        _mm256_store_ps(&a[0], _mm256_setzero_ps());
        _mm256_store_ps(&a[8], _mm256_setzero_ps());
        a[0] = 1.0f;
        a[5] = 1.0f;
        a[10] = 1.0f;
        a[15] = 1.0f;
    }
};


