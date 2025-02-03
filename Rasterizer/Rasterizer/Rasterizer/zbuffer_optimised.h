#pragma once
#include <concepts>

// Zbuffer class for managing depth values during rendering.
// This class is template-constrained to only work with floating-point types (`float` or `double`).

template<std::floating_point T> // Restricts T to be a floating-point type
class Zbuffer {
    T* buffer;                 // Pointer to the buffer storing depth values
    unsigned int width, height; // Dimensions of the Z-buffer

public:

    Zbuffer(unsigned int w, unsigned int h) {
        create(w, h);
    }

    // Default constructor for creating an uninitialized Z-buffer.
    Zbuffer() {
    }
    void create(unsigned int w, unsigned int h) {
        //!maybe we can use multithreads here?
        width = w;
        height = h;
        buffer = new T[width * height]; // Allocate memory for the buffer
    }

    T& operator () (unsigned int x, unsigned int y) {
        return buffer[(y * width) + x]; // Convert 2D coordinates to 1D index
    }

    void clear() {
        size_t total = static_cast<size_t>(width) * height;
        size_t i = 0;

        for (i; i + 8 <= total; i += 8) {
            _mm256_storeu_ps(&buffer[i], _mm256_set1_ps(1.0f));
        }

        for (i; i < total; i++) {
            buffer[i] = 1.0f;
        }
    }

        // Destructor to clean up memory allocated for the Z-buffer.
        ~Zbuffer() {
            delete[] buffer; // Free the allocated memory
        }
    
};
