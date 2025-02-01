#pragma once
#include <concepts>

// Zbuffer class for managing depth values during rendering.
// This class is template-constrained to only work with floating-point types (`float` or `double`).

template<std::floating_point T> // Restricts T to be a floating-point type
class Zbuffer {
    T* buffer;                 // Pointer to the buffer storing depth values
    unsigned int width, height; // Dimensions of the Z-buffer

public:
    // Constructor to initialize a Z-buffer with the given width and height.
    // Allocates memory for the buffer.
    // Input Variables:
    // - w: Width of the Z-buffer.
    // - h: Height of the Z-buffer.
    Zbuffer(unsigned int w, unsigned int h) {
        create(w, h);
    }

    // Default constructor for creating an uninitialized Z-buffer.
    Zbuffer() {
    }

    // Creates or reinitializes the Z-buffer with the given width and height.
    // Allocates memory for the buffer.
    // Input Variables:
    // - w: Width of the Z-buffer.
    // - h: Height of the Z-buffer.
    void create(unsigned int w, unsigned int h) {
        //!maybe we can use multithreads here
        width = w;
        height = h;
        buffer = new T[width * height]; // Allocate memory for the buffer
    }

    // Accesses the depth value at the specified (x, y) coordinate.
    // Input Variables:
    // - x: X-coordinate of the pixel.
    // - y: Y-coordinate of the pixel.
    // Returns a reference to the depth value at (x, y).
    T& operator () (unsigned int x, unsigned int y) {
        return buffer[(y * width) + x]; // Convert 2D coordinates to 1D index
    }
    //void clear() {
    //    for (unsigned int i = 0; i < width * height; i++) {
    //        buffer[i] = 1.0f; // Reset each depth value
    //    }
    //}

    void clear() {
        // calculate total elements
        size_t total = static_cast<size_t>(width) * height;
        size_t i = 0;

        // address 8 floats every time
        for (i; i + 8 <= total; i += 8) {
            _mm256_storeu_ps(&buffer[i], _mm256_set1_ps(1.0f));
        }

        // address remaining elements using origional methods.
        for (i; i < total; i++) {
            buffer[i] = 1.0f;
        }
    }

    //void clearP() {
    //    //!NEED THINK!
    //    // calculate total elements
    //    size_t total = static_cast<size_t>(width) * height;
    //    size_t blockSize = total / numThreads;
    //    unsigned int numThreads = std::thread::hardware_concurrency(); // get device threads
    //    std::vector<std::thread> threads;
    //    // address 8 floats every time
    //    for (i; i + 8 <= total; i += 8) {
    //        size_t start = t * blockSize;
    //        size_t end = (t == numThreads - 1) ? total : start + blockSize;
    //        {
    //            size_t i = start;

    //            _mm256_storeu_ps(&buffer[i], _mm256_set1_ps(1.0f));
    //            for (i; i < total; i++) {
    //                buffer[i] = 1.0f;
    //            }


    //            for (auto& th : threads) {
    //                th.join();
    //            }
    //               
    //        }
    //    }


        // Destructor to clean up memory allocated for the Z-buffer.
        ~Zbuffer() {
            delete[] buffer; // Free the allocated memory
        }
    
};
