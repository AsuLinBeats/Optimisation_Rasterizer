#pragma once

#include <iostream>
#include <chrono>
// ...
// Timing loop with conditionals
template <typename Func>
void TimeTest(Func func) {
	auto start = std::chrono::high_resolution_clock::now();
	auto sum = func();
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "Time doing some work : "
		<< std::chrono::duration<double, std::milli>(end - start).count()
		<< " ms\n";
	std::cout << sum << std::endl;
}


// Check if array is aligned to 32 byte

void CheckAlignment32(void *aligned_array) {
	std::cout << "Address: " << aligned_array << std::endl;
	if (reinterpret_cast<uintptr_t>(aligned_array) % 32 == 0) {
		std::cout << "Aligned to 32 bytes!" << std::endl;
	}
	else {
		std::cout << "Not aligned!" << std::endl;
	}
}



