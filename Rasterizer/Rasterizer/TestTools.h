#pragma once

#include <iostream>
#include <chrono>
// ...
// Timing loop with conditionals
void TimeTest(int (*func)(int), int N ) {
	auto start = std::chrono::high_resolution_clock::now();
	auto sum = func(N);
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "Time doing some work : "
		<< std::chrono::duration<double, std::milli>(end - start).count()
		<< " ms\n";
	std::cout << sum << std::endl;
}

