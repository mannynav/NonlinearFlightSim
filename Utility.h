#pragma once

#include <vector>

std::vector<double> generate_arange(double start, double stop, double step) {
	std::vector<double> result;
	for (double i = start; i < stop; i += step) {
		result.push_back(i);
	}
	return result;
}