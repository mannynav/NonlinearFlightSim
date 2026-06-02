#pragma once

#include <vector>

std::vector<double> generate_arange(double start, double stop, double step) {
	std::vector<double> result;
	for (double i = start; i < stop; i += step) {
		result.push_back(i);
	}
	return result;
}

double OneDInterp(const std::vector<double>& x,
	const std::vector<double>& y,
	double xi)
{
	if (x.size() != y.size())
		throw std::runtime_error("fastInterp1: x and y size mismatch");

	if (xi <= x.front()) return y.front();
	if (xi >= x.back())  return y.back();

	auto it = std::lower_bound(x.begin(), x.end(), xi);
	size_t i = std::distance(x.begin(), it);
	double x0 = x[i - 1], x1 = x[i];
	double y0 = y[i - 1], y1 = y[i];
	return y0 + (y1 - y0) * (xi - x0) / (x1 - x0);
}