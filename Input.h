#pragma once

#include <Eigen/Dense>
#include <Eigen/Eigenvalues> // Include Eigen's Eigenvalues module for eigenvalue and eigenvector computation


Eigen::VectorXd readVectorFromCsv(const std::string& filename) {
	std::ifstream file(filename);
	if (!file.is_open()) {
		std::cerr << "Error: Could not open file " << filename << std::endl;
		return Eigen::VectorXd(); // Return an empty vector
	}

	std::vector<double> data_vec;
	std::string line;

	while (std::getline(file, line)) {
		std::stringstream ss(line);
		std::string cell;

		// Assuming a single column, read the entire line as one value
		if (std::getline(ss, cell, ',')) { // Read until a comma (or end of line for single column)
			try {
				data_vec.push_back(std::stod(cell)); // Convert string to double
			}
			catch (const std::invalid_argument& e) {
				std::cerr << "Warning: Invalid number in CSV file '" << filename << "': " << cell << " (" << e.what() << ")" << std::endl;
				// You might choose to skip or handle this error differently
			}
			catch (const std::out_of_range& e) {
				std::cerr << "Warning: Number out of range in CSV file '" << filename << "': " << cell << " (" << e.what() << ")" << std::endl;
			}
		}
	}

	// Convert std::vector<double> to Eigen::VectorXd
	Eigen::VectorXd eigen_vector(data_vec.size());
	for (size_t i = 0; i < data_vec.size(); ++i) {
		eigen_vector(i) = data_vec[i];
	}

	return eigen_vector;
}
std::vector<double> readNumericCSVColumn(const std::string& filename, int columnIndex = 0)
{
	std::vector<double> data;
	std::ifstream file(filename);

	if (!file.is_open()) {
		throw std::runtime_error("Could not open file: " + filename);
	}



	std::string line;
	//Skip the header
	if (std::getline(file, line)) {}

	while (std::getline(file, line)) {
		std::stringstream lineStream(line);
		std::string cell;
		std::vector<std::string> rowData;

		while (std::getline(lineStream, cell, ',')) {
			rowData.push_back(cell);
		}

		if (columnIndex >= rowData.size()) {
			throw std::runtime_error("Column index out of range in row: " + line);
		}

		try {
			double value = std::stod(rowData[columnIndex]);
			data.push_back(value);
		}
		catch (const std::invalid_argument& e) {
			throw std::runtime_error("Invalid numeric value in row: " + line + ", column: " + std::to_string(columnIndex));
		}
		catch (const std::out_of_range& e) {
			throw std::runtime_error("Numeric value out of range in row: " + line + ", column: " + std::to_string(columnIndex));
		}
	}

	return data;
}
