#pragma once


/// <summary>
/// Function to output a solution matrix to a file
/// </summary>

void outputToFile(Eigen::MatrixXd& solution, std::string& fileName)
{
	// --- Write the non linear solution matrix to a CSV file ---
	std::ofstream outputFile(fileName);
	if (outputFile.is_open()) {
		for (int i = 0; i < solution.rows(); ++i) {
			for (int j = 0; j < solution.cols(); ++j) {
				outputFile << solution(i, j);
				if (j < solution.cols() - 1) {
					outputFile << ","; // Add a comma as a separator
				}
			}
			outputFile << "\n"; // Add a newline character for the next row
		}
		outputFile.close();
		std::cout << fileName << " written to file" << std::endl;
	}
	else {
		std::cerr << "Unable to open file for writing." << std::endl;
	}

}