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

void outputToFileWithHeaders(Eigen::MatrixXd& solution, std::string& fileName)
{

	std::ofstream outputFile(fileName);
	if (outputFile.is_open()) {
		outputFile << std::fixed << std::setprecision(10); // Set precision for output

		// Write header row (optional, but helpful)
		outputFile << "u_mps,v_mps,w_mps,p_rps,q_rps,r_rps,phi_rad,theta_rad,psi_rad,x_e,y_e,z_e\n";

		// Write each row of the solution matrix to the CSV file
		for (int i = 0; i < solution.cols(); ++i) {
			for (int j = 0; j < solution.rows(); ++j) {
				outputFile << solution(j, i);
				if (j < solution.rows() - 1) {
					outputFile << ","; // Add comma as separator
				}
			}
			outputFile << "\n"; // Add newline after each row
		}
		outputFile.close();
		std::cout << fileName << " written to file" << std::endl;
	}
	else {
		std::cerr << "Unable to open file for writing: solution.csv" << std::endl;

	}



}

