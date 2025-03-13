#include "File_DNA.h"

#include <vector>

//GenBank files contain metadata and sequence data. The actual DNA sequence is usually in the ORIGIN section.
string load_gen_bank_file(const string& filePath)
{
	ifstream file(filePath);
	if (!file.is_open())
	{
		cerr << "Error: Unable to open file " << filePath << endl;
		exit(1);
	}

	string line, sequence;
	bool inOrigin = false;

	while (getline(file, line))
	{
		// Look for the ORIGIN section
		if (line.find("ORIGIN") != string::npos)
		{
			inOrigin = true;
			continue;
		}

		if (inOrigin)
		{
			// Extract sequence data after ORIGIN
			if (line.find("//") != string::npos) break; // End of sequence
			for (char c : line)
			{
				if (c != 'N' && c != 'n')// Exclude 'N'
				{
					if (isalpha(c))
					{
						sequence += toupper(c);
					}
				}
			}
		}
	}

	file.close();
	return sequence;
}

std::string load_fasta_file(const std::string& filename) {
	std::ifstream fasta_file(filename);
	std::string sequence;

	if (!fasta_file.is_open()) {
		std::cerr << "Error: Could not open the file " << filename << std::endl;
		return ""; // Return an empty string on error
	}

	// Reserve an initial size for the sequence string
	//sequence.reserve(10'000'000); // Adjust this size based on expected sequence length

	std::string line;
	while (std::getline(fasta_file, line)) {
		// Skip the header line (starts with '>')
		if (line.empty() || line[0] == '>') {
			continue;
		}
		for (char c : line)
		{
			if (c != 'N' && c != 'n')// Exclude 'N'
			{
				if (isalpha(c))
				{
					sequence += toupper(c);
				}
			}
		}
		// Append the sequence line to the sequence string
	   // sequence.append(line); // Use append instead of +=
	}

	fasta_file.close();
	return sequence; // Return the complete sequence
}

void save_loaded_data_as_file(const std::string& filename, const std::string& data)
{
	// Create an output file stream
	std::ofstream outFile(filename);

	// Check if the file is open
	if (outFile.is_open()) {
		// Write the string to the file
		outFile << data;

		// Close the file
		outFile.close();
		std::cout << "String saved to output.txt successfully." << '\n';
	}
	else {
		std::cerr << "Unable to open file for writing." << '\n';
	}
}

std::string load_previously_saved_data(const std::string& filename)
{
	// Create an input file stream
	std::ifstream inFile(filename);
	std::string data;

	// Check if the file is open
	if (inFile.is_open())
	{

		// Read the string from the file
		// Using getline to read the entire line
		std::getline(inFile, data);

		// Close the file
		inFile.close();

		// Output the loaded string
		std::cout << "DNA loaded from file: " << filename << std::endl;
	}
	else
	{
		std::cerr << "Unable to open file for reading." << std::endl;
	}

	return  data;

}

std::string load_previously_saved_data2(const std::string& filePath)
{
	std::ifstream file(filePath, std::ios::in | std::ios::binary); // Open in binary mode
	if (!file) {
		std::cerr << "Error opening file: " << filePath << std::endl;
		return ""; // Return an empty string on error
	}

	// Buffer to hold chunks of data
	const std::size_t bufferSize = 8192; // 8 KB buffer
	std::vector<char> buffer(bufferSize);
	std::string result; // String to hold the file contents

	// Read the file in chunks
	while (file.read(buffer.data(), buffer.size()) || file.gcount() > 0) {
		std::size_t bytesRead = file.gcount(); // Get the number of bytes read
		result.append(buffer.data(), bytesRead); // Append chunk to result string
	}

	file.close(); // Close the file
	return result; // Return the loaded string
}

void extract_all_chromosomes(const std::string& filename)
{
	std::ifstream file(filename);
	if (!file)
	{
		std::cerr << "Error opening file: " << filename << std::endl;
		return;
	}

	std::ofstream outFile;
	std::string line;
	std::string current_chromosome;

	while (std::getline(file, line))
	{
		if (line[0] == '>') 
		{  // Found a new chromosome header
			if (outFile.is_open()) 
			{
				outFile.close();  // Close previous chromosome file
			}

			// Extract chromosome name (remove '>' and take first word)
			size_t space_pos = line.find(' ');
			current_chromosome = line.substr(1, space_pos - 1);

			// Open a new file for this chromosome
			std::string output_filename = current_chromosome + ".fna";
			outFile.open(output_filename);
			if (!outFile) 
			{
				std::cerr << "Error creating file: " << output_filename << std::endl;
				continue;
			}

			outFile << line << std::endl; // Write header to file
			std::cout << "Saving: " << output_filename << std::endl;
		}
		else if (outFile.is_open()) 
		{
			outFile << line << std::endl; // Write sequence to file
		}
	}

	// Close the last open file
	if (outFile.is_open())
	{
		outFile.close();
	}

	file.close();
}
