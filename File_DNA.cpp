#include "File_DNA.h"

#include <vector>

/// <summary>
/// Loads a DNA sequence from a GenBank file.
/// </summary>
/// <param name="filePath">Path to the GenBank file.</param>
/// <returns>The extracted DNA sequence as a string.</returns>
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

/// <summary>
/// Loads a DNA sequence from a FASTA file.
/// </summary>
/// <param name="filename">Path to the FASTA file.</param>
/// <returns>The extracted DNA sequence as a string.</returns>
std::string load_fasta_file(const std::string& filename)
{
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
			//if (c != 'N' && c != 'n')// Exclude 'N'
			//{
				if (isalpha(c))
				{
					sequence += toupper(c);
				}
			//}
		}
		// Append the sequence line to the sequence string
	   // sequence.append(line); // Use append instead of +=
	}

	fasta_file.close();
	return sequence; // Return the complete sequence
}

/// <summary>
/// Saves DNA sequence data to a file.
/// </summary>
/// <param name="filename">Path to the output file.</param>
/// <param name="data">DNA sequence to save.</param>
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

/// <summary>
/// Loads previously saved DNA sequence from a file.
/// </summary>
/// <param name="filename">Path to the file containing saved DNA sequence.</param>
/// <returns>DNA sequence as a string.</returns>
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

/// <summary>
/// Loads previously saved DNA sequence from a file using binary mode.
/// </summary>
/// <param name="filename">Path to the file.</param>
/// <returns>DNA sequence as a string.</returns>
std::string load_previously_saved_data_binary_mode(const std::string& filename)
{
	std::ifstream file(filename, std::ios::in | std::ios::binary); // Open in binary mode
	if (!file) {
		std::cerr << "Error opening file: " << filename << std::endl;
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

/// <summary>
/// Extracts all chromosomes from a FASTA file and saves them as individual files.
/// </summary>
/// <param name="filename">Path to the input file.</param>
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
	std::string outputPath = R"(C:\Braude\Projects\DNA\ncbi_dataset_mouse\Chromosomes\)";

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
			std::string output_filename = outputPath + current_chromosome + ".fna";
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

/// <summary>
/// Reads a specific chromosome from a FASTA file.
/// </summary>
/// <param name="filename">Path to the FASTA file.</param>
/// <param name="chromosome">Name of the chromosome to extract.</param>
/// <returns>DNA sequence of the chromosome as a string.</returns>
std::string read_chromosome(const std::string& filename, const std::string& chromosome)
{
	std::ifstream file(filename);
	if (!file) 
	{
		std::cerr << "Error opening file: " << filename << std::endl;
		return "";
	}

	std::string line, sequence;
	bool found = false;

	while (std::getline(file, line)) 
	{
		if (line[0] == '>')
		{  // New chromosome header found
			if (line.find(chromosome) != std::string::npos) 
			{
				found = true;
				continue;  // Skip the header itself
			}
			else if (found)
			{
				break;  // Stop reading when next chromosome starts
			}
		}
		else if (found) 
		{
			sequence += line;  // Append sequence line
		}
	}

	file.close();

	if (!found) 
	{
		std::cerr << "Chromosome " << chromosome << " not found!" << std::endl;
	}

	return sequence;
}

/// <summary>
/// Reads a chromosome from a file.
/// </summary>
/// <param name="filename">Path to the file.</param>
/// <returns>DNA sequence as a string.</returns>
std::string read_chromosome_file(const std::string& filename)
{
	std::ifstream file(filename);
	if (!file) 
	{
		std::cerr << "Error opening file: " << filename << std::endl;
		return "";
	}

	std::string line, sequence;
	bool first_line = true;

	while (std::getline(file, line)) 
	{
		if (first_line && line[0] == '>') 
		{
			first_line = false;  // Skip the header line
			continue;
		}
		sequence += line;  // Append sequence
	}

	file.close();

	return sequence;
}
