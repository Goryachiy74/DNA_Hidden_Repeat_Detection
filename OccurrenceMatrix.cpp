# include "OccurrenceMatrix.h"

std::unordered_map<char, int> Precompute_DNATab = {
	{'A', 0},
	{'C', 1},
	{'G', 2},
	{'T', 3}
};

// Function to fill the matrix with occurrences of letters at each position in words
std::vector<std::vector<int>> GenerateOccurrenceMatrix(std::string_view sequence, int word_size)
{
	// Initialize a 4 x word_size matrix with zeros
	std::vector<std::vector<int>> matrix(4, std::vector<int>(word_size, 0));

	// Iterate through the sequence, splitting into words
	for (size_t i = 0; i + word_size <= sequence.size(); ++i)
	{
		std::string_view word = sequence.substr(i, word_size); // Extract a word of size 'word_size'

		// Update the matrix for each position in the word
		for (int j = 0; j < word_size; ++j)
		{
			char nucleotide = word[j];
			// Check if the nucleotide is in the Precompute_DNATab
			if (Precompute_DNATab.find(nucleotide) != Precompute_DNATab.end())
			{
				// Increment the count in the matrix
				matrix[Precompute_DNATab[nucleotide]][j]++;
			}
		}
	}

	return matrix;
}


// Function to sum two matrices
std::vector<std::vector<int>> sumMatrices(const std::vector<std::vector<int>>& matrixA,
	const std::vector<std::vector<int>>& matrixB)
{
	// Check if the dimensions of the matrices are the same
	if (matrixA.size() != matrixB.size() || matrixA[0].size() != matrixB[0].size())
	{
		throw std::invalid_argument("Matrices must have the same dimensions for addition.");
	}

	// Initialize the result matrix with the same dimensions
	std::vector<std::vector<int>> result(matrixA.size(), std::vector<int>(matrixA[0].size(), 0));

	// Perform element-wise addition
	for (size_t i = 0; i < matrixA.size(); ++i)
	{
		for (size_t j = 0; j < matrixA[i].size(); ++j) 
		{
			result[i][j] = matrixA[i][j] + matrixB[i][j];
		}
	}

	return result; // Return the resulting matrix
}

// Function to subtract two matrices
std::vector<std::vector<int>> subtractMatrices(const std::vector<std::vector<int>>& matrixA,
	const std::vector<std::vector<int>>& matrixB)
{
	// Check if the dimensions of the matrices are the same
	if (matrixA.size() != matrixB.size() || matrixA[0].size() != matrixB[0].size())
	{
		throw std::invalid_argument("Matrices must have the same dimensions for subtraction.");
	}

	// Initialize the result matrix with the same dimensions
	std::vector<std::vector<int>> result(matrixA.size(), std::vector<int>(matrixA[0].size(), 0));

	// Perform element-wise subtraction
	for (size_t i = 0; i < matrixA.size(); ++i) 
	{
		for (size_t j = 0; j < matrixA[i].size(); ++j) 
		{
			result[i][j] = matrixA[i][j] - matrixB[i][j];
		}
	}

	return result; // Return the resulting matrix
}


vector<vector<int>> AddSequenceToOccurrenceMatrix(const std::vector<std::vector<int>>& prevMatrix, std::string_view sequence, int wordSize)
{
	auto matrixToAdd = GenerateOccurrenceMatrix(sequence, wordSize);
	auto matrixAfterAdd = sumMatrices(prevMatrix, matrixToAdd);

	return matrixAfterAdd;
}

vector<vector<int>> RemoveSequenceFromOccurrenceMatrix(const std::vector<std::vector<int>>& prevMatrix, std::string_view sequence, int wordSize)
{
	auto matrixToSub = GenerateOccurrenceMatrix(sequence, wordSize);
	auto matrixAfterSub = subtractMatrices(prevMatrix, matrixToSub);

	return matrixAfterSub;
}

std::pair<double, std::string> CalculatePercentageSumAndWord(const std::vector<std::vector<int>>& matrix)
{
	double totalSum = 0.0;           // Variable to store the total percentage sum
	std::string representativeWord(matrix[0].size(), ' ');            // String to store the word with the best scores
	int numColumns = matrix[0].size(); // Number of columns (word size)

	char DNATabReverse[4] = { 0 };
	DNATabReverse[0] = 'A';
	DNATabReverse[1] = 'T';
	DNATabReverse[2] = 'C';
	DNATabReverse[3] = 'G';

	// Iterate through each column
	for (int col = 0; col < numColumns; ++col)
	{
		int maxValue = 0;     // Track the highest value in the column
		int columnTotal = 0;  // Sum of values in the column
		char bestLetter = 'A'; // Letter with the highest score in this column

		// Calculate the total and find the max value in the column
		for (int row = 0; row < 4; ++row)
		{
			columnTotal += matrix[row][col]; // Sum up values in the column
			if (matrix[row][col] > maxValue)
			{
				maxValue = matrix[row][col];
				bestLetter = DNATabReverse[row]; // Update the best letter for this column
			}
		}

		// Add the percentage of the max value to the total sum
		if (columnTotal > 0)
		{ // Avoid division by zero
			totalSum += static_cast<double>(maxValue) / columnTotal;
		}

		// Append the best letter for this column to the word
		representativeWord[col] = bestLetter;
	}

	return { totalSum, representativeWord }; // Return both the total sum and the constructed word
}

