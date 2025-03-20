#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>

using namespace std;

/// <summary>
/// Generates an occurrence matrix for DNA sequences.
/// </summary>
/// <param name="sequence">The DNA sequence to analyze.</param>
/// <param name="word_size">The size of the word to split the sequence into.</param>
/// <returns>A 4 x word_size matrix representing nucleotide occurrences.</returns>
std::vector<std::vector<int>> GenerateOccurrenceMatrix(std::string_view sequence, int word_size);

/// <summary>
/// Adds two matrices element-wise.
/// </summary>
/// <param name="matrixA">The first matrix.</param>
/// <param name="matrixB">The second matrix to add.</param>
/// <returns>The resulting matrix after addition.</returns>
std::vector<std::vector<int>> sumMatrices(const std::vector<std::vector<int>>& matrixA,
	const std::vector<std::vector<int>>& matrixB);

/// <summary>
/// Subtracts two matrices element-wise.
/// </summary>
/// <param name="matrixA">The first matrix.</param>
/// <param name="matrixB">The matrix to subtract.</param>
/// <returns>The resulting matrix after subtraction.</returns>
std::vector<std::vector<int>> subtractMatrices(const std::vector<std::vector<int>>& matrixA,
	const std::vector<std::vector<int>>& matrixB);

/// <summary>
/// Adds a new DNA sequence's occurrence matrix to an existing one.
/// </summary>
/// <param name="prevMatrix">The previous occurrence matrix.</param>
/// <param name="sequence">The new DNA sequence to add.</param>
/// <param name="wordSize">The word size for the matrix.</param>
/// <returns>The updated occurrence matrix.</returns>
vector<vector<int>> AddSequenceToOccurrenceMatrix(const std::vector<std::vector<int>>& prevMatrix, std::string_view sequence, int wordSize);

/// <summary>
/// Removes a DNA sequence's occurrence matrix from an existing one.
/// </summary>
/// <param name="prevMatrix">The previous occurrence matrix.</param>
/// <param name="sequence">The DNA sequence to remove.</param>
/// <param name="wordSize">The word size for the matrix.</param>
/// <returns>The updated occurrence matrix after removal.</returns>
vector<vector<int>> RemoveSequenceFromOccurrenceMatrix(const std::vector<std::vector<int>>& prevMatrix, std::string_view sequence, int wordSize);

/// <summary>
/// Calculates the total percentage sum of the maximum nucleotide occurrences per column
/// and constructs the representative word based on highest frequency nucleotides.
/// </summary>
/// <param name="matrix">The occurrence matrix.</param>
/// <returns>A pair consisting of the total percentage sum and the representative word.</returns>
std::pair<double, std::string> CalculatePercentageSumAndWord(const std::vector<std::vector<int>>& matrix);
