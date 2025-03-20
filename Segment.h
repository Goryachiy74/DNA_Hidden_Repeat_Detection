#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <functional>
#include <thread>
#include <mutex>
#include <future>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include <tuple>
#include <iomanip> // For std::setprecision and std::fixed
#include <algorithm>
#include "OccurrenceMatrix.h"

/// <summary>
/// Segments a DNA sequence based on calculated costs and best words.
/// </summary>
/// <param name="sequence">The DNA sequence to segment.</param>
/// <param name="minSegmentSize">Minimum size of each segment (in words).</param>
/// <param name="wordSize">Size of each word in nucleotides.</param>
/// <param name="lookaheadSize">Number of steps to look ahead when searching for optimal segments.</param>
/// <returns>A vector of tuples containing start, end, cost, and best word for each segment.</returns>
std::vector<std::tuple<uint64_t, uint64_t, double, std::string>> SegmentDNACostAndWord(
	const std::string& sequence,
	int minSegmentSize,
	int wordSize,
	int lookaheadSize/*,
	const std::function<std::pair<double, std::string>(std::string_view, int)>& costFunction*/
);

/// <summary>
/// Saves segmented DNA data to a CSV file.
/// </summary>
/// <param name="segments">Vector of tuples containing segment data.</param>
/// <param name="filename">Path to the output CSV file.</param>
void saveSegmentsToCSV(const std::vector<std::tuple<uint64_t, uint64_t, double, std::string>>& segments, const std::string& filename);

/// <summary>
/// Loads segmented DNA data from a CSV file.
/// </summary>
/// <param name="filePath">Path to the CSV file.</param>
/// <returns>Vector of tuples containing segment data.</returns>
std::vector<std::tuple<uint64_t, uint64_t, double, std::string>> loadSegmentsFromCSV(const std::string& filePath);

/// <summary>
/// Checks if two DNA sequences are cyclic rotations of each other.
/// </summary>
/// <param name="a">First DNA sequence.</param>
/// <param name="b">Second DNA sequence.</param>
/// <returns>True if b is a cyclic rotation of a; otherwise, false.</returns>
std::vector<std::tuple<uint64_t, uint64_t, double, std::string>> MergeSimilarSegments(
	const std::vector<std::tuple<uint64_t, uint64_t, double, std::string>>& segments,
	const std::string& sequence,
	int wordSize);
