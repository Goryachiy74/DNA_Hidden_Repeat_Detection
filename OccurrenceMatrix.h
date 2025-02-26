#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>

using namespace std;

std::vector<std::vector<int>> GenerateOccurrenceMatrix(std::string_view sequence, int word_size);

std::vector<std::vector<int>> sumMatrices(const std::vector<std::vector<int>>& matrixA,
	const std::vector<std::vector<int>>& matrixB);

std::vector<std::vector<int>> subtractMatrices(const std::vector<std::vector<int>>& matrixA,
	const std::vector<std::vector<int>>& matrixB);

vector<vector<int>> AddSequenceToOccurrenceMatrix(const std::vector<std::vector<int>>& prevMatrix, std::string_view sequence, int wordSize);

vector<vector<int>> RemoveSequenceFromOccurrenceMatrix(const std::vector<std::vector<int>>& prevMatrix, std::string_view sequence, int wordSize);

std::pair<double, std::string> CalculatePercentageSumAndWord(const std::vector<std::vector<int>>& matrix);
