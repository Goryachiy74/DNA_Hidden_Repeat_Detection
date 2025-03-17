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

std::vector<std::tuple<uint64_t, uint64_t, double, std::string>> SegmentDNACostAndWord(
	const std::string& sequence,
	int minSegmentSize,
	int wordSize,
	int lookaheadSize/*,
	const std::function<std::pair<double, std::string>(std::string_view, int)>& costFunction*/
);

void saveSegmentsToCSV(const std::vector<std::tuple<uint64_t, uint64_t, double, std::string>>& segments, const std::string& filename);
std::vector<std::tuple<uint64_t, uint64_t, double, std::string>> loadSegmentsFromCSV(const std::string& filePath);

std::vector<std::tuple<uint64_t, uint64_t, double, std::string>> MergeSimilarSegments(
	const std::vector<std::tuple<uint64_t, uint64_t, double, std::string>>& segments,
	const std::string& sequence,
	int wordSize);
