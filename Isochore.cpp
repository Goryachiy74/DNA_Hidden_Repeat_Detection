#include "Isochore.h"

#include "File_DNA.h"

std::mutex isochoreMtx; // Mutex for thread safety
uint64_t isochoreProgress = static_cast<uint64_t>(0.0);; // Shared progress variable
bool isochoreRunning = true; // Control variable for the progress thread
uint64_t isochoreTotalsize = 100;

/// <summary>
/// Function to update and display progress at regular intervals.
/// </summary>
static void updateProgress()
{
	auto startTime = std::chrono::high_resolution_clock::now(); // Start time

	while (isochoreRunning)
	{
		std::this_thread::sleep_for(std::chrono::seconds(1)); // Update every second

		// Lock the mutex to safely access the progress variable
		std::lock_guard<std::mutex> lock(isochoreMtx);
		double percentage = (static_cast<double>(isochoreProgress) / isochoreTotalsize) * 100;

		// Calculate elapsed time
		auto currentTime = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = currentTime - startTime;

		// Convert elapsed time to hours, minutes, and seconds
		int totalSeconds = static_cast<int>(elapsed.count());
		int hours = totalSeconds / 3600;
		int minutes = (totalSeconds % 3600) / 60;
		int seconds = totalSeconds % 60;

		// Output the current progress and elapsed time
		std::cout << "\rProgress: " << std::fixed << std::setprecision(4)
			<< percentage << "% completed. Elapsed time: "
			<< hours << "h:" << minutes << "m:" << seconds << "s." << std::flush;
	}
}

/// <summary>
/// Calculates whether a base is G or C.
/// </summary>
/// <param name="base">The base character to check.</param>
/// <returns>1 if the base is G or C; otherwise, 0.</returns>
int calculateBaseGC(char base)
{
	return (base == 'G' || base == 'C') ? 1 : 0;
}

/// <summary>
/// Checks if a base is unknown (not A, T, G, or C).
/// </summary>
/// <param name="base">The base character to check.</param>
/// <returns>1 if the base is unknown; otherwise, 0.</returns>
int isUnknownBase(char base)
{
	return (base != 'A' && base != 'T' && base != 'G' && base != 'C');
}

/// <summary>
/// Detects isochores in a genome sequence using a sliding window approach.
/// </summary>
/// <param name="genomeSequence">The DNA sequence to analyze.</param>
/// <param name="OutputFolder">Folder to save output files.</param>
/// <param name="windowSize">Size of the sliding window.</param>
/// <param name="stepSize">Step size to slide the window.</param>
void detect_isochores_optimized(const std::string& genomeSequence, const std::string& OutputFolder, uint64_t windowSize, uint64_t stepSize)
{
	std::string fileName = (fs::path(OutputFolder) /
		("isochores_output_" + std::to_string(windowSize) + "_" + std::to_string(stepSize) + ".csv")).string();

	std::ofstream outfile(fileName);
	outfile << "Start,End,GC_Content\n";

	uint64_t genomeSize = genomeSequence.size();
	isochoreTotalsize = genomeSize;
	// Initialize GC content for the first window
	int gcCount = 0;
	int unknownCount = 0;

	// Initialize GC content and unknown characters for the first window
	for (uint64_t i = 0; i < windowSize; i++)
	{
		gcCount += calculateBaseGC(genomeSequence[i]);
		unknownCount += isUnknownBase(genomeSequence[i]);
	}

	// Only normalize by valid bases (exclude unknown characters)
	double gcPercentage = (unknownCount < static_cast<int>(windowSize))
		? (gcCount / static_cast<double>(windowSize - unknownCount)) * 100.0
		: 0.0;

	outfile << 0 << "," << windowSize << "," << gcPercentage << "\n";

	for (uint64_t pos = stepSize; pos + windowSize <= genomeSize; pos += stepSize)
	{
		// Subtract the base exiting the left side of the window
		for (uint64_t i = pos - stepSize; i < pos; i++)
		{
			gcCount -= calculateBaseGC(genomeSequence[i]);
			unknownCount -= isUnknownBase(genomeSequence[i]);
		}

		// Add the base entering the right side of the window
		for (uint64_t i = pos + windowSize - stepSize; i < pos + windowSize; i++)
		{
			gcCount += calculateBaseGC(genomeSequence[i]);
			unknownCount += isUnknownBase(genomeSequence[i]);
		}

		// Normalize only over valid bases (exclude unknown characters)
		double gcContentPercentage = (unknownCount < windowSize)
			? (gcCount / (double)(windowSize - unknownCount)) * 100.0
			: 0.0;

		outfile << pos << "," << (pos + windowSize) << "," << gcContentPercentage << "\n";

		// Optionally update progress here
		if (pos % (STEP_SIZE * 1000) == 0)
		{
			std::lock_guard<std::mutex> lock(isochoreMtx);
			isochoreProgress = pos;
			isochoreTotalsize = genomeSize;
		}
	}

	outfile.close();
}

/// <summary>
/// Runs the isochore detection in a separate thread and tracks progress.
/// </summary>
/// <param name="genomeSequence">The DNA sequence to analyze.</param>
/// <param name="outputFolder">Folder to save output files.</param>
/// <param name="windowSize">Size of the sliding window.</param>
/// <param name="stepSize">Step size to slide the window.</param>
void runIsochoreDetection(const std::string& genomeSequence,
	const std::string& outputFolder,
	uint64_t windowSize, uint64_t stepSize)
{
	// Start progress thread
	std::thread progressThread(updateProgress);

	detect_isochores_optimized(genomeSequence, outputFolder, windowSize, stepSize);

	// Stop progress thread
	isochoreRunning = false;
	progressThread.join();
}

/// <summary>
/// Merges segments with GC content calculated from the DNA sequence.
/// </summary>
/// <param name="sequence">The full DNA sequence as a string.</param>
/// <param name="segments">Vector of segments (start, end, cost, best word).</param>
/// <returns>
/// A new vector containing segments with an additional GC content field.
/// </returns>
std::vector<std::tuple<uint64_t, uint64_t, double, std::string, double>> mergeSegmentsWithGCContent(
	const std::string& sequence,
	const std::vector<std::tuple<uint64_t, uint64_t, double, std::string>>& segments)
{
	std::vector<std::tuple<uint64_t, uint64_t, double, std::string, double>> result;

	for (const auto& seg : segments) 
	{
		uint64_t start = std::get<0>(seg);
		uint64_t end = std::get<1>(seg);
		double cost = std::get<2>(seg);
		std::string bestWord = std::get<3>(seg);
		int gcCount = 0;
		int unknownCount = 0;
		uint64_t windowSize = end - start;
		// Initialize GC content and unknown characters for the first window
		for (uint64_t i = start; i < end; i++)
		{
			gcCount += calculateBaseGC(sequence[i]);
			unknownCount += isUnknownBase(sequence[i]);
		}

		double gcPercentage = (unknownCount < static_cast<int>(windowSize))
			? (gcCount / static_cast<double>(windowSize - unknownCount)) * 100.0
			: 0.0;

		// Add to result
		result.emplace_back(start, end, cost, bestWord, gcPercentage);
	}

	return result;
}

// Function to find overlap between isochores and segments
std::vector<Overlap> findIsochoreSegmentOverlap(
	const std::vector<Isochore>& isochores,
	const std::vector<std::tuple<uint64_t, uint64_t, double, std::string>>& segments,
	int& singleSegmentIsochores,
	double& singleSegmentGCSum,
	double& totalCostSum,
	double& maxCost,
	double& minCost,
	std::map<std::string, int>& wordFrequency)
{
	std::vector<Overlap> overlaps;
	singleSegmentIsochores = 0;
	singleSegmentGCSum = 0;
	totalCostSum = 0;
	maxCost = -1;
	minCost = 1e9;

	for (const auto& iso : isochores) {
		// Find overlapping segments
		std::vector<std::tuple<uint64_t, uint64_t, double, std::string>> overlappingSegments;
		for (const auto& seg : segments)
		{
			uint64_t seg_start = get<0>(seg);
			uint64_t seg_end = get<1>(seg);

			if (seg_start < iso.end && seg_end > iso.start)
			{
				overlappingSegments.push_back(seg);
			}
		}

		// If exactly one segment creates this isochore
		if (overlappingSegments.size() == 1)
		{
			singleSegmentIsochores++;
			singleSegmentGCSum += iso.gc_content;
		}

		// Add all overlaps to result
		for (const auto& seg : overlappingSegments)
		{
			uint64_t seg_start = get<0>(seg);
			uint64_t seg_end = get<1>(seg);
			double seg_cost = get<2>(seg);
			std::string best_word = get<3>(seg);

			uint64_t overlapStart = std::max(iso.start, seg_start);
			uint64_t overlapEnd = std::min(iso.end, seg_end);
			uint64_t overlapLength = overlapEnd - overlapStart;

			overlaps.push_back({
				iso.start, iso.end, iso.gc_content,
				seg_start, seg_end, seg_cost,
				best_word, overlapLength
				});

			// Track cost stats
			totalCostSum += seg_cost;
			maxCost = std::max(maxCost, seg_cost);
			minCost = std::min(minCost, seg_cost);

			// Track best word frequency
			wordFrequency[best_word]++;
		}
	}

	return overlaps;
}

/// <summary>
/// Saves isochores to a CSV file.
/// </summary>
/// <param name="isochores">Vector of isochores.</param>
/// <param name="window_size">Size of the sliding window.</param>
/// <param name="gc_threshold">GC content threshold.</param>
void saveIsochoresToCsv(std::vector<Isochore> isochores, size_t window_size, double gc_threshold)
{
	// Save results to CSV file
	std::string fileName = "isochores_"
		+ std::to_string(window_size) + "_"
		+ std::to_string(gc_threshold) + ".csv";


	std::ofstream csv_file(fileName);
	if (csv_file.is_open())
	{
		csv_file << "Start,End,GC_Content\n"; // CSV header
		for (const Isochore& isochore : isochores)
		{
			csv_file << isochore.start << "," << isochore.end << "," << isochore.gc_content << "\n";
		}
		csv_file.close();
		std::cout << "Isochores saved to isochores.csv" << std::endl;
	}
	else
	{
		std::cerr << "Unable to open file for writing." << std::endl;
	}
}

/// <summary>
/// Loads isochores from a CSV file.
/// </summary>
/// <param name="filename">Path to the CSV file.</param>
/// <returns>Vector of isochores.</returns>
vector<Isochore> loadIsochores(const string& filename)
{
	vector<Isochore> isochores;
	ifstream file(filename);

	if (!file.is_open()) {
		cerr << "Error: Unable to open file " << filename << "\n";
		return isochores;
	}

	string line;
	getline(file, line); // Skip header line

	while (getline(file, line))
	{
		stringstream ss(line);
		string field;
		Isochore iso;

		getline(ss, field, ',');
		iso.start = stoull(field);

		getline(ss, field, ',');
		iso.end = stoull(field);

		getline(ss, field, ',');
		iso.gc_content = stod(field);

		isochores.push_back(iso);
	}

	file.close();
	return isochores;
}


// Development section

// ---- 2. Save Overlaps to CSV ----
void saveOverlapsToCSV(const string& filename, const vector<Overlap>& overlaps)
{
	ofstream file(filename);

	if (!file.is_open())
	{
		cerr << "Error: Unable to create file " << filename << "\n";
		return;
	}

	// Write header
	file << "Isochore Start,Isochore End,Isochore GC,Segment Start,Segment End,Segment Cost,Best Word,Overlap Length\n";

	// Write data
	for (const auto& o : overlaps) {
		file << o.isochore_start << ","
			<< o.isochore_end << ","
			<< o.isochore_gc << ","
			<< o.segment_start << ","
			<< o.segment_end << ","
			<< o.segment_cost << ","
			<< o.best_word << ","
			<< o.overlap_length << "\n";
	}

	file.close();
	cout << "Overlaps saved to " << filename << "\n";
}

// ---- 3. Save Statistics to File ----
void saveStatisticsToFile(
	const string& filename,
	int totalIsochores,
	int singleSegmentIsochores,
	double avgGCContentSingleSegment,
	double avgCost,
	double maxCost,
	double minCost,
	const string& mostFrequentWord,
	int wordFrequency)
{
	ofstream file(filename);

	if (!file.is_open()) {
		cerr << "Error: Unable to create file " << filename << "\n";
		return;
	}

	file << "--- Statistics ---\n";
	file << "Total Isochores: " << totalIsochores << "\n";
	file << "Isochores Created by Single Segment: " << singleSegmentIsochores << "\n";
	file << "Percentage of Single Segment Isochores: "
		<< (singleSegmentIsochores * 100.0) / totalIsochores << "%\n";
	file << "Average GC Content of Single Segment Isochores: " << avgGCContentSingleSegment << "\n";
	file << "Average Segment Cost: " << avgCost << "\n";
	file << "Max Segment Cost: " << maxCost << "\n";
	file << "Min Segment Cost: " << minCost << "\n";
	file << "Most Frequent Best Word: " << mostFrequentWord
		<< " (" << wordFrequency << " times)\n";

	file.close();
	cout << "Statistics saved to " << filename << "\n";
}


double calculate_gc_content(const std::string& sequence)
{
	double gc_count = 0;
	for (char nucleotide : sequence)
	{
		if (nucleotide == 'G' || nucleotide == 'C')
		{
			gc_count++;
		}
	}
	return gc_count / sequence.size();
}

// Function to calculate GC content
double calculateGCContent(std::string_view sequence) {
	int gcCount = 0;
	for (char base : sequence)
	{
		if (base == 'G' || base == 'C')
		{
			gcCount++;
		}
	}
	return (sequence.length() > 0) ? (gcCount * 100.0 / sequence.length()) : 0;
}

int calculateGCContent(const char* arr, uint64_t pos, uint64_t windowSize)
{
	int gcCount = 0;
	for (uint64_t i = pos; i < pos + windowSize; i++)
	{
		char base = arr[i];
		if (base == 'G' || base == 'C' || base == 'g' || base == 'c')
		{
			gcCount++;
		}
	}
	return  gcCount;
}

// Function to calculate GC and CG pairs in a given sequence
std::pair<double, size_t> calculate_gc_pairs(const std::string& sequence)
{
	size_t gc_count = 0; // Count of GC and CG pairs
	size_t total_pairs = 0; // Total pairs in the sequence

	for (size_t i = 0; i < sequence.length() - 1; ++i)
	{
		if ((sequence[i] == 'G' && sequence[i + 1] == 'C') ||
			(sequence[i] == 'C' && sequence[i + 1] == 'G'))
		{
			gc_count++;
		}
		total_pairs++;
	}

	return { static_cast<double>(gc_count), total_pairs };
}

std::vector<Isochore> detect_isochores(const std::string& dna_sequence, size_t window_size, double gc_threshold)
{
	std::vector<Isochore> isochores;
	size_t length = dna_sequence.length();
	std::thread progressThread(updateProgress);


	// Initialize GC pairs count for the first window
	auto [gc_count, total_pairs] = calculate_gc_pairs(dna_sequence.substr(0, window_size));



	for (size_t i = 0; i <= length - window_size; ++i)
	{
		// Update progress every 1000 iterations to reduce mutex locking
		if (i % 1000 == 0)
		{
			std::lock_guard<std::mutex> lock(isochoreMtx);
			isochoreProgress = i; // Update progress
			isochoreTotalsize = length - window_size;
		}


		// Calculate GC content for the current window
		double gc_content = (total_pairs > 0) ? (gc_count / total_pairs) : 0.0;


		if (gc_content >= gc_threshold)
		{
			// Check if this is the start of a new isochore
			if (isochores.empty() || isochores.back().end < i)
			{
				isochores.push_back({ i, i + window_size - 1, gc_content });
			}
			else
			{
				// Extend the last isochore
				isochores.back().end = i + window_size - 1;
			}
		}

		// Slide the window: update GC pairs count
		if (i + window_size < length) // Ensure we don't go out of bounds
		{
			// Remove the pair going out of the window
			if (dna_sequence[i] == 'G' && dna_sequence[i + 1] == 'C')
				gc_count--;
			else if (dna_sequence[i] == 'C' && dna_sequence[i + 1] == 'G')
				gc_count--;

			// Add the new pair coming into the window
			if (dna_sequence[i + window_size - 1] == 'G' && dna_sequence[i + window_size] == 'C')
				gc_count++;
			else if (dna_sequence[i + window_size - 1] == 'C' && dna_sequence[i + window_size] == 'G')
				gc_count++;

			// Update total pairs (only if the window is valid)
			total_pairs = window_size - 1; // Total pairs in the window
		}
	}

	// Stop the progress thread
	isochoreRunning = false;
	progressThread.join(); // Wait for the progress thread to finish
	return isochores;
}

void detect_isochores(const std::string& genomeSequence, const std::string& OutputFolder)
{
	std::thread progressThread(updateProgress);
	int gcContentCount = calculateGCContent(genomeSequence.c_str(), 0, WINDOW_SIZE);

	std::string fileName = OutputFolder + "isochores_"
		+ std::to_string(WINDOW_SIZE) + "_"
		+ std::to_string(STEP_SIZE) + ".csv";
	std::ofstream csv_file(fileName);
	csv_file << "Start,End,GC_Content\n"; // CSV header
	char buf[128]{ 0 };

	for (size_t i = 0; i + WINDOW_SIZE + STEP_SIZE <= genomeSequence.size(); i += STEP_SIZE)
	{
		int gcPrefix = calculateGCContent(genomeSequence.c_str(), i, STEP_SIZE);
		int gcSuffix = calculateGCContent(genomeSequence.c_str(), i + WINDOW_SIZE, STEP_SIZE);

		gcContentCount = gcContentCount - gcPrefix + gcSuffix;
		int gcContent = (gcContentCount * 100) / WINDOW_SIZE;


		// Update progress every 1000 iterations to reduce mutex locking
		if (i % 1000 == 0)
		{
			std::lock_guard<std::mutex> lock(isochoreMtx);
			isochoreProgress = i; // Update progress
			isochoreTotalsize = genomeSequence.size();
		}

		//part of the C++ standard library and cross-platform
		snprintf(buf, sizeof(buf), "%" PRIu64 ", %" PRIu64 ", %d\n", i, i + WINDOW_SIZE, gcContent);

		csv_file << buf;

	}
	// Stop the progress thread
	isochoreRunning = false;
	progressThread.join(); // Wait for the progress thread to finish
}

// Function to detect isochores
std::vector<Isochore> detect_isochores(const std::string& genomeSequence)
{
	std::vector<Isochore> isochores;
	std::thread progressThread(updateProgress);
	int gcContentCount = calculateGCContent(genomeSequence.c_str(), 0, WINDOW_SIZE);
	int gcContent = (gcContentCount * 100) / WINDOW_SIZE;
	Isochore iso;
	iso.start = 0;
	iso.end = WINDOW_SIZE;
	iso.gc_content = gcContent;
	isochores.push_back(iso);

	std::string fileName = "isochores_"
		+ std::to_string(WINDOW_SIZE) + "_"
		+ std::to_string(STEP_SIZE) + ".csv";
	std::ofstream csv_file(fileName);
	csv_file << "Start,End,GC_Content\n"; // CSV header
	char buf[128]{ 0 };

	for (size_t i = 0; i + WINDOW_SIZE + STEP_SIZE <= genomeSequence.size(); i += STEP_SIZE)
	{
		int gcPrefix = calculateGCContent(genomeSequence.c_str(), i, STEP_SIZE);
		int gcSuffix = calculateGCContent(genomeSequence.c_str(), i + WINDOW_SIZE, STEP_SIZE);

		gcContentCount = gcContentCount - gcPrefix + gcSuffix;
		gcContent = (gcContentCount * 100) / WINDOW_SIZE;


		// Update progress every 1000 iterations to reduce mutex locking
		if (i % 1000 == 0)
		{
			std::lock_guard<std::mutex> lock(isochoreMtx);
			isochoreProgress = i; // Update progress
			isochoreTotalsize = genomeSequence.size();
		}

		// part of the C++ standard library and cross-platform
		snprintf(buf, sizeof(buf), "%" PRIu64 ", %" PRIu64 ", %d\n", i, i + WINDOW_SIZE, gcContent);
		csv_file << buf;

	}
	// Stop the progress thread
	isochoreRunning = false;
	progressThread.join(); // Wait for the progress thread to finish
	return isochores;
}

// Function to calculate GC content in a DNA sequence
double calculateGCContent2(const std::string& sequence)
{
	int gcCount = 0, totalCount = 0;
	for (char base : sequence)
	{
		if (base == 'G' || base == 'C' || base == 'g' || base == 'c')
		{
			gcCount++;
		}
		if (isalpha(base))
		{ // Ignore non-DNA characters
			totalCount++;
		}
	}
	return totalCount > 0 ? (gcCount * 100.0 / totalCount) : 0.0;
}

// Function to process FASTA file and compute GC content
void processFASTA2(const std::string& filename, int windowSize, const std::string& outputCSV) {
	std::ifstream file(filename);
	std::ofstream outputFile(outputCSV);

	if (!file.is_open() || !outputFile.is_open()) {
		std::cerr << "Error opening file!" << std::endl;
		return;
	}

	std::string line, sequence = "";
	long long position = 0;
	isochoreTotalsize = 0;

	// Compute total genome size
	while (getline(file, line)) {
		if (!line.empty() && line[0] != '>') { // Ignore headers
			isochoreTotalsize += line.size();
		}
	}
	file.clear();
	file.seekg(0, std::ios::beg); // Reset file pointer

	// Write CSV Header
	outputFile << "Position,GC_Content\n";

	// Start progress tracking thread
	std::thread progressThread(updateProgress);

	while (getline(file, line)) {
		if (line.empty()) continue;
		if (line[0] == '>') continue; // Skip header

		sequence += line; // Append sequence data

		while (sequence.size() >= static_cast<size_t>(windowSize))
		{
			std::string window = sequence.substr(0, windowSize);
			double gcContent = calculateGCContent(window);

			outputFile << position << "," << gcContent << "\n";

			position += windowSize;

			// Update progress
			{
				std::lock_guard<std::mutex> lock(isochoreMtx);
				isochoreProgress += windowSize;
			}

			sequence = sequence.substr(windowSize);
		}
	}

	// Stop progress thread
	isochoreRunning = false;
	progressThread.join();

	file.close();
	outputFile.close();
	std::cout << "\nProcessing complete! Output saved in " << outputCSV << std::endl;
}

// Function to process the DNA sequence
void processSequence(const std::string& sequence, int windowSize, const std::string& outputCSV)
{
	std::ofstream outputFile(outputCSV);
	if (!outputFile.is_open()) {
		std::cerr << "Error opening output file!" << std::endl;
		return;
	}

	long long position = 0;
	isochoreTotalsize = sequence.size(); // Set total size for progress tracking

	// Write CSV Header
	outputFile << "Position,GC_Content\n";

	// Start progress tracking thread
	std::thread progressThread(updateProgress);

	for (size_t i = 0; i + windowSize <= sequence.size(); i += windowSize) {
		std::string window = sequence.substr(i, windowSize);
		double gcContent = calculateGCContent(window);

		outputFile << position << "," << gcContent << "\n";
		position += windowSize;

		// Update progress
		{
			std::lock_guard<std::mutex> lock(isochoreMtx);
			isochoreProgress = position;
		}
	}

	// Stop progress thread
	isochoreRunning = false;
	progressThread.join();

	outputFile.close();
	std::cout << "\nProcessing complete! Output saved in " << outputCSV << std::endl;
}