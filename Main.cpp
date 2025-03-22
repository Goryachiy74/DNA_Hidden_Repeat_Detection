#include <iostream> // Required for std::cout and std::cerr
#include <string>   // Required for std::string
#include <cstdlib>  // Required for std::atoi
#include "File_DNA.h"

#include "Isochore.h"
#include "Segment.h"

// Default values for optional parameters
const uint64_t DEFAULT_WINDOW_SIZE = 10000; // Recommended: 50-100 kb
const uint64_t DEFAULT_STEP_SIZE = 10;      // Recommended: ~1 kb

// Function prototypes
void processFullDna(const std::string& filePath, int minSegmentSize, int wordSize, int lookaheadSize, uint64_t windowSize, uint64_t stepSize, const std::string& outputPath);
void processChromosome(const std::string& filePath, int minSegmentSize, int wordSize, int lookaheadSize, uint64_t windowSize, uint64_t stepSize, const std::string& outputPath);

// ======================== Helper Functions ========================
void clearInputBuffer()
{
	std::cin.clear(); // Clear error flags
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Clear buffer
}

std::string getValidatedString(const std::string& prompt)
{
	std::string value;
	while (true)
	{
		std::cout << prompt;
		std::getline(std::cin, value);
		if (!value.empty()) break;
		std::cout << "Invalid input. Try again.\n";
	}
	return value;
}

int getValidatedInt(const std::string& prompt)
{
	int value;
	while (true)
	{
		std::cout << prompt;
		std::cin >> value;
		if (!std::cin.fail()) break;
		std::cout << "Invalid input. Try again.\n";
		clearInputBuffer();
	}
	clearInputBuffer();
	return value;
}

uint64_t getValidatedUInt(const std::string& prompt, uint64_t defaultValue)
{
	uint64_t value;
	std::cout << prompt << " (Default = " << defaultValue << "): ";
	std::string input;
	std::getline(std::cin, input);

	if (input.empty())
	{
		return defaultValue;
	}

	try
	{
		value = std::stoull(input);
		return value;
	}
	catch (...)
	{
		std::cout << "Invalid input. Using default value: " << defaultValue << "\n";
		return defaultValue;
	}
}

// ======================== Help Function ========================
void printHelp(const std::string& programName)
{
	std::cout << "\nUsage: " << programName
		<< " <file_path> <inputType> <minSegmentSize> <wordSize> <lookaheadSize> [windowSize] [stepSize] [outputPath]\n"
		<< "\nParameters:\n"
		<< "  file_path       - Path to the input file (FASTA format)\n"
		<< "  inputType       - Type of input ('fullDna' or 'chromosome')\n"
		<< "  minSegmentSize  - Minimum size of a segment to be processed\n"
		<< "  wordSize        - Size of the word used for analysis\n"
		<< "  lookaheadSize   - Size of the lookahead window\n"
		<< "  windowSize      - Size of sliding window (default = " << DEFAULT_WINDOW_SIZE << ")\n"
		<< "  stepSize        - Step size for sliding window (default = " << DEFAULT_STEP_SIZE << ")\n"
		<< "  outputPath      - (Optional) Path to save the output results\n"
		<< "\nExamples:\n"
		<< "  " << programName << " input_fullDna.fasta fullDna 100 5 10\n"
		<< "  " << programName << " input_chromosome.fasta chromosome 100 5 10 50000 1000 /output/folder\n"
		<< "\nOptions:\n"
		<< "  -h, --help      - Display this help message\n"
		<< std::endl;
}

// ======================== Main Function ========================
int main(int argc, char** argv)
{
	std::string filePath;
	std::string inputType;
	int minSegmentSize = 0;
	int wordSize = 0;
	int lookaheadSize = 0;
	uint64_t windowSize = DEFAULT_WINDOW_SIZE;
	uint64_t stepSize = DEFAULT_STEP_SIZE;
	std::string outputPath;

	// Check for help option
	if (argc == 2 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help")) {
		printHelp(argv[0]);
		return 0;
	}

	// Read provided parameters
	if (argc >= 2) filePath = argv[1];
	if (argc >= 3) inputType = argv[2];
	if (argc >= 4) minSegmentSize = std::atoi(argv[3]);
	if (argc >= 5) wordSize = std::atoi(argv[4]);
	if (argc >= 6) lookaheadSize = std::atoi(argv[5]);
	if (argc >= 7) windowSize = std::stoull(argv[6]); // Optional with default
	if (argc >= 8) stepSize = std::stoull(argv[7]);   // Optional with default
	if (argc == 9) outputPath = argv[8];              // Optional

	// Ask for missing required values
	if (filePath.empty()) filePath = getValidatedString("Enter file path: ");
	if (inputType != "fullDna" && inputType != "chromosome") {
		inputType = getValidatedString("Enter input type ('fullDna' or 'chromosome'): ");
	}
	if (minSegmentSize == 0) minSegmentSize = getValidatedInt("Enter minimum segment size: ");
	if (wordSize == 0) wordSize = getValidatedInt("Enter word size: ");
	if (lookaheadSize == 0) lookaheadSize = getValidatedInt("Enter lookahead size: ");
	if (outputPath.empty()) outputPath = getValidatedString("Enter output path (or press Enter to skip): ");

	// Output the received parameters for verification
	std::cout << "\n=== Parameters ===\n";
	std::cout << "File Path: " << filePath << std::endl;
	std::cout << "Input Type: " << inputType << std::endl;
	std::cout << "Minimum Segment Size: " << minSegmentSize << std::endl;
	std::cout << "Word Size: " << wordSize << std::endl;
	std::cout << "Lookahead Size: " << lookaheadSize << std::endl;
	std::cout << "Window Size: " << windowSize << std::endl;
	std::cout << "Step Size: " << stepSize << std::endl;

	if (inputType == "fullDna") 
	{
		processFullDna(filePath, minSegmentSize, wordSize, lookaheadSize, windowSize, stepSize, outputPath);
	}
	else 
	{
		processChromosome(filePath, minSegmentSize, wordSize, lookaheadSize, windowSize, stepSize, outputPath);
	}

	return 0;
}

// ======================== Sample Processing Functions ========================
void processFullDna(const std::string& filePath, int minSegmentSize, int wordSize, int lookaheadSize, uint64_t windowSize, uint64_t stepSize, const std::string& outputPath)
{
	std::cout << "\n[Processing Full DNA] -> File: " << filePath << std::endl;
	std::cout << "Output Path: " << (outputPath.empty() ? "Not provided" : outputPath) << std::endl;

	string dnaSequence = load_fasta_file(filePath);

	std::cout << "DNA loaded! Size of sequence is  : " << dnaSequence.size() << std::endl;

	std::cout << "Isochore Detection started : " << dnaSequence.size() << std::endl;
	std::cout << "Window size is : " << windowSize << std::endl;
	std::cout << "Step size is : " << stepSize << std::endl;
	runIsochoreDetection(dnaSequence, outputPath, windowSize, stepSize);

	std::cout << "\nThe Word Size is  : " << wordSize << std::endl;
	std::cout << "The Minimum Segment Size is  : " << minSegmentSize * wordSize << " nucleotides" << std::endl;
	std::cout << "The lookahead Size is  : " << lookaheadSize * wordSize << " nucleotides" << std::endl;

	auto segments = SegmentDNACostAndWord(dnaSequence, minSegmentSize, wordSize, lookaheadSize);

	std::string fileName = outputPath + "segments_output_"
		+ std::to_string(minSegmentSize) + "_"
		+ std::to_string(wordSize) + "_"
		+ std::to_string(lookaheadSize) + ".csv";
	saveSegmentsToCSV(segments, fileName);

	std::cout << "Segments saved successfully!: " << fileName << std::endl;

	std::cout << "Number of segments before merge is  : " << segments.size() << std::endl;

	std::cout << "Merge Segments started " << std::endl;

	auto merged = MergeSimilarSegments(segments, dnaSequence, wordSize);

	std::cout << "Number of segments after merge is : " << merged.size() << std::endl;

	std::string mergedFileName = outputPath + "merged_segments_output_"
		+ std::to_string(minSegmentSize) + "_"
		+ std::to_string(wordSize) + "_"
		+ std::to_string(lookaheadSize) + ".csv";

	saveSegmentsToCSV(merged, mergedFileName);

	std::cout << "Merged Segments saved successfully!: " << mergedFileName << std::endl;

	auto result = mergeSegmentsWithGCContent(dnaSequence, merged);

	std::string resultFileName = (fs::path(outputPath) /
		("segments_GcContent_output_" + std::to_string(minSegmentSize) + "_" + std::to_string(wordSize) + "_" + std::to_string(lookaheadSize) + ".csv")).string();

	saveSegmentsGcContentToCsv(result, resultFileName);

	std::cout << "Merged Segments with GC Content saved successfully!: " << resultFileName << std::endl;
}

void processChromosome(const std::string& chromosomeFile, int minSegmentSize, int wordSize, int lookaheadSize, uint64_t windowSize, uint64_t stepSize, const std::string& outputPath)
{
	std::cout << "\n[Processing Chromosome] -> File: " << chromosomeFile << std::endl;
	std::cout << "Output Path: " << (outputPath.empty() ? "Not provided" : outputPath) << std::endl;

	std::cout << "Loading of Chromosome Started from file : " << chromosomeFile << std::endl;
	auto chromosome = read_chromosome_file(chromosomeFile);

	std::cout << "Chromosome loaded! Size of Chromosome is  : " << chromosome.size() << std::endl;

	std::cout << "Isochore Detection started : " << chromosome.size() << std::endl;
	std::cout << "Window size is : " << windowSize << std::endl;
	std::cout << "Step size is : " << stepSize << std::endl;
	runIsochoreDetection(chromosome, outputPath, windowSize, stepSize);


	std::cout << "\nThe Word Size is  : " << wordSize << std::endl;
	std::cout << "The Minimum Segment Size is  : " << minSegmentSize * wordSize << " nucleotides" << std::endl;
	std::cout << "The lookahead Size is  : " << lookaheadSize * wordSize << " nucleotides" << std::endl;

	auto segments = SegmentDNACostAndWord(chromosome, minSegmentSize, wordSize, lookaheadSize);

	std::string fileName = (fs::path(outputPath) /
		("segments_output_" + std::to_string(minSegmentSize) + "_" + std::to_string(wordSize) + "_" + std::to_string(lookaheadSize) + ".csv")).string();

	saveSegmentsToCSV(segments, fileName);

	std::cout << "Segments saved successfully!: " << fileName << std::endl;

	std::cout << "Number of segments before merge is  : " << segments.size() << std::endl;

	std::cout << "Merge Segments started " << std::endl;

	auto merged = MergeSimilarSegments(segments, chromosome, wordSize);

	std::cout << "\nNumber of segments after merge is : " << merged.size() << std::endl;

	std::string mergedFileName = (fs::path(outputPath) /
		("merged_segments_output_" + std::to_string(minSegmentSize) + "_" + std::to_string(wordSize) + "_" + std::to_string(lookaheadSize) + ".csv")).string();

	saveSegmentsToCSV(merged, mergedFileName);

	std::cout << "Merged Segments saved successfully!: " << mergedFileName << std::endl;

	auto result = mergeSegmentsWithGCContent(chromosome, merged);

	std::string resultFileName = (fs::path(outputPath) /
		("segments_GcContent_output_" + std::to_string(minSegmentSize) + "_" + std::to_string(wordSize) + "_" + std::to_string(lookaheadSize) + ".csv")).string();

	saveSegmentsGcContentToCsv(result, resultFileName);

	std::cout << "Merged Segments with GC Content saved successfully!: " << resultFileName << std::endl;

	// Section for creation overlaps and statistics disabled

	//int singleSegmentIsochores = 0;
	//double singleSegmentGCSum = 0;
	//double totalCostSum = 0;
	//double maxCost = 0;
	//double minCost = 0;
	//map<string, int> wordFrequency;

	//std::string isochoresFileName = (fs::path(outputPath) /
	//	("isochores_output_" + std::to_string(windowSize) + "_" + std::to_string(stepSize) + ".csv")).string();

	//auto isochores = loadIsochores(isochoresFileName);

	//// Find overlaps
	//vector<Overlap> overlaps = findIsochoreSegmentOverlap(
	//	isochores, segments, singleSegmentIsochores,
	//	singleSegmentGCSum, totalCostSum, maxCost,
	//	minCost, wordFrequency
	//);

	//std::string overlapsFileName = (fs::path(outputPath) /
	//	("overlaps_isochores_output_" + std::to_string(windowSize) + "_" + std::to_string(stepSize) + ".csv")).string();

	//// Save to CSV
	//saveOverlapsToCSV(overlapsFileName, overlaps);

	//std::string statisticsFileName = (fs::path(outputPath) /
	//	("statistics_overlaps_" + std::to_string(windowSize) + "_" + std::to_string(stepSize) + ".csv")).string();

	//// Save statistics
	//saveStatisticsToFile(statisticsFileName, isochores.size(),
	//	singleSegmentIsochores, singleSegmentGCSum / singleSegmentIsochores,
	//	totalCostSum / overlaps.size(), maxCost, minCost,
	//	wordFrequency.begin()->first, wordFrequency.begin()->second);

}

