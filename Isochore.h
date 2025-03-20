#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <mutex>
#include <cinttypes>
#include <filesystem>
#include <map>
#include <thread>
#include <chrono> 
#include <cstdio>
using namespace std;
namespace fs = std::filesystem;

#ifdef _MSC_VER
#pragma warning(disable : 4244) // Disable int-to-char conversion warning
#pragma warning(disable : 4267) // Disable size_t-to-int conversion warning
#endif

struct Isochore
{
    size_t start;
    size_t end;
    double gc_content;
};

// Struct for Overlap Result
struct Overlap
{
    uint64_t isochore_start;
    uint64_t isochore_end;
    double isochore_gc;
    uint64_t segment_start;
    uint64_t segment_end;
    double segment_cost;
    std:: string best_word;
    uint64_t overlap_length;
};

// Parameters
const uint64_t WINDOW_SIZE = 10000; // or 100000 (50-100 kb recommended)

const uint64_t STEP_SIZE = 10;     // ~1kb steps for precise boundaries detection

/// <summary>
/// Detects isochores in a genome sequence using a sliding window approach.
/// </summary>
/// <param name="genomeSequence">The DNA sequence to analyze.</param>
/// <param name="OutputFolder">Folder to save output files.</param>
/// <param name="windowSize">Size of the sliding window.</param>
/// <param name="stepSize">Step size to slide the window.</param>
void detect_isochores_optimized(const std::string& genomeSequence, const std::string& OutputFolder, uint64_t windowSize, uint64_t stepSize);

/// <summary>
/// Runs the isochore detection in a separate thread and tracks progress.
/// </summary>
/// <param name="genomeSequence">The DNA sequence to analyze.</param>
/// <param name="outputFolder">Folder to save output files.</param>
/// <param name="windowSize">Size of the sliding window.</param>
/// <param name="stepSize">Step size to slide the window.</param>
void runIsochoreDetection(const std::string& genomeSequence,
    const std::string& outputFolder,
    uint64_t windowSize, uint64_t stepSize);

/// <summary>
/// Saves isochores to a CSV file.
/// </summary>
/// <param name="isochores">Vector of isochores.</param>
/// <param name="window_size">Size of the sliding window.</param>
/// <param name="gc_threshold">GC content threshold.</param>
void saveIsochoresToCsv(std::vector<Isochore> isochores, size_t window_size, double gc_threshold);

/// <summary>
/// Loads isochores from a CSV file.
/// </summary>
/// <param name="filename">Path to the CSV file.</param>
/// <returns>Vector of isochores.</returns>
vector<Isochore> loadIsochores(const string& filename);

//Development section
std::vector<Isochore> detect_isochores(const std::string& dna_sequence, size_t window_size, double gc_threshold);

void saveOverlapsToCSV(const string& filename, const vector<Overlap>& overlaps);

std::vector<Isochore> detect_isochores(const std::string& genomeSequence);

void processFASTA2(const std::string& filename, int windowSize, const std::string& outputCSV);

void processSequence(const std::string& sequence, int windowSize, const std::string& outputCSV);

void detect_isochores(const std::string& genomeSequence, const std::string& OutputFolder);

std::vector<Overlap> findIsochoreSegmentOverlap(
    const std::vector<Isochore>& isochores,
    const std::vector<std::tuple<uint64_t, uint64_t, double, std::string>>& segments,
    int& singleSegmentIsochores,
    double& singleSegmentGCSum,
    double& totalCostSum,
    double& maxCost,
    double& minCost,
    std::map<std::string, int>& wordFrequency);

void saveStatisticsToFile(
    const string& filename,
    int totalIsochores,
    int singleSegmentIsochores,
    double avgGCContentSingleSegment,
    double avgCost,
    double maxCost,
    double minCost,
    const string& mostFrequentWord,
    int wordFrequency);


