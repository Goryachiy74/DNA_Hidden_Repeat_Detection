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
using namespace std;
namespace fs = std::filesystem;

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

std::vector<Isochore> detect_isochores(const std::string& dna_sequence, size_t window_size, double gc_threshold);

void saveIsochoresToCsv(std::vector<Isochore> isochores, size_t window_size, double gc_threshold);

std::vector<Isochore> detect_isochores(const std::string& genomeSequence);

void processFASTA2(const std::string& filename, int windowSize, const std::string& outputCSV);

void processSequence(const std::string& sequence, int windowSize, const std::string& outputCSV);

void detect_isochores(const std::string& genomeSequence, const std::string& OutputFolder);

void detect_isochores_optimized(const std::string& genomeSequence, const std::string& OutputFolder, uint64_t windowSize, uint64_t stepSize);

void runIsochoreDetection(const std::string& genomeSequence,
    const std::string& outputFolder,
    uint64_t windowSize, uint64_t stepSize);

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

void saveOverlapsToCSV(const string& filename, const vector<Overlap>& overlaps);

vector<Isochore> loadIsochores(const string& filename);
