#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <mutex>
#include <cinttypes>
#include <filesystem>
namespace fs = std::filesystem;

struct Isochore {
    size_t start;
    size_t end;
    double gc_content;
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

