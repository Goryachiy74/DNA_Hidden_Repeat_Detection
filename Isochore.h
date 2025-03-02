#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <mutex>

struct Isochore {
    size_t start;
    size_t end;
    double gc_content;
};

std::vector<Isochore> detect_isochores(const std::string& dna_sequence, size_t window_size, double gc_threshold);
void saveIsochoresToCsv(std::vector<Isochore> isochores, size_t window_size, double gc_threshold);
std::vector<Isochore> detect_isochores2(const std::string& genomeSequence);

