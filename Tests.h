#pragma once
#include "File_DNA.h"

void CompareLoadSpeed();
void TestFasta();
void MatrixTest();
void SegmentTest(std::string filePath, int minSegmentSize, int wordSize, int lookaheadSize);
void DetectIsochores(std::string filePath, size_t window_size, double gc_threshold);
void DetectIsochores2(std::string filePath, size_t window_size, double gc_threshold);
void processFASTATest(const std::string& inputFile, int windowSize, const std::string& outputFile);
void testZakhariaCode();
void MergeSegmentTest(std::string segmentFilePath, std::string sequanceFilePath, int minSegmentSize, int wordSize, int lookaheadSize);
void ExtractChromosomesTest();
void ReadChromosome(std::string chromosomeFile, std::string outputFolder, int minSegmentSize, int wordSize, int lookaheadSize);
void DetectIsochoresInChromosome(std::string chromosomeFile, std::string outputPath);
void DetectIsochoresInChromosome2(std::string chromosomeFile, std::string outputPath);
void DetectIsochoresInChromosome3(std::string chromosomeFile, std::string outputPath);



