#include "Tests.h"
#include "OccurrenceMatrix.h"
#include "Segment.h"
#include "Isochore.h"
#include "GenomDNA.h"

// Function to measure execution time
template <typename Func>
double measureExecutionTime(Func func) {
	auto start = std::chrono::high_resolution_clock::now(); // Start time
	func(); // Call the function
	auto end = std::chrono::high_resolution_clock::now(); // End time

	// Calculate duration in milliseconds
	std::chrono::duration<double, std::milli> duration = end - start;
	return duration.count(); // Return duration in milliseconds
}

// Function to display the occurrence matrix
void displayMatrix(const std::vector<std::vector<int>>& matrix) {
	std::cout << "Occurrence Matrix:" << std::endl;

	// Iterate through each row
	for (const auto& row : matrix) {
		// Iterate through each column in the row
		for (int count : row) {
			std::cout << count << " "; // Print the count
		}
		std::cout << std::endl; // Move to the next line after each row
	}
}

void TestFasta()
{
	string fastaFile = R"(C:\Braude\Projects\ncbi_dataset\ncbi_dataset\data\GCF_000001405.40\GCF_000001405.40_GRCh38.p14_genomic.fna)";

	std::cout << "Loading of DNA Started from file : " << fastaFile << std::endl;

	string dnaSequence = load_fasta_file(fastaFile);

	std::cout << "DNA loaded! Size of sequence is  : " << dnaSequence.size() << std::endl;

	string fastaFileToSave = R"(C:\Braude\Projects\DNA\GCF_000001405.40_GRCh38.p14_genomic.txt)";

	std::cout << "Saving of DNA file started to : " << fastaFileToSave << std::endl;

	save_loaded_data_as_file(fastaFileToSave, dnaSequence);

	std::cout << "DNA saved! " << std::endl;

	std::cout << "Loading of DNA Started from file : " << fastaFileToSave << std::endl;

	string dnaSequenceFromSaved = load_previously_saved_data(fastaFileToSave);

	std::cout << "DNA loaded! Size of sequence is  : " << dnaSequenceFromSaved.size() << std::endl;

	if (dnaSequence.size() != dnaSequenceFromSaved.size())
	{
		std::cout << "Sequences are different" << std::endl;
	}
	if (dnaSequence[100] != dnaSequenceFromSaved[100])
	{
		std::cout << "Sequences are different" << std::endl;
	}
}

void ExtractChromosomesTest()
{
	string fastaFile = R"(C:\Braude\Projects\ncbi_dataset\ncbi_dataset\data\GCF_000001405.40\GCF_000001405.40_GRCh38.p14_genomic.fna)";

	std::cout << "Loading of DNA Started from file : " << fastaFile << std::endl;

	extract_all_chromosomes(fastaFile);
}

void ReadChromosome()
{
	string fastaFile = R"(C:\Braude\Projects\DNA\Chromosome\NC_000001.11.fna)";
	std::cout << "Loading of Chromosome Started from file : " << fastaFile << std::endl;
	auto seq = read_chromosome_file(fastaFile);

	std::cout << "DNA loaded! Size of sequence is  : " << seq.size() << std::endl;

	//int minSegmentSize = 2000;
	//int word_size = 5;
	//int lookaheadSize = 20000;
	//int seqSize = 4000000;
	//std::string firstSeq = dnaSequenceFromSaved.substr(0, seqSize);
	std::cout << "The Sequence size is  : " << seq.size() << std::endl;
	std::cout << "The Word Size is  : " << 3 << std::endl;
	std::cout << "The Minimum Segment Size is  : " << 6000 << " nucleotides" << std::endl;
	std::cout << "The lookahead Size is  : " << 60000 << " nucleotides" << std::endl;

	auto segments = SegmentDNACostAndWord(seq, 2000, 3, 20000);

	// Print the resulting segments, costs, and best words
	std::cout << "Segmented DNA sequence with costs and best words:\n";
	for (const auto& [start, end, cost, bestWord] : segments) {
		std::cout << "Segment: [" << start << ", " << end << ") -> "
			<< " len: " << end - start
			<< " | Cost: " << cost
			<< " | Best Word: " << bestWord << "\n";
	}

	std::string fileName = "C:\\Braude\\Projects\\DNA\\segments_output_"
		+ std::to_string(2000) + "_"
		+ std::to_string(3) + "_"
		+ std::to_string(20000) + ".csv";
	saveSegmentsToCSV(segments, fileName);



}

void CompareLoadSpeed()
{
	string fastaFileToSave = R"(C:\Braude\Projects\DNA\GCF_000001405.40_GRCh38.p14_genomic.txt)";

	std::cout << "(1)Loading of DNA Started from file : " << fastaFileToSave << std::endl;
	auto start = std::chrono::high_resolution_clock::now(); // Start time
	string dnaSequenceFromSaved = load_previously_saved_data(fastaFileToSave);
	auto end = std::chrono::high_resolution_clock::now(); // End time
	std::chrono::duration<double, std::milli> duration = end - start;
	double timeA = duration.count();
	std::cout << "(1)DNA loaded! Size of sequence is  : " << dnaSequenceFromSaved.size() << std::endl;

	std::cout << "(2)Loading of DNA Started from file : " << fastaFileToSave << std::endl;
	auto start2 = std::chrono::high_resolution_clock::now(); // Start time
	string dnaSequenceFromSaved2 = load_previously_saved_data2(fastaFileToSave);
	auto end2 = std::chrono::high_resolution_clock::now(); // End time
	std::chrono::duration<double, std::milli> duration2 = end2 - start2;
	double timeB = duration2.count();
	std::cout << "(2)DNA loaded! Size of sequence is  : " << dnaSequenceFromSaved2.size() << std::endl;


	if (dnaSequenceFromSaved.size() != dnaSequenceFromSaved2.size())
	{
		std::cout << "Sequences are different" << std::endl;
	}
	if (dnaSequenceFromSaved[100] != dnaSequenceFromSaved2[100])
	{
		std::cout << "Sequences are different" << std::endl;
	}

	std::cout << "Execution time of functionA: " << timeA << " ms" << std::endl;


	std::cout << "Execution time of functionB: " << timeB << " ms" << std::endl;

	// Compare execution times
	if (timeA < timeB) {
		std::cout << "functionA is faster." << std::endl;
	}
	else if (timeA > timeB) {
		std::cout << "functionB is faster." << std::endl;
	}
	else {
		std::cout << "Both functions have the same execution time." << std::endl;
	}
}


void MatrixTest()
{
	string fastaFile = R"(C:\Braude\Projects\DNA\GCF_000001405.40_GRCh38.p14_genomic.txt)";
	std::cout << "(1)Loading of DNA Started from file : " << fastaFile << std::endl;

	string dnaSequenceFromSaved = load_previously_saved_data2(fastaFile);

	std::cout << "(1)DNA loaded! Size of sequence is  : " << dnaSequenceFromSaved.size() << std::endl;

	int word_size = 5;
	std::string firstSeq = dnaSequenceFromSaved.substr(0, 20);
	std::cout << "The First Sequence is  : " << firstSeq << std::endl;
	auto matrix1 = GenerateOccurrenceMatrix(firstSeq, word_size);
	std::cout << "Matrix1 : " << std::endl;
	displayMatrix(matrix1);

	std::string secondSeq = dnaSequenceFromSaved.substr(word_size, 20);
	std::cout << "The Second Sequence is  : " << secondSeq << std::endl;
	auto matrix2 = GenerateOccurrenceMatrix(secondSeq, word_size);
	std::cout << "Matrix2 : " << std::endl;
	displayMatrix(matrix2);

	std::string word = dnaSequenceFromSaved.substr(0, word_size);
	std::cout << "Word To Remove : " << word << std::endl;
	std::string seqToRemove = dnaSequenceFromSaved.substr(0, word_size);
	std::cout << "Section To Remove : " << seqToRemove << std::endl;
	auto matrixToSub = GenerateOccurrenceMatrix(seqToRemove, word_size);
	auto matrixAfterSub = RemoveSequenceFromOccurrenceMatrix(matrix1, seqToRemove, word_size);
	std::cout << "matrixAfterSub : " << std::endl;
	displayMatrix(matrixAfterSub);


	std::string wordToAdd = dnaSequenceFromSaved.substr(20, word_size);
	std::cout << "Word To Add : " << wordToAdd << std::endl;
	std::string seqToAdd = dnaSequenceFromSaved.substr(20, word_size);
	std::cout << "Section Added: " << seqToAdd << std::endl;
	auto matrixAfterAdd = AddSequenceToOccurrenceMatrix(matrixAfterSub, seqToAdd, word_size);
	std::cout << "matrixAfterAdd : " << std::endl;
	displayMatrix(matrixAfterAdd);

}

void SegmentTest(std::string filePath, int minSegmentSize, int wordSize, int lookaheadSize)
{
	//string fastaFile = R"(C:\Braude\Projects\DNA\GCF_000001405.40_GRCh38.p14_genomic.txt)";
	std::cout << "Loading of DNA Started from file : " << filePath << std::endl;

	string dnaSequenceFromSaved = load_previously_saved_data2(filePath);

	std::cout << "DNA loaded! Size of sequence is  : " << dnaSequenceFromSaved.size() << std::endl;

	//int minSegmentSize = 2000;
	//int word_size = 5;
	//int lookaheadSize = 20000;
	//int seqSize = 4000000;
	//std::string firstSeq = dnaSequenceFromSaved.substr(0, seqSize);
	std::cout << "The Sequence size is  : " << dnaSequenceFromSaved.size() << std::endl;
	std::cout << "The Word Size is  : " << wordSize << std::endl;
	std::cout << "The Minimum Segment Size is  : " << minSegmentSize * wordSize << " nucleotides" << std::endl;
	std::cout << "The lookahead Size is  : " << lookaheadSize * wordSize << " nucleotides" << std::endl;

	auto segments = SegmentDNACostAndWord(dnaSequenceFromSaved, minSegmentSize, wordSize, lookaheadSize);

	// Print the resulting segments, costs, and best words
	std::cout << "Segmented DNA sequence with costs and best words:\n";
	for (const auto& [start, end, cost, bestWord] : segments) {
		std::cout << "Segment: [" << start << ", " << end << ") -> "
			<< " len: " << end - start
			<< " | Cost: " << cost
			<< " | Best Word: " << bestWord << "\n";
	}

	std::string fileName = "segments_output_"
		+ std::to_string(minSegmentSize) + "_"
		+ std::to_string(wordSize) + "_"
		+ std::to_string(lookaheadSize) + ".csv";
	saveSegmentsToCSV(segments, fileName);

}

void DetectIsochores(std::string filePath, size_t window_size, double gc_threshold)
{
	std::cout << "Loading of DNA Started from file : " << filePath << std::endl;

	string dnaSequenceFromSaved = load_previously_saved_data2(filePath);

	std::cout << "DNA loaded! Size of sequence is  : " << dnaSequenceFromSaved.size() << std::endl;
	std::cout << "The Window Size is  : " << window_size << std::endl;
	std::cout << "The GC threshold is  : " << gc_threshold * 100 << " percents" << std::endl;
	auto isochores = detect_isochores(dnaSequenceFromSaved, window_size, gc_threshold);
	saveIsochoresToCsv(isochores, window_size, gc_threshold);
}

void DetectIsochores2(std::string filePath, size_t window_size, double gc_threshold)
{
	std::cout << "Loading of DNA Started from file : " << filePath << std::endl;

	string dnaSequenceFromSaved = load_previously_saved_data2(filePath);

	std::cout << "DNA loaded! Size of sequence is  : " << dnaSequenceFromSaved.size() << std::endl;
	std::cout << "The Window Size is  : " << window_size << std::endl;
	std::cout << "The GC threshold is  : " << gc_threshold * 100 << " percents" << std::endl;
	auto isochores = detect_isochores2(dnaSequenceFromSaved);
	saveIsochoresToCsv(isochores, 10000, 10000);
}

void processFASTATest(const std::string& inputFile, int windowSize, const std::string& outputFile)
{
	string sequence = load_previously_saved_data2(inputFile);

	if (sequence.empty()) {
		cerr << "Error: Loaded sequence is empty!" << endl;
	}

	// Remove FASTA headers (lines starting with '>')
	sequence.erase(remove_if(sequence.begin(), sequence.end(), [](char c) { return c == '>' || c == '\n'; }), sequence.end());

	// Process sequence and compute GC content
	processSequence(sequence, windowSize, outputFile);
}

void testZakhariaCode()
{
	std::string filename = R"(C:\Projects\DNA_Hidden_Repeat_Detection\x64\Debug\ncbi_dataset\ncbi_dataset\data\GCF_000001405.40\GCF_000001405.40_GRCh38.p14_genomic.fna)";

	GenomDNA genom;
	genom.Isoch1_TriplCorel(filename);
}

void MergeSegmentTest(std::string segmentFilePath, std::string sequanceFilePath, int minSegmentSize, int wordSize, int lookaheadSize)
{
	std::cout << "Loading of DNA Started from file : " << sequanceFilePath << std::endl;

	string dnaSequenceFromSaved = load_previously_saved_data2(sequanceFilePath);

	std::cout << "DNA loaded! Size of sequence is  : " << dnaSequenceFromSaved.size() << std::endl;

	std::cout << "The Word Size is  : " << wordSize << std::endl;
	std::cout << "The Minimum Segment Size is  : " << minSegmentSize * wordSize << " nucleotides" << std::endl;
	std::cout << "The lookahead Size is  : " << lookaheadSize * wordSize << " nucleotides" << std::endl;

	std::cout << "Loading of Segments started from file : " << segmentFilePath << std::endl;

	auto segments = loadSegmentsFromCSV(segmentFilePath);

	std::cout << "Segments loaded! Number of segments before merge is  : " << segments.size() << std::endl;

	std::cout << "Merge Segments started " << std::endl;
	auto merged = MergeSimilarSegments(segments, dnaSequenceFromSaved, wordSize);

	std::string fileName = "merged_segments_output_"
		+ std::to_string(minSegmentSize) + "_"
		+ std::to_string(wordSize) + "_"
		+ std::to_string(lookaheadSize) + ".csv";

	saveSegmentsToCSV(merged, fileName);

}