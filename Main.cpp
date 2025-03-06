#include <iostream> // Required for std::cout and std::cerr
#include <string>   // Required for std::string
#include <cstdlib>  // Required for std::atoi
#include "File_DNA.h"
#include "Tests.h"

int main(int argc, char** argv)
{
 //   // Check if the correct number of arguments is provided
 //   if (argc != 5)
 //   {
 //       std::cerr << "Usage: " << argv[0] << " <file_path> <minSegmentSize> <wordSize> <lookaheadSize>" << std::endl;
 //       return 1; // Return an error code
 //   }

 //   // Parse the command-line arguments
 //   std::string file_path = argv[1]; // File path
 //   int minSegmentSize = std::atoi(argv[2]); // Convert to int
 //   int wordSize = std::atoi(argv[3]); // Convert to int
 //   int lookaheadSize = std::atoi(argv[4]); // Convert to int

 //   // Output the received parameters for verification
 //   std::cout << "File Path: " << file_path << std::endl;
 //   std::cout << "Minimum Segment Size: " << minSegmentSize << std::endl;
 //   std::cout << "Word Size: " << wordSize << std::endl;
 //   std::cout << "Lookahead Size: " << lookaheadSize << std::endl;

 //   // Call your functions here
 //   // CompareLoadSpeed();
 //   // MatrixTest();
 //   //DetectIsochores2(file_path, minSegmentSize, 0.2);
 //  // SegmentTest(file_path, minSegmentSize, wordSize, lookaheadSize);
 //   string inputFile = R"(C:\Braude\Projects\DNA_Hidden_Repeat_Detection\GCF_000001405.40_GRCh38.p14_genomic.fna)";
 //   string fastaFileToSave = R"(C:\Braude\Projects\DNA\GCF_000001405.40_GRCh38.p14_genomic.txt)";
 //   string outputFile = "gc_content500.csv";
 //   int windowSize = 500000; // 100kb window
	//processFASTATest(fastaFileToSave, windowSize, outputFile);

  testZakhariaCode();
    return 0; // Successful execution
}
