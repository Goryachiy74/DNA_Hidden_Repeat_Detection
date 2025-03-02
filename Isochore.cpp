#include "Isochore.h"



std::mutex mtxI; // Mutex for thread safety
uint64_t progressI = 0.0; // Shared progress variable
bool runningI = true; // Control variable for the progress thread
uint64_t totalsizeI = 100;

// Parameters
const int WINDOW_SIZE = 300000;  // 300 kb
const int STEP_SIZE = 100000;    // 100 kb

static void updateProgress()
{
    auto startTime = std::chrono::high_resolution_clock::now(); // Start time

    while (runningI)
    {
        std::this_thread::sleep_for(std::chrono::seconds(1)); // Update every 10 seconds

        // Lock the mutex to safely access the progress variable
        std::lock_guard<std::mutex> lock(mtxI);
        double percentage = (static_cast<double>(progressI) / totalsizeI) * 100;

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
double calculateGCContent(const std::string& sequence) {
    int gcCount = 0;
    for (char base : sequence) {
        if (base == 'G' || base == 'C') {
            gcCount++;
        }
    }
    return (sequence.length() > 0) ? (gcCount * 100.0 / sequence.length()) : 0;
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
            std::lock_guard<std::mutex> lock(mtxI);
            progressI = i; // Update progress
            totalsizeI = length - window_size;
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
    runningI = false;
    progressThread.join(); // Wait for the progress thread to finish
    return isochores;
}

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

// Function to detect isochores
std::vector<Isochore> detect_isochores2(const std::string& genomeSequence) {
    std::vector<Isochore> isochores;
    std::thread progressThread(updateProgress);

    for (size_t i = 0; i + WINDOW_SIZE <= genomeSequence.size(); i += STEP_SIZE) {
        std::string windowSeq = genomeSequence.substr(i, WINDOW_SIZE);
        double gcContent = calculateGCContent(windowSeq);

        // Update progress every 1000 iterations to reduce mutex locking
        if (i % 1000 == 0)
        {
            std::lock_guard<std::mutex> lock(mtxI);
            progressI = i; // Update progress
            totalsizeI = genomeSequence.size();
        }

        Isochore iso;
        iso.start = i;
        iso.end = i + WINDOW_SIZE;
        iso.gc_content = gcContent;

        isochores.push_back(iso);
    }
    // Stop the progress thread
    runningI = false;
    progressThread.join(); // Wait for the progress thread to finish
    return isochores;
}

// Function to calculate GC content in a DNA sequence
double calculateGCContent2(const std::string& sequence) {
    int gcCount = 0, totalCount = 0;
    for (char base : sequence) {
        if (base == 'G' || base == 'C' || base == 'g' || base == 'c') {
            gcCount++;
        }
        if (isalpha(base)) { // Ignore non-DNA characters
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
    totalsizeI = 0;

    // Compute total genome size
    while (getline(file, line)) {
        if (!line.empty() && line[0] != '>') { // Ignore headers
            totalsizeI += line.size();
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

        while (sequence.size() >= windowSize) {
            std::string window = sequence.substr(0, windowSize);
            double gcContent = calculateGCContent(window);

            outputFile << position << "," << gcContent << "\n";

            position += windowSize;

            // Update progress
            {
                std::lock_guard<std::mutex> lock(mtxI);
                progressI += windowSize;
            }

            sequence = sequence.substr(windowSize);
        }
    }

    // Stop progress thread
    runningI = false;
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
    totalsizeI = sequence.size(); // Set total size for progress tracking

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
            std::lock_guard<std::mutex> lock(mtxI);
            progressI = position;
        }
    }

    // Stop progress thread
    runningI = false;
    progressThread.join();

    outputFile.close();
    std::cout << "\nProcessing complete! Output saved in " << outputCSV << std::endl;
}