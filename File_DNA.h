#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <thread>
#include <mutex>

using namespace std;

/// <summary>
/// Loads a DNA sequence from a GenBank file.
/// </summary>
/// <param name="filePath">Path to the GenBank file.</param>
/// <returns>The extracted DNA sequence as a string.</returns>
string load_gen_bank_file(const string& filePath);

/// <summary>
/// Loads a DNA sequence from a FASTA file.
/// </summary>
/// <param name="filename">Path to the FASTA file.</param>
/// <returns>The extracted DNA sequence as a string.</returns>
std::string load_fasta_file(const std::string& filename);

/// <summary>
/// Saves DNA sequence data to a file.
/// </summary>
/// <param name="filename">Path to the output file.</param>
/// <param name="data">DNA sequence to save.</param>
void save_loaded_data_as_file(const std::string& filename, const std::string& data);

/// <summary>
/// Loads previously saved DNA sequence from a file.
/// </summary>
/// <param name="filename">Path to the file containing saved DNA sequence.</param>
/// <returns>DNA sequence as a string.</returns>
std::string load_previously_saved_data(const std::string& filename);

/// <summary>
/// Loads previously saved DNA sequence from a file using binary mode.
/// </summary>
/// <param name="filename">Path to the file.</param>
/// <returns>DNA sequence as a string.</returns>
std::string load_previously_saved_data_binary_mode(const std::string& filename);

/// <summary>
/// Extracts all chromosomes from a FASTA file and saves them as individual files.
/// </summary>
/// <param name="filename">Path to the input file.</param>
void extract_all_chromosomes(const std::string& filename);

/// <summary>
/// Reads a specific chromosome from a FASTA file.
/// </summary>
/// <param name="filename">Path to the FASTA file.</param>
/// <param name="chromosome">Name of the chromosome to extract.</param>
/// <returns>DNA sequence of the chromosome as a string.</returns>
std::string read_chromosome(const std::string& filename, const std::string& chromosome);

/// <summary>
/// Reads a chromosome from a file.
/// </summary>
/// <param name="filename">Path to the file.</param>
/// <returns>DNA sequence as a string.</returns>
std::string read_chromosome_file(const std::string& filename);



