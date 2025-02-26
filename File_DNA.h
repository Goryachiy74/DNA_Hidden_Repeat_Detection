#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <thread>
#include <mutex>

using namespace std;

///Function to parse a GenBank file and extract the sequence
string load_gen_bank_file(const string& filePath);

std::string load_fasta_file(const std::string& filename);

void save_loaded_data_as_file(const std::string& filename, const std::string& data);

std::string load_previously_saved_data(const std::string& filename);

std::string load_previously_saved_data2(const std::string& filename);

