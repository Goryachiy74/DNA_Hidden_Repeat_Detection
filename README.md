
# DNA Hidden Repeat Detection

A new generation of DNA hidden repeat detection algorithm and its application for isochore research.

## 📖 Overview
This project implements a C++ algorithm to detect hidden repeats in human DNA and analyze segment similarity to identify isochores. The goal is to better understand genomic structure and composition through repeat detection and GC content analysis.

## 🚀 Features
- Detect hidden repeats in DNA sequences
- Merge segments based on similarity
- Analyze GC content to identify isochores
- Output detailed segment distribution and isochore charts

## 🛠️ Installation
1. Clone the repository:
```bash
git clone https://github.com/Goryachiy74/DNA_Hidden_Repeat_Detection.git
```
2. Open the solution file in **Visual Studio**:
   - `DNA_Hidden_Repeat_Detection.sln`
3. Build the project.

## 📌 Usage
1. Run the executable:
```bash
./DNA_Hidden_Repeat_Detection <input_file>
```
2. Output files will be generated in the `/output` folder:
   - `segments_output_*.txt` – Detected segments
   - `merged_segments_output_*.txt` – Merged segments
   - `isochores_output_*.txt` – Identified isochores

## 🧪 Example
Example input file:
```
ATGCGTACGATCGATCGTACGATCGTAGCTAGCTACGATCG
```
Example output:
- `segments_output_*.txt` → List of repeated segments
- `merged_segments_output_*.txt` – List Merged segments
- `isochores_output_*.txt` → GC content and isochore data

## 📊 Results
The program generates:
- Summary of detected segments
- GC content plots
- Visual analysis of segment distribution

## 📂 CSV Output Compatibility
The program generates a CSV file that can be used as input for the [DNA Chart App](https://github.com/Goryachiy74/dna_chart_app) to visualize segment data and isochore analysis.

## 🤝 Contributing
Contributions are welcome! Open an issue or submit a pull request.

## 📄 License
MIT License. See `LICENSE` for details.
