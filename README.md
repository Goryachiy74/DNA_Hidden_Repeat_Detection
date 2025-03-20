
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

---

# 🧬 Biological Purpose of the Application

## 🌍 Biological Background
The application addresses several key problems in genomics and bioinformatics:

### ✅ 1. DNA Hidden Repeat Detection
- Detects hidden repeats within DNA sequences that are not apparent through direct alignment.
- Hidden repeats can reveal important biological signals, such as functional elements, regulatory sequences, or evolutionary patterns.

### ✅ 2. Segment Similarity Analysis
- Merges similar segments based on cyclic rotation and similarity.
- Similar segments often correspond to functionally significant regions like enhancers, promoters, or repeats.

### ✅ 3. Isochore Detection
- Isochores are long stretches of DNA with relatively homogeneous GC content.
- Isochore structures are linked to gene density, replication timing, and chromatin organization.

### ✅ 4. Overlap Analysis Between Isochores and Repeats
- Finds overlaps between detected repeats and isochores.
- This can reveal functionally significant regions (e.g., active regulatory regions, conserved elements).

---

## 💡 Why This Application Is Good

### 🚀 High Performance and Scalability
- Efficient handling of large genomic sequences using multithreading.
- Sliding window approach minimizes memory usage and improves processing time.

### 🔬 Biological Relevance
- Hidden repeats and isochores are biologically significant:
  - Hidden repeats can mark transposon activity, genomic rearrangements, and regulatory elements.
  - Isochores influence gene expression and chromatin structure.

### 🎯 Detection of Functional Elements
- Hidden repeats are often found near:
  - **Promoters** – Regions where transcription is initiated.
  - **Enhancers** – Regions that increase gene expression.
  - **Replication origins** – Regions where DNA replication starts.

### 🧬 Insights Into Genome Evolution
- Isochores reflect evolutionary pressures and genome organization:
  - GC-rich isochores are often gene-rich and correspond to euchromatin (actively transcribed regions).
  - AT-rich isochores are linked to heterochromatin and gene-poor regions.

### 🏥 Link to Disease and Medical Research
- Hidden repeats and GC content are linked to:
  - **Cancer** – Repeats and mutations in GC-rich regions are linked to genomic instability.
  - **Neurodegenerative diseases** – Triplet repeats are involved in diseases like Huntington's disease.

---

## 🔬 Interesting Research Questions

### 1. Are Hidden Repeats Conserved Across Species?
- By comparing hidden repeats across species, you could identify evolutionarily conserved regulatory regions.

### 2. Do Isochore Boundaries Overlap With Functional Elements?
- Isochore boundaries might align with transcription start sites, replication origins, or chromatin modification sites.

### 3. Role of GC Content in Genome Organization
- GC content is linked to DNA stability, mutation rates, and recombination.

### 4. Mutational Hotspots and Hidden Repeats
- Do hidden repeats align with known fragile sites or oncogenic hotspots?

### 5. Are Hidden Repeats Linked to Epigenetic Regulation?
- Certain repeats may recruit proteins involved in DNA methylation and chromatin remodeling.

### 6. Functional Role of Overlapping Repeats and Isochores
- Are repeats preferentially located within isochores?
- Are repeat-rich isochores associated with higher evolutionary rates?

---

## 🚀 Future Directions

### 🔹 Functional Annotation
- Link detected repeats and isochores to known genes and regulatory elements.

### 🔹 Comparative Genomics
- Compare repeat and isochore structures across species.

### 🔹 Machine Learning
- Train a model to classify functional elements based on repeat and isochore data.

### 🔹 Clinical Applications
- Use repeat and isochore patterns to identify cancer biomarkers.

---

## 🌟 What Makes It Powerful
✅ Efficient handling of large genome files.  
✅ Combination of repeat detection and isochore analysis.  
✅ Parallelized for high performance.  
✅ Potential for discovery of novel functional elements.  
✅ Practical for comparative genomics and medical research.  

---

## 🤝 Contributing
Contributions are welcome! Open an issue or submit a pull request.

## 📄 License
MIT License. See `LICENSE` for details.

---

## 🌍 Author
**[Maxim Goryachev]** – [GitHub Profile](https://github.com/Goryachiy74)  
📧 **Contact:** goryachiy74@gmail.com  

---