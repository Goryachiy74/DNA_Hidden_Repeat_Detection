#include "GenomDNA.h"

void GenomDNA::Isoch1_TriplCorel(const std::string& filename) {
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << '\n';
        return;
    }

    std::ofstream fout;
    std::string descr, trip1 = "AA", trip2 = "CG", sequence;
    constexpr int WINDOW_SIZE = 100000, MAX_DISTANCE = 500;
    std::vector<CorType1> PosCorMas1;
    
    // Initialize correlation matrices
    std::array<std::array<CorStruc1, MAX_DISTANCE>, 5> Corel1{}, Corel2{};
    
    char a;
    while (fin >> a) {
        if (a == '>') {  // Read description line
            std::getline(fin, descr);
            continue;
        }

        sequence.push_back(a);
        if (sequence.size() < trip1.size()) continue;

        // Shift sequence window
        if (sequence.size() > trip1.size()) sequence.erase(sequence.begin());

        // Check triplets and record positions
        CorType1 ct;
        if (sequence == trip1) ct = {static_cast<long>(sequence.size()), 1};
        else if (sequence == trip2) ct = {static_cast<long>(sequence.size()), 2};
        else continue;

        PosCorMas1.push_back(ct);
    }
    
    fin.close();

    // Compute correlations
    for (int i = 0; i < PosCorMas1.size() - 1; ++i) {
        for (int j = i + 1; j < PosCorMas1.size() && (PosCorMas1[j].pos - PosCorMas1[i].pos) < MAX_DISTANCE; ++j) {
            int dist = PosCorMas1[j].pos - PosCorMas1[i].pos;
            int n = std::min(4, dist / 10);

            if (PosCorMas1[i].type == 1 && PosCorMas1[j].type == 1) Corel1[n][dist].t1++;
            if (PosCorMas1[i].type == 2 && PosCorMas1[j].type == 2) Corel1[n][dist].t2++;
            if (PosCorMas1[i].type != PosCorMas1[j].type) Corel1[n][dist].t3++;
        }
    }

    // Smooth correlations
    for (int j = 0; j < 5; j++) {
        for (int i = 1; i < MAX_DISTANCE - 1; i++) {
            Corel2[j][i].t1 = (Corel1[j][i - 1].t1 + Corel1[j][i].t1 + Corel1[j][i + 1].t1) / 3.0;
            Corel2[j][i].t2 = (Corel1[j][i - 1].t2 + Corel1[j][i].t2 + Corel1[j][i + 1].t2) / 3.0;
            Corel2[j][i].t3 = (Corel1[j][i - 1].t3 + Corel1[j][i].t3 + Corel1[j][i + 1].t3) / 3.0;
        }
    }

    // Save results
    std::string output_file = "Corel_Isoch_" + trip1 + "_" + trip2 + ".dat";
    fout.open(output_file);
    if (!fout) {
        std::cerr << "Error: Could not create " << output_file << '\n';
        return;
    }

    fout << "num\tAA0\tCG0\tcros0\n";
    for (int i = 1; i < MAX_DISTANCE; i++) {
        fout << i;
        for (int j = 0; j < 5; j++) {
            fout << '\t' << Corel1[j][i].t1 << '\t' << Corel1[j][i].t2 << '\t' << Corel1[j][i].t3;
        }
        fout << '\n';
    }

    fout.close();
    std::cout << "Processing complete. Output saved in " << output_file << '\n';
}
