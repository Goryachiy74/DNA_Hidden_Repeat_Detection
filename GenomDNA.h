#ifndef GENOMDNA_H
#define GENOMDNA_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <array>

struct CorStruc1 {
  double t1 = 0.0, t2 = 0.0, t3 = 0.0;
};

struct CorType1 {
  long pos;
  int type;
};

class GenomDNA {
public:
  void Isoch1_TriplCorel(const std::string& filename);
};

#endif // GENOMDNA_H
