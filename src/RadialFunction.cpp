
#include "RadialFunction.h"

#include <iostream>
#include <fstream>


int WriteToFile(const RadialFunction& rf, const std::string& filename)
{
  std::ofstream outfile;
  outfile.open(filename);

  for (size_t i=0; i<rf.radii.size(); ++i) {
    outfile << rf.radii[i] << "  " << rf.data[i] << "\n";
  }

  outfile.close();
  return 0;
}


