
#ifndef TFDH_UTILS_H
#define TFDH_UTILS_H

#include <ostream>
#include <vector>


template <typename T>
std::ostream& operator<<(std::ostream& s, const std::vector<T>& v) {
  // do nothing if v is empty
  if (v.size()==1) {
    s << v[0];
  }
  else if (v.size() > 1) {
    s << "(";
    for (size_t i=0; i<v.size(); ++i) {
      s << v[i];
      if (i != v.size()-1)
        s << ", ";
    }
    s << ")";
  }
  return s;
}


#endif // TFDH_UTILS_H
