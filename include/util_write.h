//#ifndef UTIL_WRITE_H
//#define UTIL_WRITE_H

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iomanip>

using namespace std;

template<class InputIterator, class InputIterator2>
void writedat(string filename,
              InputIterator begin1, InputIterator end1,
              InputIterator2 begin2, InputIterator2 end2,
              int precision1=10, int precision2=10){
  std::ofstream f;
  f.exceptions(std::ofstream::failbit | std::ofstream::badbit);
  f.open(filename);
  for ( ; begin1 != end1 and begin2 != end2; ++begin1, ++begin2)
    f << std::setprecision(precision1) << *begin1 << '\t'
      << std::setprecision(precision2) << *begin2 << '\n';
}

template<class InputIterator, class InputIterator2, class InputIterator3>
void writedat(string filename,
              InputIterator begin1, InputIterator end1,
              InputIterator2 begin2, InputIterator2 end2,
              InputIterator3 begin3, InputIterator3 end3,
              int precision1=10, int precision2=10, int precision3=10){
  std::ofstream f;
  f.exceptions(std::ofstream::failbit | std::ofstream::badbit);
  f.open(filename);
  for ( ; begin1 != end1 and begin2 != end2 and begin3 != end3; ++begin1, ++begin2, ++begin3)
    f << std::setprecision(precision1) << *begin1 << '\t'
      << std::setprecision(precision2) << *begin2 << '\t'
      << std::setprecision(precision3) << *begin3 << '\n';
}

template<class InputIterator, class InputIterator2, class InputIterator3, class InputIterator4>
void writedat(string filename,
              InputIterator begin1, InputIterator end1,
              InputIterator2 begin2, InputIterator2 end2,
              InputIterator3 begin3, InputIterator3 end3,
              InputIterator4 begin4, InputIterator4 end4,
              int precision1=10, int precision2=10, int precision3=10, int precision4=10){
  std::ofstream f;
  f.exceptions(std::ofstream::failbit | std::ofstream::badbit);
  f.open(filename);
  for ( ; begin1 != end1 and begin2 != end2 and begin3 != end3 and begin4 != end4; ++begin1, ++begin2, ++begin3, ++begin4)
    f << std::setprecision(precision1) << *begin1 << '\t'
      << std::setprecision(precision2) << *begin2 << '\t'
      << std::setprecision(precision3) << *begin3 << '\t'
      << std::setprecision(precision4) << *begin4 << '\n';
}

template<class InputIterator, class InputIterator2, class InputIterator3, class InputIterator4, class InputIterator5>
void writedat(string filename,
              InputIterator begin1, InputIterator end1,
              InputIterator2 begin2, InputIterator2 end2,
              InputIterator3 begin3, InputIterator3 end3,
              InputIterator4 begin4, InputIterator4 end4,
              InputIterator5 begin5, InputIterator5 end5,
              int precision1=10, int precision2=10, int precision3=10, int precision4=10, int precision5=10){
  std::ofstream f;
  f.exceptions(std::ofstream::failbit | std::ofstream::badbit);
  f.open(filename);
  for ( ; begin1 != end1 and begin2 != end2 and begin3 != end3 and begin4 != end4 and begin5 != end5; ++begin1, ++begin2, ++begin3, ++begin4, ++begin5)
    f << std::setprecision(precision1) << *begin1 << '\t'
      << std::setprecision(precision2) << *begin2 << '\t'
      << std::setprecision(precision3) << *begin3 << '\t'
      << std::setprecision(precision4) << *begin4 << '\t'
      << std::setprecision(precision5) << *begin5 << '\n';
}


template<class InputIterator, class InputIterator2, class InputIterator3, class InputIterator4, class InputIterator5, class InputIterator6>
void writedat(string filename,
              InputIterator begin1, InputIterator end1,
              InputIterator2 begin2, InputIterator2 end2,
              InputIterator3 begin3, InputIterator3 end3,
              InputIterator4 begin4, InputIterator4 end4,
              InputIterator5 begin5, InputIterator5 end5,
              InputIterator6 begin6, InputIterator6 end6,
              int precision1=10, int precision2=10, int precision3=10, int precision4=10, int precision5=10, int precision6=10){
  std::ofstream f;
  f.exceptions(std::ofstream::failbit | std::ofstream::badbit);
  f.open(filename);
  for ( ; begin1 != end1 and begin2 != end2 and begin3 != end3 and begin4 != end4 and begin5 != end5 and begin6 != end6; ++begin1, ++begin2, ++begin3, ++begin4, ++begin5, ++begin6)
    f << std::setprecision(precision1) << *begin1 << '\t'
      << std::setprecision(precision2) << *begin2 << '\t'
      << std::setprecision(precision3) << *begin3 << '\t'
      << std::setprecision(precision4) << *begin4 << '\t'
      << std::setprecision(precision5) << *begin5 << '\t'
      << std::setprecision(precision6) << *begin6 << '\n';
}

//#endif
