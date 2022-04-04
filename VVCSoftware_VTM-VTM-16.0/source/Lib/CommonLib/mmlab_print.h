#pragma once
#ifndef __MMLABPRINT__
#define __MMLABPRINT__

#include <chrono>
#include <string>
#define xcompressNO 0
#define xcheckRDcostintraNO 1

class MMLAB_printtool
{
private:
  std::chrono::high_resolution_clock::time_point start_time;
  int                                            enterNO[2] = { 0 };

public:
  MMLAB_printtool();
  void print(std::string s);
};
#endif