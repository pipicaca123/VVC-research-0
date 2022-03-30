#include <chrono>
#include <string>
#include <iostream>
#include "mmlab_print.h"

MMLAB_printtool::MMLAB_printtool()
{
  start_time = std::chrono::high_resolution_clock::now();
}
void MMLAB_printtool::print(std::string s)
{
  std::cout<<s;
}
