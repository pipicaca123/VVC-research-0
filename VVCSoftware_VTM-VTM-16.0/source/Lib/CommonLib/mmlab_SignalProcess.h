#ifndef __MMLABSIGNALPROCESS__
#define __MMLABSIGNALPROCESS__

#include <opencv2/core/mat.hpp>
#include <sys/stat.h>
#include <string>
#include <iostream>
#include "Picture.h"

class MMlab_SignalProcess{
private:
    static std::vector<cv::Mat> inputs;

public:
    MMlab_SignalProcess(){  }; //do nothing
    //get the orignal pixel values and transform to CV format.
    static void GetPixelInputs(const CodingStructure &cs);
};
#endif