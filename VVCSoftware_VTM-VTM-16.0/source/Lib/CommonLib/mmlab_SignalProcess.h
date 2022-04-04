#ifndef __MMLABSIGNALPROCESS__
#define __MMLABSIGNALPROCESS__

#include <opencv2/opencv.hpp>
#include <sys/stat.h>
#include <string>
#include <iostream>
#include "Picture.h"

class MMlab_SignalProcess{
private:
    std::vector<cv::Mat> inputs;
    void ImageShow(cv::Mat& img); 
public:
    MMlab_SignalProcess(){  }; //do nothing
    void GetPixelInputs(const CodingStructure &cs); //get the orignal pixel values and transform to CV format.
    void ImageProcessing(void);                       // apply Canny algorithm
};
#endif