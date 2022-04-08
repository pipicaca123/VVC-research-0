#ifndef __MMLABSIGNALPROCESS__
#define __MMLABSIGNALPROCESS__

#include <opencv2/opencv.hpp>
#include <sys/stat.h>
#include <string>
#include <iostream>
#include "Picture.h"
#include "EncoderLib/EncModeCtrl.h"
class MMlab_SignalProcess{
private:
    std::vector<cv::Mat> inputs;
    cv::Mat processedCTU_img;
    int CTUx,CTUy,CTUh,CTUw;
    void CUStatistics(int posx, int posy, int width, int height);
    //debug function
    void ImageShow(cv::Mat& img); 
public:
    MMlab_SignalProcess(){  }; //do nothing
    void GetPixelInputs(const CodingStructure &cs); //get the orignal pixel values and transform to CV format.
    void ImageProcessing(void); // apply Canny algorithm
    void EarlyStopAlgorithm(EncModeCtrl& m_modeCtrl,int posx, int posy, int width, int height);
    void ReleaseprocessedCTUimg(){    inputs.pop_back();}
    
};
#endif