// #pragma comment(lib, "opencv_world455d.lib")
#include "mmlab_SignalProcess.h"
// #include <opencv2/core/mat.hpp>
#include <opencv2/opencv.hpp>
#include <string>
#include <iostream>

void MMlab_SignalProcess::GetPixelInputs(const CodingStructure &cs){
        //get the CS picture size & position
        int cs_x = cs.area.lx(); //the definition is col->x, row->y
        int cs_y = cs.area.ly();
        int cs_w = cs.area.lwidth();
        int cs_h = cs.area.lheight();
        int stride = cs.picture->m_bufs[PIC_TRUE_ORIGINAL].bufs[0].stride;
        std::cout<<cs_x<<","<<cs_y<<","<<cs_w<<","<<cs_h<<","<<stride<<std::endl;
        cv::Mat cs_img = cv::Mat::zeros(cs_h,cs_w,CV_16U);
        // std::cout<<"img size: "<<cs_img.size().width<<","<<cs_img.size().height<<std::endl;
        for(int y=cs_x;y<(cs_w+cs_x);y++){
            for(int x=cs_y;x<(cs_h+cs_y);x++){
                int imgy = x-cs_y;
                int imgx = y-cs_x;
               // std::cout << imgx << "," << imgy << std::endl;
                //in vtm the signal has been transform to 10 bits (guess: the program shift left 2 bits of orignal signal)
                cs_img.at<uint16_t>(imgy, imgx) = cs.picture->m_bufs[PIC_TRUE_ORIGINAL].bufs[0].buf[x * stride + y];
                cs_img.at<uint16_t>(imgy, imgx) >>= 2; //shift right 2 bits
                cs_img.at<uint16_t>(imgy, imgx) *= ((2<<8)/4); //exchange to 16 bits
            }
        }
        // ImageShow(cs_img);
        inputs.push_back(cs_img);
}


void MMlab_SignalProcess::ImageProcessing(void){
    cv::Mat img = inputs.back();
    // cvCanny only supports single-channel 8-bit images
    img.convertTo(img,CV_8UC1,1/256.0);
    // std::cout<<"image size is: "<<img.size().height<<","<<img.size().width<<std::endl;
    cv::GaussianBlur(img,img,cv::Size(3,3),1.5);
    cv::Canny(img,img,5,45);
    ImageShow(img);
    inputs.pop_back();
}

void MMlab_SignalProcess::ImageShow(cv::Mat& img){ // use for checking the matrix data is right or not.
    // std::cout << "M=" << std::endl << img << std::endl; 
    cv::imshow("test",img);
    cv::waitKey(0);
    cv::destroyAllWindows();
    
}