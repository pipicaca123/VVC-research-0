// #pragma comment(lib, "opencv_world455d.lib")
#include "mmlab_SignalProcess.h"
// #include <opencv2/core/mat.hpp>
#include <opencv2/opencv.hpp>
#include <string>
#include <iostream>

void MMlab_SignalProcess::GetPixelInputs(const CodingStructure &cs){
        //get the CS picture size & position
        int cs_x = cs.area.lx();
        int cs_y = cs.area.ly();
        int cs_w = cs.area.lwidth();
        int cs_h = cs.area.lheight();
        int stride = cs.picture->m_bufs[PIC_TRUE_ORIGINAL].bufs[0].stride;
        std::cout<<cs_x<<","<<cs_y<<","<<cs_w<<","<<cs_h<<std::endl;
        cv::Mat cs_img = cv::Mat::zeros(cs_w,cs_h,CV_16U);
        // std::cout<<"img size: "<<cs_img.size().width<<","<<cs_img.size().height<<std::endl;
        for(int y=cs_y;y<(cs_h+cs_y);y++){
            for(int x=cs_x;x<(cs_w+cs_x);x++){
                int imgx = x-cs_x;
                int imgy = y-cs_y;
                //in vtm the signal has been transform to 10 bits (guess: the program shift left 2 bits of orignal signal)
                cs_img.at<uint16_t>(imgx, imgy) = cs.picture->m_bufs[PIC_TRUE_ORIGINAL].bufs[0].buf[x * stride + y];
                cs_img.at<uint16_t>(imgx, imgy) >>= 2; //shift right 2 bits
                cs_img.at<uint16_t>(imgx, imgy) *= ((2<<8)/4); //exchange to 16 bits
            }
        }
        /* use for checking the matrix data is right.*/
        // std::cout << "M=" << std::endl << cs_img << std::endl; 
        cv::imshow("test",cs_img);
        cv::waitKey(0);
        cv::destroyAllWindows();
        
        /*-------------------------------------------*/
}
