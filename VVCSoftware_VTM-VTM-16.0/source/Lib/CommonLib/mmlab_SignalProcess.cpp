// #pragma comment(lib, "opencv_world455d.lib")
#include "mmlab_SignalProcess.h"
// #include <opencv2/core/mat.hpp>
#include <opencv2/opencv.hpp>
#include <string>
#include <iostream>

//debugger parameter
#define COMPUTER_USE 0 // 1: I'm using lab's computer, 0: I'm using my own laptop
#define IMAGE_DEBUG // use for activate the image show bebugger

void MMlab_SignalProcess::GetPixelInputs(const CodingStructure &cs){
        //get the CS picture size & position
        CTUx = cs.area.lx(); //the definition is col->x, row->y
        CTUy = cs.area.ly();
        CTUw = cs.area.lwidth();
        CTUh = cs.area.lheight();
        int stride = cs.picture->m_bufs[PIC_TRUE_ORIGINAL].bufs[0].stride;
        std::cout<<CTUx<<","<<CTUy<<","<<CTUw<<","<<CTUh<<","<<stride<<std::endl;
        cv::Mat cs_img = cv::Mat::zeros(CTUh,CTUw,CV_16U);
        // std::cout<<"img size: "<<cs_img.size().width<<","<<cs_img.size().height<<std::endl;
        for(int y=CTUx;y<(CTUw+CTUx);y++){
            for(int x=CTUy;x<(CTUh+CTUy);x++){
                int imgy = x-CTUy;
                int imgx = y-CTUx;
               // std::cout << imgx << "," << imgy << std::endl;
                //in vtm the signal has been transform to 10 bits (guess: the program shift left 2 bits of orignal signal)
                cs_img.at<uint16_t>(imgy, imgx) = cs.picture->m_bufs[PIC_TRUE_ORIGINAL].bufs[0].buf[x * stride + y];
                cs_img.at<uint16_t>(imgy, imgx) >>= 2; //shift right 2 bits
                cs_img.at<uint16_t>(imgy, imgx) *= ((2<<8)/4); //exchange to 16 bits
            }
        }
        #ifdef IMAGE_DEBUG
            // ImageShow(img);
        #endif
        inputs.push_back(cs_img);
}


void MMlab_SignalProcess::ImageProcessing(void){
    cv::Mat img = inputs.back();
    // cvCanny only supports single-channel 8-bit images
    img.convertTo(img,CV_8UC1,1/256.0);
    // std::cout<<"image size is: "<<img.size().height<<","<<img.size().width<<std::endl;
    cv::GaussianBlur(img,processedCTU_img,cv::Size(3,3),1.5);
    cv::Canny(img,processedCTU_img,5,45);
    // for(int x=0;x<128;x++){
    //     for(int y=0;y<128;y++){  
    //         std::cout<<(int)(processedCTU_img.at<uchar>(x,y)/255)<<" ";
    //     }
    //     std::cout<<std::endl;
    // }
    #ifdef IMAGE_DEBUG
        ImageShow(processedCTU_img);
    #endif
    // inputs.pop_back();
}


void MMlab_SignalProcess::EarlyStopAlgorithm(EncModeCtrl& m_modeCtrl,int posx, int posy, int width, int height){
    // m_modeCtrl.printMode();
  CUStatistics(posx, posy, width, height);
}

void MMlab_SignalProcess::CUStatistics(int posx, int posy, int width, int height){
    std::vector<int> ver_edge(width,0);
    std::vector<int> hor_edge(height,0);
    // std::cout<<posx-CTUx<<","<<posy-CTUy<<","<<width<<","<<height<<std::endl;
    cv::Rect rect(posx-CTUx,posy-CTUy,width,height); // 4 parameters corresponse to x,y,w,h
    cv::Mat CU_img = cv::Mat(processedCTU_img,rect).clone();
    // if(width>=64 && height>=32)
    //     ImageShow(CU_img);

    //get horizontal edge statistics information
    for(int y=0;y<height;y++){
        for(int x=0;x<width;x++){
            hor_edge[y] += (int)(CU_img.at<uchar>(x,y)/255);
        }
        std::cout<<hor_edge[y]<<",";
    }
    std::cout<<std::endl;
    for(int x=0;x<width;x++){
        for(int y=0;y<height;y++){
            ver_edge[x] += (int)(CU_img.at<uchar>(x,y)/255);
        }
        std::cout<<ver_edge[x]<<",";
    }
    std::cout<<std::endl;
}


void MMlab_SignalProcess::ImageShow(cv::Mat& img){ // use for checking the matrix data is right or not.
    // std::cout << "M=" << std::endl << img << std::endl; 
    cv::imshow("test",img);
    #if (COMPUTER_USE) // lab computer
        cv::waitKey(10);
    #else//my own computer is ok for the below format.
        cv::waitKey(0);
        cv::destroyAllWindows();
    #endif
    
}