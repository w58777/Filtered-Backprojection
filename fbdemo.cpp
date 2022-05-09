#include "fbdemo.h"
#include "ui_fbdemo.h"

#include <iostream>

#include <QPushButton>
#include <QFileDialog>
#include <QString>
#include <QImage>
#include <QPixmap>
#include <QDebug>
#include <QColor>

#define PI acos(-1)

FBDemo::FBDemo(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::FBDemo)
{
    ui->setupUi(this);
    connect(ui->loadBtn,&QPushButton::clicked,this,&FBDemo::onClickedLoadBtn);
    connect(ui->parallelBtn,&QPushButton::clicked,this,&FBDemo::onClickedParallelBtn);
    connect(ui->fanBtn,&QPushButton::clicked,this,&FBDemo::onClickedParallelBtn);
    connect(ui->slFilterBtn,&QPushButton::clicked,this,&FBDemo::onClickedSLFilterBtn);
    connect(ui->rlFilterBtn,&QPushButton::clicked,this,&FBDemo::onClickedRLFilterBtn);
    connect(ui->reconBtn,&QPushButton::clicked,this,&FBDemo::onClickedReconBtn);
}

FBDemo::~FBDemo()
{
    delete []FBDemo::pixelValue;
    delete []FBDemo::reconImg;
    delete []FBDemo::sinogram;
    delete []FBDemo::filtSino;
    delete ui;
}

void FBDemo::onClickedLoadBtn()
{
    // Load image from computer
    QString arg("Img(*.jpg *.png *.bmp)");
    FBDemo::fileName = QFileDialog::getOpenFileName(this, tr("Load Image"),"",tr("Img(*.jpg *.png *.bmp)"), &arg);

    if(FBDemo::fileName == nullptr)return;
    FBDemo::img.load(FBDemo::fileName);

    // Get size of image then show it
    int w = FBDemo::img.width();
    int h = FBDemo::img.height();
    ui->line_width->setText(QString::number(w));
    ui->line_height->setText(QString::number(h));

    // Convert format to gray scale 8-bit
    FBDemo::img = FBDemo::img.convertToFormat(QImage::Format_Grayscale8);

    // Upsampling or downsampling to 512x512
    double **tmpRed = new double*[512];
    double **tmpGreen = new double*[512];
    double **tmpBlue = new double*[512];
    for(int i = 0; i < 512; ++i){
        tmpRed[i] = new double[512];
        tmpGreen[i] = new double[512];
        tmpBlue[i] = new double[512];
        for(int j = 0; j < 512; ++j){
            tmpRed[i][j] = qRed(FBDemo::img.pixel(i*w/512,j*h/512));
            tmpGreen[i][j] = qGreen(FBDemo::img.pixel(i*w/512,j*h/512));
            tmpBlue[i][j] = qBlue(FBDemo::img.pixel(i*w/512,j*h/512));
        }
    }

    // Create local QImage variable and set pixel value
    QImage adjImg(512,512,QImage::Format_Grayscale8);
    for(int i = 0; i < 512; ++i){
        for(int j = 0; j < 512; ++j){
            adjImg.setPixel(i,j,qRgb(tmpRed[i][j],tmpGreen[i][j],tmpBlue[i][j]));
        }
    }

    // Delete the pointer, release memory
    delete[] tmpRed;
    delete[] tmpGreen;
    delete[] tmpBlue;
    FBDemo::img = adjImg;

    // Initial 2D array to store pixel value
    FBDemo::pixelValue = new int*[512];
    for(int i = 0; i < 512; ++i){
        FBDemo::pixelValue[i] = new int[512];
        for(int j = 0; j < 512; ++j){
            FBDemo::pixelValue[i][j] = qGray(FBDemo::img.pixel(i,j));
        }
    }

    // Show the image
    QPixmap pixmap(QPixmap::fromImage(FBDemo::img));
    ui->label_img->setPixmap(pixmap);
}

void FBDemo::onClickedParallelBtn()
{
    // Initial sinogram 180x724x2
    FBDemo::sinogram = new double*[360]; // From 0 to 180 degrees, delta theta is 1 degree
    for(int i = 0; i < 360; ++i){
        FBDemo::sinogram[i] = new double[sqrt(2)*512];
        for(int j = 0; j < sqrt(2)*512; ++j){
            FBDemo::sinogram[i][j] = 0.0;
        }
    }

    // Projection process, every pixel project to the sinogram at every angle of the source
    for(double theta = 0; theta < 180; theta += 0.5){
        for(int i = 0; i < 512; ++i){
            for(int j = 0; j < 512; ++j){
                int sinopos =  361.5 - (255.5 - i)*cos(theta*PI/180) - (255.5 - j)*sin(theta*PI/180);
                FBDemo::sinogram[int(theta*2)][sinopos] += FBDemo::pixelValue[i][j];
            }
        }
    }

    // Local 2D array is created to be aimed at show the sinogram
    double **tmp = new double*[512];
    for(int i = 0; i < 512; ++i){
        tmp[i] = new double[512];
        for(int j = 0; j < 512; ++j){
            tmp[i][j] = FBDemo::sinogram[i*360/512][j*724/512]; // Resize
        }
    }

    int max = 0;
    for(int i = 0; i < 512; ++i){ // Find the max value
        for(int j = 0; j < 512; ++j){
            if(tmp[i][j]>max)max = tmp[i][j];
        }
    }
    for(int i = 0; i < 512; ++i){ // Normalize the gray level range to 0~255
        for(int j = 0; j < 512; ++j){
            tmp[i][j] = tmp[i][j]*255/max;
        }
    }

    // Create local QImage variable to show the sinogram
    QImage sinoImg(512,512,QImage::Format_Grayscale8);
    for(int i = 0; i < 512; ++i){
        for(int j = 0; j < 512; ++j){
            sinoImg.setPixel(i,j,qRgb(tmp[i][j]*11/32,tmp[i][j]*16/32,tmp[i][j]*5/32)); // Every rgb portion approximately match the qcolor transform standard
        }
    }

    delete []tmp;

    // Show the sinogram
    QPixmap pixmap(QPixmap::fromImage(sinoImg));
    ui->label_sino->setPixmap(pixmap);
}

void FBDemo::onClickedFanBtn()
{
    // Initial sinogram 180x724x2
    FBDemo::sinogram = new double*[360];
    FBDemo::fanSino = new double*[360];
    for(int i = 0; i < 360; ++i){
        FBDemo::sinogram[i] = new double[sqrt(2)*512];
        FBDemo::fanSino[i] = new double[sqrt(2)*512];
        for(int j = 0; j < sqrt(2)*512; ++j){
            FBDemo::sinogram[i][j] = 0.0;
            FBDemo::fanSino[i][j] = 0.0;
        }
    }

    // Projection process, every pixel project to the sinogram at every angle of the source
    // Assuming R equals to 256*sqrt(6), fan angle equals to 60 degrees
    // gamma is the angle between the raythrough (x, y) and the center line of the associatedfan
    double R = 256*sqrt(6);
    for(int beta = 0; beta < 360; ++beta){
        for(int i = 0; i < 512; ++i){
            for(int j = 0; j < 512; ++j){
                int gammapos =  ((beta%90) + atan((R*cos(beta)-i-255.5)/(R*sin(beta)-j-255.5))*180/PI) * 724 / 60;
                FBDemo::fanSino[beta][gammapos] += FBDemo::pixelValue[i][j];
                // Transform to parallel space
                // theta = gamma + beta; gamma = beta -arctan((R*cos(beta)-i-255.5)/(R*sin(beta)-j-255.5))
                int theta = 2 * (beta%90) -  atan((R*cos(beta)-i-255.5)/(R*sin(beta)-j-255.5))*180/PI;
                // r = R * sin(gamma)
                double r = sin(atan((R*cos(beta)-i-255.5)/(R*sin(beta)-j-255.5))) * R + 361;
                if(r < 0 || r >= 724){}// out of range
                FBDemo::sinogram[theta*2][int(r)] += FBDemo::pixelValue[i][j];
            }
        }
    }

    // Local 2D array is created to be aimed at show the sinogram
    double **tmp = new double*[512];
    for(int i = 0; i < 512; ++i){
        tmp[i] = new double[512];
        for(int j = 0; j < 512; ++j){
            tmp[i][j] = FBDemo::fanSino[i*360/512][j*724/512]; // Resize
        }
    }

    int max = 0;
    for(int i = 0; i < 512; ++i){ // Find the max value
        for(int j = 0; j < 512; ++j){
            if(tmp[i][j]>max)max = tmp[i][j];
        }
    }
    for(int i = 0; i < 512; ++i){ // Normalize the gray level range to 0~255
        for(int j = 0; j < 512; ++j){
            tmp[i][j] = tmp[i][j]*255/max;
        }
    }

    // Create local QImage variable to show the sinogram
    QImage sinoImg(512,512,QImage::Format_Grayscale8);
    for(int i = 0; i < 512; ++i){
        for(int j = 0; j < 512; ++j){
            sinoImg.setPixel(i,j,qRgb(tmp[i][j]*11/32,tmp[i][j]*16/32,tmp[i][j]*5/32)); // Every rgb portion approximately match the qcolor transform standard
        }
    }

    delete []tmp;

    // Show the sinogram
    QPixmap pixmap(QPixmap::fromImage(sinoImg));
    ui->label_sino->setPixmap(pixmap);
}

void FBDemo::onClickedSLFilterBtn()
{
    FBDemo::filtSino = new double*[360];
    for(int i = 0; i < 360; ++i){
        FBDemo::filtSino[i] = new double[sqrt(2)*512];
        for(int j = 0; j < sqrt(2)*512; ++j){
            FBDemo::filtSino[i][j] = 0.0;
        }
    }

    // S-L filter
    double *slfilter = new double[2*sqrt(2)*512];
    for(int i = 0; i < sqrt(2)*512-1; ++i){
        slfilter[722-i] = -2 / ((PI*PI)*(4*(i+1)*(i+1) - 1));
    }
    slfilter[723] = 2 / (PI*PI);
    for(int i = 0; i < sqrt(2)*512; ++i){
        slfilter[i+724] = -2 / ((PI*PI)*(4*(i+1)*(i+1) - 1));
    }

    // filt
    for(int theta = 0; theta < 360; ++theta){
        for(int i = 0; i < 724; ++i){
            for(int k = 0; k < 724; ++k){
                FBDemo::filtSino[theta][i] += FBDemo::sinogram[theta][k] * slfilter[i-k+723];
            }
        }
    }

    // Local 2D array is created to be aimed at show the filtered sinogram
    double **tmp = new double*[512];
    for(int i = 0; i < 512; ++i){
        tmp[i] = new double[512];
        for(int j = 0; j < 512; ++j){
            tmp[i][j] = FBDemo::filtSino[i*360/512][j*724/512]; // Resize
        }
    }

    int max = tmp[0][0];
    int min = tmp[0][0];
    for(int i = 0; i < 512; ++i){ // Find the max and min value
        for(int j = 0; j < 512; ++j){
            if(tmp[i][j]>max)max = tmp[i][j];
            if(tmp[i][j]<min)min = tmp[i][j];
        }
    }

    for(int i = 0; i < 512; ++i){ // Normalize the gray level range to 0~255
        for(int j = 0; j < 512; ++j){
            tmp[i][j] = (tmp[i][j]-min)*255/(max-min);
        }
    }

    // Create local QImage variable to show the sinogram
    QImage sinoImg(512,512,QImage::Format_Grayscale8);
    for(int i = 0; i < 512; ++i){
        for(int j = 0; j < 512; ++j){
            sinoImg.setPixel(i,j,qRgb(tmp[i][j]*11/32,tmp[i][j]*16/32,tmp[i][j]*5/32)); // Every rgb portion approximately match the qcolor transform standard
        }
    }

    delete []tmp;

    // Show the sinogram
    QPixmap pixmap(QPixmap::fromImage(sinoImg));
    ui->label_sino->setPixmap(pixmap);
}

void FBDemo::onClickedRLFilterBtn()
{
    FBDemo::filtSino = new double*[360];
    for(int i = 0; i < 360; ++i){
        FBDemo::filtSino[i] = new double[sqrt(2)*512];
        for(int j = 0; j < sqrt(2)*512; ++j){
            FBDemo::filtSino[i][j] = 0.0;
        }
    }

    // R-L filter
    double *rlfilter = new double[2*sqrt(2)*512];
    for(int i = 0; i < sqrt(2)*512-1; ++i){
        if(i % 2 == 0) rlfilter[722-i] = -4 / ((PI*PI)*(i+1)*(i+1));
        else rlfilter[722-i] = 0;
    }
    rlfilter[723] = 1;
    for(int i = 0; i < sqrt(2)*512; ++i){
        if(i % 2 == 0) rlfilter[i+724] = -4 / ((PI*PI)*(i+1)*(i+1));
        else rlfilter[i+724] = 0;
    }

    // filt
    for(int theta = 0; theta < 360; ++theta){
        for(int i = 0; i < 724; ++i){
            for(int k = 0; k < 724; ++k){
                FBDemo::filtSino[theta][i] += FBDemo::sinogram[theta][k] * rlfilter[i-k+723];
            }
        }
    }

    // Local 2D array is created to be aimed at show the filtered sinogram
    double **tmp = new double*[512];
    for(int i = 0; i < 512; ++i){
        tmp[i] = new double[512];
        for(int j = 0; j < 512; ++j){
            tmp[i][j] = FBDemo::filtSino[i*360/512][j*724/512]; // Resize
        }
    }

    int max = tmp[0][0];
    int min = tmp[0][0];
    for(int i = 0; i < 512; ++i){ // Find the max and min value
        for(int j = 0; j < 512; ++j){
            if(tmp[i][j]>max)max = tmp[i][j];
            if(tmp[i][j]<min)min = tmp[i][j];
        }
    }

    for(int i = 0; i < 512; ++i){ // Normalize the gray level range to 0~255
        for(int j = 0; j < 512; ++j){
            tmp[i][j] = (tmp[i][j]-min)*255/(max-min);
        }
    }

    // Create local QImage variable to show the sinogram
    QImage sinoImg(512,512,QImage::Format_Grayscale8);
    for(int i = 0; i < 512; ++i){
        for(int j = 0; j < 512; ++j){
            sinoImg.setPixel(i,j,qRgb(tmp[i][j]*11/32,tmp[i][j]*16/32,tmp[i][j]*5/32)); // Every rgb portion approximately match the qcolor transform standard
        }
    }

    delete []tmp;

    // Show the sinogram
    QPixmap pixmap(QPixmap::fromImage(sinoImg));
    ui->label_sino->setPixmap(pixmap);
}

void FBDemo::onClickedReconBtn()
{
    // Initialize the int array of reconstructed image
    FBDemo::reconImg = new int*[512];
    for(int i = 0; i < 512; ++i){
        FBDemo::reconImg[i] = new int[512];
        for(int j = 0; j < sqrt(2)*512; ++j){
            FBDemo::reconImg[i][j] = 0;
        }
    }

    // Reconstruction process, pixel driven backprojection, every pixel find the value of filtered sinogram at every angle
    for(int i = 0; i < 512; ++i){
        for(int j = 0; j < 512; ++j){
            for(double theta = 0; theta < 180; theta +=0.5){
                int reconpos =  361.5 - (255.5 - i)*cos(theta*PI/180) - (255.5 - j)*sin(theta*PI/180);
                FBDemo::reconImg[i][j] += FBDemo::filtSino[int(theta*2)][reconpos];
            }
        }
    }

    int max = FBDemo::reconImg[0][0];
    int min = FBDemo::reconImg[0][0];
    for(int i = 0; i < 512; ++i){ // Find the max and min value
        for(int j = 0; j < 512; ++j){
            if(FBDemo::reconImg[i][j]>max)max = FBDemo::reconImg[i][j];
            if(FBDemo::reconImg[i][j]<min)min = FBDemo::reconImg[i][j];
        }
    }

    for(int i = 0; i < 512; ++i){ // Normalize the gray level range to 0~255
        for(int j = 0; j < 512; ++j){
            FBDemo::reconImg[i][j] = (FBDemo::reconImg[i][j]-min)*255/(max-min);
        }
    }

    // Create local QImage variable and set the pixel value after reconstruction
    QImage reImg(512,512,QImage::Format_Grayscale8);
    for(int i = 0; i < 512; ++i){
        for(int j = 0; j < 512; ++j){
            reImg.setPixel(i,j,qRgb(FBDemo::reconImg[i][j]*11/32,FBDemo::reconImg[i][j]*16/32,FBDemo::reconImg[i][j]*5/32));
        }
    }

    // Show the reconstructed image
    QPixmap pixmap(QPixmap::fromImage(reImg));
    ui->label_recon->setPixmap(pixmap);
}
