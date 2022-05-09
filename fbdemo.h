#ifndef FBDEMO_H
#define FBDEMO_H

#include <complex>
#include <iostream>
#include <valarray>

#include <QMainWindow>
#include <QString>

QT_BEGIN_NAMESPACE
namespace Ui { class FBDemo; }
QT_END_NAMESPACE

class FBDemo : public QMainWindow
{
    Q_OBJECT

public:
    FBDemo(QWidget *parent = nullptr);
    ~FBDemo();

private slots:
    void onClickedLoadBtn();
    void onClickedParallelBtn();
    void onClickedFanBtn();
    void onClickedSLFilterBtn();
    void onClickedRLFilterBtn();
    void onClickedReconBtn();

private:
    QString fileName;
    QImage img;

    int **pixelValue;
    int **reconImg;
    double **sinogram;
    double **fanSino;
    double **filtSino;

private:
    Ui::FBDemo *ui;
};
#endif // FBDEMO_H
