#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSize.h"

#include "location.h"

typedef itk::Image<unsigned char,3> ImageType_u8;

#define V(a,b,c)  p_im[(c)*xysize+(b)*width+(a)]
#define Vin(a,b,c)  incell[(c)*xysize+(b)*width+(a)]

namespace Ui {
    class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void cell_tiffer();
    void cellFileSelecter();
    void tiffFileSelecter();
    void dataChanged();
    void checkReady();
private:
    Ui::MainWindow *ui;

public:
    FILE *fpout;
    int readCellData(const char *cellFile); //, LOCATION *locations, int *ncells);
    int readTiff(const char *, int *, int *, int *);
    int createTiffData();   //LOCATION *locations, int ncells);
    int createTiff(const char *, unsigned char *, int, int, int);
//    bool inSphere(LOCATION p);

    QString cellFileName;
    QString tiffFileName;
    bool is_celldata;
    bool is_tiff;
    bool cell_read;
    bool tiff_read;
    bool ready;
    bool is_network;
    ImageType_u8::Pointer im_u8;
    int width, height, depth, xysize;
    int margin;
    unsigned char *incell;
    unsigned char *p_im;
    double voxelsize;
    int ncells;
    double delta_x;
    LOCATION *location;
    double zfraction, zthickness;
    float sphereCentre[3], sphereRadius;
};

#endif // MAINWINDOW_H
