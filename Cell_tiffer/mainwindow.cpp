 //Qt
#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <QFileInfoList>
#include <QTextStream>
#include <QDir>
#include <iostream>
//ITK
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSize.h"

using namespace std;

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    ui->textEdit->setReadOnly(true);
    QString infoFile = QCoreApplication::applicationDirPath() + "/cell_tiffer_info.txt";
    QFile file(infoFile);
    bool ok = file.open(QIODevice::ReadOnly | QIODevice::Text);
    if (!ok) {
        ui->textEdit->append("The information file is missing:");
        ui->textEdit->append(infoFile);
    } else {
        QTextStream in(&file);
        QString line = in.readLine();
        while (!line.isNull()) {
            ui->textEdit->append(line);
            line = in.readLine();
        }
        file.close();
        ui->textEdit->moveCursor(QTextCursor::Start);
    }

    fpout = fopen("cell_tiffer.out","w");
    is_celldata = false;
    is_tiff = false;
    cell_read = false;
//    voxelsize = 1;
//    zfraction = 0.75;
//    zthickness = 10;    // um
    checkReady();
}

MainWindow::~MainWindow()
{
    delete ui;
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
void MainWindow::cellFileSelecter()
{
    ui->labelResult->setText("");
    cellFileName = QFileDialog::getOpenFileName(this,
        tr("Input cell location file"), ".", tr("Cell Files (*.out)"));
    if (cellFileName.compare(ui->labelCellFile->text()) != 0) {
        cell_read = false;
    }
    ui->labelCellFile->setText(cellFileName);
    is_celldata = true;
    checkReady();
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
void MainWindow::tiffFileSelecter()
{
    ui->labelResult->setText("");
    tiffFileName = QFileDialog::getSaveFileName(this,
        tr("Output TIFF file"), ".", tr("TIFF Files (*.tif)"));
    if (tiffFileName.compare(ui->labelTiffFile->text()) != 0) {
        printf("%s\n",tiffFileName.toAscii().constData());
        printf("%s\n",ui->labelTiffFile->text().toAscii().constData());
    }
    ui->labelTiffFile->setText(tiffFileName);
    is_tiff = true;
    checkReady();
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
void MainWindow::dataChanged()
{
    checkReady();
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
void MainWindow::checkReady()
{
    bool voxelOK, zfractionOK, zthicknessOK;
    voxelsize = ui->lineEdit_voxel->text().toDouble();
//    voxelsize[1] = ui->lineEdit_yvoxel->text().toFloat();
//    voxelsize[2] = ui->lineEdit_zvoxel->text().toFloat();
//    margin = ui->lineEditMargin->text().toFloat();

    if (voxelsize > 0) // && voxelsize[1] > 0 && voxelsize[2] > 0)
        voxelOK = true;
    else
        voxelOK = false;

    zfraction = ui->lineEdit_zfraction->text().toDouble();
    zfractionOK = (zfraction > -1 && zfraction < 1);
    zthickness = ui->lineEdit_zthickness->text().toDouble();
    zthicknessOK = (zthickness > 0);
    if (voxelOK && zfractionOK && zthicknessOK && is_celldata && is_tiff) {
        ready = true;
        ui->pushButtonGo->setEnabled(true);
    } else {
        ready = false;
        ui->pushButtonGo->setEnabled(false);
    }
}


//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
void MainWindow::cell_tiffer()
{
    int err;
    QString resultstr, numstr;

    ui->labelResult->setText("");
    if (!ready) {
        ui->labelResult->setText("Not ready: select files");
        return;
    }
    if (!cell_read) {
        resultstr = "reading cell location file...";
        ui->labelResult->setText(resultstr);
        QCoreApplication::processEvents();
        err = readCellData(cellFileName.toAscii().constData());
        if (err != 0) {
            resultstr = "FAILED: readCellData";
            ui->labelResult->setText(resultstr);
            return;
        }
    }
    resultstr = "creating image file...";
    ui->labelResult->setText(resultstr);
    QCoreApplication::processEvents();
    err = createTiffData();
    if (err != 0) {
        resultstr = "FAILED: createTiffData";
        ui->labelResult->setText(resultstr);
        return;
    }
    resultstr = "Created tiff dat...";
    ui->labelResult->setText(resultstr);
    fprintf(fpout,"Created tiff data\n");
    err = createTiff(tiffFileName.toAscii().constData(),p_im,width,height,depth);
    if (err != 0) {
        resultstr = "FAILED: createTiff";
        ui->labelResult->setText(resultstr);
        return;
    }
    if (p_im) free(p_im);
    checkReady();
    resultstr = "SUCCESS";
    ui->labelResult->setText(resultstr);
    return;
}

