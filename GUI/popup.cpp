#include "mainwindow.h"
#include "QMessageBox"
#include "QFile"
#include <QDebug>

#include "plotwin.h"

#define NPLOT 100

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::setupPopup()
{
    connect(pushButton_SF_1,SIGNAL(clicked()),this,SLOT(pushButton_clicked()));
    connect(pushButton_SF_2,SIGNAL(clicked()),this,SLOT(pushButton_clicked()));
    connect(pushButton_glucose_0,SIGNAL(clicked()),this,SLOT(pushButton_clicked()));
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::pushButton_clicked()
{

    int HMI_SCALE = 1;
    int cellType;
    QString title, plotType, plotName;

    QObject *senderObj = sender(); // This will give Sender object
    QString senderObjName = senderObj->objectName();
//    qDebug() << senderObjName << "\n";
    plotwin = new PlotWin(this);
    QWidget *cw = plotwin->centralWidget();
    QFrame *plotFrame = cw->findChild<QFrame *>("plotFrame");
    QCustomPlot *popup_plot = new QCustomPlot(plotFrame);
    popup_plot->setObjectName("popup_plot");
    popup_plot->setGeometry(QRect(5*HMI_SCALE, 5*HMI_SCALE, 660*HMI_SCALE, 380*HMI_SCALE));

    // generate data:
    // Extract plot type and cell type from senderObjName
    QStringList list = senderObjName.split("_");
    if (list.size() != 3) return;
    plotType = list[1];
    cellType = list[2].toInt();
    if (plotType == "SF") {
        plotName = "Survival Fraction cell type " + list[2];
        plotwin->setWindowTitle(plotName);
        double C_O2;
        double maxdose = 50;
        C_O2 = 0;
        QVector<double> x0(NPLOT), y0(NPLOT); // initialize with entries 0..100
        makeSFPlot(list[2], C_O2, maxdose, &x0, &y0);
        // create graph and assign data to it:
        popup_plot->addGraph();
        popup_plot->graph(0)->setData(x0, y0);
        popup_plot->graph(0)->setPen(QPen(Qt::red));
        popup_plot->graph(0)->setName("O2 = 0%");

        C_O2 = 0.18;
        QVector<double> x1(NPLOT), y1(NPLOT); // initialize with entries 0..100
        makeSFPlot(list[2], C_O2, maxdose, &x1, &y1);
        // create graph and assign data to it:
        popup_plot->addGraph();
        popup_plot->graph(1)->setData(x1, y1);
        popup_plot->graph(1)->setPen(QPen(Qt::blue));
        popup_plot->graph(1)->setName("O2 = 20%");

        // give the axes some labels:
        popup_plot->xAxis->setLabel("Dose (Gy)");
        popup_plot->yAxis->setLabel("SF");
        // set axes ranges
        popup_plot->xAxis->setRange(0, maxdose);
//        popup_plot->xAxis->setAutoSubTicks(false);
//        popup_plot->xAxis->setSubTickCount(5);
//        popup_plot->xAxis->setAutoTickCount(6);
//        popup_plot->xAxis->setTickStep(5);
        popup_plot->yAxis->setRange(1.e-5, 1);
        popup_plot->yAxis->setScaleType(QCPAxis::stLogarithmic);
        popup_plot->yAxis->setScaleLogBase(10);
        popup_plot->yAxis->setNumberFormat("eb"); // e = exponential, b = beautiful decimal powers
        popup_plot->yAxis->setNumberPrecision(0); // makes sure "1*10^4" is displayed only as "10^4"
        popup_plot->legend->setVisible(true);
    } else if (plotType == "glucose") {
        plotName = "Medium Glucose Depletion";
        plotwin->setWindowTitle(plotName);
        double ndays;
        double C_O2;
        C_O2 = 0.18;
        QVector<double> x0(NPLOT), y0(NPLOT); // initialize with entries 0..100
        makeGlucosePlot(&ndays, &x0, &y0);
        // create graph and assign data to it:
        popup_plot->addGraph();
        popup_plot->graph(0)->setData(x0, y0);
        popup_plot->graph(0)->setPen(QPen(Qt::blue));
        // give the axes some labels:
        popup_plot->xAxis->setLabel("Day");
        popup_plot->yAxis->setLabel("Concentration (mM)");
        // set axes ranges
        popup_plot->xAxis->setRange(0, ndays);
        popup_plot->yAxis->setRange(0, 10);
//        popup_plot->yAxis->setScaleType(QCPAxis::stLogarithmic);
//        popup_plot->yAxis->setScaleLogBase(10);
//        popup_plot->yAxis->setNumberFormat("eb"); // e = exponential, b = beautiful decimal powers
//        popup_plot->yAxis->setNumberPrecision(0); // makes sure "1*10^4" is displayed only as "10^4"
//        popup_plot->legend->setVisible(true);
//        popup_plot->graph(0)->setName("O2 = 0%");

    }

    //    for (int i=0; i<101; ++i)
    //    {
    //      x[i] = i/50.0 - 1; // x goes from -1 to 1
    //      y[i] = x[i]*x[i];  // let's plot a quadratic function
    //    }

    plotwin->show();
}

//--------------------------------------------------------------------------------------------------------
// (NPLOT-1)*nt*dt = ndays*24*60*60 ==> number of time steps nt
//--------------------------------------------------------------------------------------------------------
void MainWindow::makeGlucosePlot(double *ndays, QVector<double> *x, QVector<double> *y)
{
    QLineEdit *line;
    double dt=100;  // sec
    double C, vol_cm3, t, metab, dCdt;
    double MM_C0, max_cell_rate;
    int ncells, Ng, nt;

    // Get parameter values from the GUI fields for glucose
    line = findChild<QLineEdit *>("lineEdit_glucose_ncells");
    ncells = line->text().toInt();
    line = findChild<QLineEdit *>("lineEdit_glucose_ndays");
    *ndays = line->text().toDouble();
    line = findChild<QLineEdit *>("line_MEDIUM_VOLUME");
    vol_cm3 = line->text().toDouble();
    line = findChild<QLineEdit *>("line_GLUCOSE_BDRY_CONC");
    C = line->text().toDouble();
    line = findChild<QLineEdit *>("line_GLUCOSE_MM_KM");
    MM_C0 = line->text().toDouble();
    line = findChild<QLineEdit *>("line_GLUCOSE_HILL_N");
    Ng = line->text().toInt();
    line = findChild<QLineEdit *>("line_GLUCOSE_CONSUMPTION");
    max_cell_rate = line->text().toDouble();

    max_cell_rate *= 1.0e6;     // mol/cell/s -> mumol/cell/s
    nt = ((*ndays)*24*60*60)/((NPLOT-1)*dt);
    (*x)[0] = 0;
    (*y)[0] = C;
    t = 0;
    for (int i=1; i<NPLOT; ++i)
    {
        for (int k=0; k<nt; k++) {
            metab = pow(C,Ng)/(pow(MM_C0,Ng) + pow(C,Ng));
//            qDebug("C: %f Ng: %d MM_C0: %f metab: %f",C,Ng,MM_C0,metab);
            dCdt = (-metab*max_cell_rate)/vol_cm3;	// convert mass rate (mol/s) to concentration rate (mM/s)
            dCdt = ncells*dCdt;
            C = C + dCdt*dt;
            if (C < 0) C = 0;
            t = t + dt;
        }
        (*x)[i] = t/(24*60*60);     // sec -> days
        (*y)[i] = C;
    }
}

//--------------------------------------------------------------------------------------------------------
// Assume that SER = 1 (no drug sensitisation)
//--------------------------------------------------------------------------------------------------------
void MainWindow::makeSFPlot(QString cellTypeStr, double C_O2, double maxdose, QVector<double> *x, QVector<double> *y)
{
    QLineEdit *line;
    double dose;

    // Get parameter values from the GUI fields for celltype
    QString objAlphaName = "line_RADIATION_ALPHA_H_" + cellTypeStr;
    QString objBetaName = "line_RADIATION_BETA_H_" + cellTypeStr;
    QString objOERAlphaName = "line_RADIATION_OER_ALPHA_" + cellTypeStr;
    QString objOERBetaName = "line_RADIATION_OER_BETA_" + cellTypeStr;
    QString objKmName = "line_RADIATION_KM_" + cellTypeStr;
    line = findChild<QLineEdit *>(objAlphaName);
    double LQ_alpha_H = line->text().toDouble();
    line = findChild<QLineEdit *>(objBetaName);
    double LQ_beta_H = line->text().toDouble();
    line = findChild<QLineEdit *>(objOERAlphaName);
    double LQ_OER_am = line->text().toDouble();
    line = findChild<QLineEdit *>(objOERBetaName);
    double LQ_OER_bm = line->text().toDouble();
    line = findChild<QLineEdit *>(objKmName);
    double LQ_K_ms = line->text().toDouble();

    double SER = 1;
    for (int i=0; i<NPLOT; ++i)
    {
        dose = (maxdose*i)/NPLOT;
        double OER_alpha_d = dose*(LQ_OER_am*C_O2 + LQ_K_ms)/(C_O2 + LQ_K_ms);
        double OER_beta_d = dose*(LQ_OER_bm*C_O2 + LQ_K_ms)/(C_O2 + LQ_K_ms);

        OER_alpha_d = OER_alpha_d*SER;
        OER_beta_d = OER_beta_d*SER;

        double expon = LQ_alpha_H*OER_alpha_d + LQ_beta_H*pow(OER_beta_d,2);
        double SF = exp(-expon);

        (*x)[i] = dose;
        (*y)[i] = SF;
    }

//    for (int i=0; i<101; ++i)
//    {
//      (*x)[i] = i/50.0 - 1; // x goes from -1 to 1
//      (*y)[i] = (*x)[i]*(*x)[i];  // let's plot a quadratic function
//    }

}

