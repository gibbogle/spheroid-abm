#ifndef FIELD_H
#define FIELD_H

#include <QDialog>
#include <QMessageBox>
#include <QtGui>
#include "plot.h"

#define CANVAS_WIDTH 696
#define MAX_CONC 6

struct field_data {
    int site[3];
    int state;
    double dVdt;
    double volume;
    double conc[MAX_CONC+1];
};

typedef field_data FIELD_DATA;

#define X_AXIS 1
#define Y_AXIS 2
#define Z_AXIS 3

#define OXYGEN 0
#define GLUCOSE 1
#define DRUG_A 2
#define DRUG_A_METAB 3
#define DRUG_B 4
#define DRUG_B_METAB 5
#define GROWTH_RATE 6       // we pretend that this is a concentration

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

extern "C" {
    void get_fieldinfo(int *, int *, double *, int *, int *, int *, int *);
    void get_fielddata(int *, double *, int *, int *, FIELD_DATA *, int *);
//    void get_concdata(int *, double *, double *);
}

class Field : public QMainWindow
{
public:
    Field(QWidget *);
    ~Field();
    void chooseParameters();
    void displayField(int, int *);
    void displayField1();
    void setSliceChanged();
    void chooseFieldColor(double c, double cmin, double cmax, bool use_log, int rgbcol[]);
    void chooseRateColor(double fr, int rgbcol[]);
    void getTitle(QString *);
    bool isConcPlot();
    void setConcPlot(bool);
    void makeConcPlot(QMdiArea *);
    void updateConcPlot();
    bool isVolPlot();
    void setVolPlot(bool);
    void makeVolPlot(QMdiArea *);
    void updateVolPlot();
    bool isOxyPlot();
    void setOxyPlot(bool);
    void makeOxyPlot(QMdiArea *);
    void updateOxyPlot();
    void selectConstituent();
    void setExecuting(bool);
    void setSaveImages(bool);
    void setUseLogScale(bool);

    QWidget *field_page;
    bool save_images;
    bool use_log;
    QGraphicsView* view;
    int NX;
    int axis;
    double fraction;
    int hour;
    int ifield;
    int nsites, nconst, const_used[16];
    QString const_name[16];
    QString constituentText;
    int constituent;
    bool slice_changed;
	bool constituent_changed;
    bool useConcPlot;
    bool useVolPlot;
    bool useOxyPlot;
    FIELD_DATA *data;
    Plot *pGconc, *pGvol, *pGoxy;
    bool executing;
    char msg[1024];

//public slots:
    void setConstituent(QAbstractButton* button);
    void setPlane(QAbstractButton* button);
	void setFraction(QString text);
};

#endif // FIELD_H
