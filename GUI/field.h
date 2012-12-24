#ifndef FIELD_H
#define FIELD_H

#include <QDialog>
#include <QMessageBox>
#include <QtGui>
#include "plot.h"

#define CANVAS_WIDTH 620

struct field_data {
    int site[3];
    int state;
    double volume;
    double conc[4];
};

typedef field_data FIELD_DATA;

#define X_AXIS 1
#define Y_AXIS 2
#define Z_AXIS 3

#define OXYGEN 0
#define GLUCOSE 1
#define DRUG_A 2
#define DRUG_B 3

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

extern "C" {
    void get_fieldinfo(int *, int *, double *, int *, int *, int *);
    void get_fielddata(int *, double *, int *, int *, FIELD_DATA *);
//    void get_concdata(int *, double *, double *);
}

class Field : public QMainWindow
{
public:
    Field(QWidget *);
    ~Field();
    void chooseParameters();
	void displayField();
    void setSliceChanged();
    void chooseColor(double fr, int rgbcol[]);
    void makeConcPlot(QMdiArea *);
    void updateConcPlot();

    QWidget *field_page;
    int NX;
    int axis;
    double fraction;
    int MAX_CHEMO;
    int nsites, nconst, cused[10];
    int constituent;
    bool slice_changed;
	bool constituent_changed;
	FIELD_DATA *data;
    Plot *pG;
    char msg[1024];

//public slots:
    void setConstituent(QAbstractButton* button);
    void setPlane(QAbstractButton* button);
	void setFraction(QString text);
};

#endif // FIELD_H
