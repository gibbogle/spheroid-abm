#ifndef FIELD_H
#define FIELD_H

#include <QDialog>
#include <QMessageBox>
#include <QtGui>
//#include <QTranslator>

struct field_data {
    int site[3];
    int state;
    double volume;
    double conc[4];
};

typedef field_data FIELD_DATA;

#define X_AXIS 0
#define Y_AXIS 1
#define Z_AXIS 2

#define OXYGEN 0
#define GLUCOSE 1
#define DRUG_A 2
#define DRUG_B 3

extern "C" {
    void get_fieldinfo(int *, int *);
    void get_fielddata(int *, FIELD_DATA *);
}

class Field : public QMainWindow
{
public:
    Field();
    ~Field();
    void chooseParameters();
    void displayField(QWidget *);

    int axis;
    double fraction;
    int constituent;
    bool slice_changed;
    FIELD_DATA *data;
    char msg[1024];

public slots:
    void setConstituent(QAbstractButton* button);
    void setPlane(QAbstractButton* button);
};

#endif // FIELD_H
