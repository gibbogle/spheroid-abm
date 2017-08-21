#ifndef QMYGROUPBOX_H
#define QMYGROUPBOX_H

#include <qgroupbox.h>
#include <QMouseEvent>

#include "log.h"
LOG_USE();

class QMyGroupBox: public QGroupBox
{
    Q_OBJECT

public:
    QMyGroupBox(QWidget *parent = 0);

signals:
    void groupBoxClicked(QString text);

private:

    void mousePressEvent (QMouseEvent *event);

};

#endif

