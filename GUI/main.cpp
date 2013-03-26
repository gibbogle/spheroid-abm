/****************************************************************************
****************************************************************************/

//! [0]
#include <QApplication>

#include "mainwindow.h"
#include "log.h"

LOG_DECLARE;

int main(int argc, char *argv[])
{
    int res;
    char msg[128];
	//initialize file logger
	LOG_INIT("zzz.log");

    QApplication app(argc, argv);

    MainWindow mainWin;
    mainWin.show();

    LOG_MSG("A");
    res = app.exec();
    sprintf(msg,"Result code: %d",res);
    LOG_MSG(msg);
    return res;
}
//! [0]
