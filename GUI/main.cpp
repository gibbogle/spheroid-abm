/****************************************************************************
****************************************************************************/

//! [0]
#include <QApplication>

#include "mainwindow.h"
#include "log.h"

LOG_DECLARE;

int main(int argc, char *argv[])
{
	//initialize file logger
	LOG_INIT("zzz.log");

    QApplication app(argc, argv);

    MainWindow mainWin;
    mainWin.show();

    return app.exec();
}
//! [0]
