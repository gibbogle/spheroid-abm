#-------------------------------------------------
#
# Project created by QtCreator 2011-04-18T11:11:25
#
#-------------------------------------------------

QT       += core gui

LFLAGS += -staticlib

TARGET = cell_tiffer
TEMPLATE = app

INCLUDEPATH  += "C:/Program Files (x86)/ITK/include/itk-4.7"

SOURCES += main.cpp\
           mainwindow.cpp\
		   tiff.cpp\
		   cell.cpp

HEADERS  += mainwindow.h\
            location.h

FORMS    += mainwindow.ui

win32: LIBS += -L"C:/Program Files (x86)/ITK/lib" -lITKIOTiff-4.7 -litktiff-4.7 -litksys-4.7 -litkzlib-4.7 \
        -lITKIOImagebase-4.7 -lITKPath-4.7 -lITKIOJPEG-4.7 -litkjpeg-4.7 -litkvnl-4.7 -litkvnl_algo-4.7 \
        -litkNetlibSlatec-4.7 -litkv3p_netlib-4.7

LIBS += "C:/Program Files (x86)/ITK/lib/ITKCommon-4.7.lib"

win32: LIBS += advapi32.lib shell32.lib

#DESTDIR = ../../bin
