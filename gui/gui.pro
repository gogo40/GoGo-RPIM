#-------------------------------------------------
#
# Project created by QtCreator 2012-06-28T23:05:50
#
#-------------------------------------------------

QT       += core gui

TARGET = gui
TEMPLATE = app


SOURCES += main.cpp\
mainwindow.cpp \
genDom.cpp \
runSimulation.cpp \
maxwellT.cpp \
../domains/supportdomain.cc \
../domains/domains.cc \
../rpim/rpim.cc \
../lane/lane.cc \
../math/point.cc \
../math/matrix.cc \
../maxwell/maxwell.cc \
    genVideo.cpp \
    saveResVideo.cpp \
    savevideoform.cpp

HEADERS  += mainwindow.h \
genDom.h \
runSimulation.h \
../domains/supportdomain.h \
../domains/domains.h \
../rpim/rpim.h \
../lane/lane.h \
../math/point.h \
../math/matrix.h \
../maxwell/maxwell.h \
    genVideo.h \
    saveResVideo.h \
    savevideoform.h

FORMS    += mainwindow.ui \
    savevideoform.ui


INCLUDEPATH += ../
DEPENDPATH += ../
