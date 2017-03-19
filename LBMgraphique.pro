#-------------------------------------------------
#
# Project created by QtCreator 2017-02-24T18:50:04
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = LBMgraphique
TEMPLATE = app


SOURCES += main.cpp\
        grid.cpp

HEADERS  += grid.h \
    lbm.hpp \
    Debug.hpp

FORMS    += grid.ui

CONFIG += c++11

LIBS += -lboost_system
