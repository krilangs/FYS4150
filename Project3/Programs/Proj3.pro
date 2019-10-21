TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        main.cpp \
        weights.cpp

INCLUDEPATH += C:\Qt\5.13.0\mingw73_64\include
DEPENDPATH += C:\Qt\5.13.0\mingw73_64\include

QMAKE_CXXFLAGS += -fopenmp
LIBS += -fopenmp

HEADERS += \
    weights.h
