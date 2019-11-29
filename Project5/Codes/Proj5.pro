TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        SolarSystem.cpp \
        classes.cpp \
        main.cpp \
        utils.cpp

HEADERS += \
    classes.h \
    planetsLib.h \
    utils.h
